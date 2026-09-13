#!/usr/bin/env python3
"""Derive an audited POLYSOLVER wrapper for a BAM contig convention.

This script intentionally uses Python 3.4-compatible syntax because the pinned
POLYSOLVER container provides an old Python runtime. It changes only the frozen
hg38 interval tokens and the frozen temporary-directory assignment.
"""

from __future__ import print_function

import argparse
import hashlib
import json
import os
import sys


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def transform(source, output, audit, spec_path, contig, tmp_dir):
    with open(spec_path, "r") as handle:
        spec = json.load(handle)
    if spec.get("schema_version") != "champhla-polysolver-wrapper-patch-1":
        raise ValueError("unsupported wrapper-patch specification")
    if contig not in ("6", "chr6"):
        raise ValueError("contig must be exactly 6 or chr6")
    observed_source_hash = file_sha256(source)
    if observed_source_hash != spec.get("source_sha256"):
        raise ValueError("POLYSOLVER wrapper source SHA-256 mismatch")
    with open(source, "rb") as handle:
        text = handle.read().decode("utf-8")

    tmp_token = spec["temporary_directory_token"]
    expected_tmp = int(spec["expected_temporary_directory_substitutions"])
    observed_tmp = text.count(tmp_token)
    if observed_tmp != expected_tmp:
        raise ValueError("temporary-directory token count was %d, expected %d" %
                         (observed_tmp, expected_tmp))
    text = text.replace(tmp_token, "TMP_DIR=" + tmp_dir)

    interval_tokens = spec["hg38_interval_tokens"]
    expected_contig = (int(spec["expected_hg38_contig_substitutions_for_chr6"])
                       if contig == "chr6" else 0)
    observed_contig = 0
    for token in interval_tokens:
        count = text.count(token)
        if count != 1:
            raise ValueError("hg38 interval token %s occurred %d times, expected 1" %
                             (token, count))
        if contig == "chr6":
            text = text.replace(token, "chr" + token, 1)
            observed_contig += 1
    if observed_contig != expected_contig:
        raise ValueError("hg38 contig substitution count was %d, expected %d" %
                         (observed_contig, expected_contig))

    encoded = text.encode("utf-8")
    with open(output, "wb") as handle:
        handle.write(encoded)
    os.chmod(output, 0o755)
    derived_hash = file_sha256(output)
    record = {
        "schema_version": "champhla-polysolver-wrapper-transform-audit-1",
        "transformation_version": spec["transformation_version"],
        "source_sha256": observed_source_hash,
        "derived_wrapper_sha256": derived_hash,
        "contig_convention": contig,
        "temporary_directory": tmp_dir,
        "temporary_directory_substitutions": observed_tmp,
        "hg38_contig_substitutions": observed_contig,
        "expected_hg38_contig_substitutions": expected_contig,
        "spec_sha256": file_sha256(spec_path),
    }
    with open(audit, "w") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return record


def verify(source, output, audit, spec_path, contig, tmp_dir):
    with open(spec_path, "r") as handle:
        spec = json.load(handle)
    with open(audit, "r") as handle:
        record = json.load(handle)
    checks = {
        "schema_version": "champhla-polysolver-wrapper-transform-audit-1",
        "transformation_version": spec["transformation_version"],
        "source_sha256": file_sha256(source),
        "derived_wrapper_sha256": file_sha256(output),
        "contig_convention": contig,
        "temporary_directory": tmp_dir,
        "spec_sha256": file_sha256(spec_path),
    }
    for key, expected in checks.items():
        if record.get(key) != expected:
            raise ValueError("derived-wrapper audit drift for %s" % key)
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--spec", required=True)
    parser.add_argument("--contig", required=True)
    parser.add_argument("--tmp-dir", required=True)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    try:
        if args.verify_only:
            record = verify(args.source, args.output, args.audit, args.spec,
                            args.contig, args.tmp_dir)
        else:
            record = transform(args.source, args.output, args.audit, args.spec,
                               args.contig, args.tmp_dir)
    except (IOError, OSError, ValueError, KeyError) as error:
        print("POLYSOLVER wrapper transformation failed: %s" % error, file=sys.stderr)
        return 2
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
