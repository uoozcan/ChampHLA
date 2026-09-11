"""The Nextflow workflow must remain parseable by Groovy.

A `\"\"\"` script block is a Groovy GString. Every `$` meant for the shell has to be
escaped, and `$(` is not valid Groovy at all, so the workflow failed to compile on its
first ever execution:

    token recognition error at: '(' @ line 24, column 16.
        ln -s $(realpath ${fastq1}) R1.fastq

A backslash starts a Groovy escape sequence too, so one meant for the shell or a regex
must be doubled. These checks are static stand-ins for compilation, which needs Nextflow
and so cannot run in CI.
"""
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / "workflow"
MODULE_FILES = sorted(WORKFLOW.glob("modules/*.nf"))
NF_FILES = MODULE_FILES + [WORKFLOW / "main.nf"]

# Names Nextflow resolves. Anything else inside ${...} in a script block is a shell
# variable and must be escaped, or Groovy will interpolate it away.
NEXTFLOW_NAMES = {
    "sample_id", "bam", "cram", "fastq1", "fastq2", "reference_fasta",
    "type_flag", "build", "preset", "reference_name", "kourami_db", "kourami_dir",
    "modality", "requiredColumn", "message", "projectDir",
}
VALID_GROOVY_ESCAPES = set("nrtbf\\'\"$u")


def test_workflow_files_exist():
    assert NF_FILES, "no Nextflow files found"


@pytest.mark.parametrize("path", NF_FILES, ids=lambda p: p.name)
def test_no_unescaped_command_substitution(path: Path):
    """`$(` aborts the Groovy lexer outright."""
    offenders = [
        f"{path.name}:{number}: {line.strip()[:90]}"
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1)
        if re.search(r"(?<!\\)\$\(", line)
    ]
    assert not offenders, offenders


@pytest.mark.parametrize("path", MODULE_FILES, ids=lambda p: p.name)
def test_shell_variables_are_escaped(path: Path):
    """A shell variable left unescaped is silently interpolated to nothing by Groovy.

    Only modules are scanned. main.nf is the orchestrator and its `${...}` occurrences are
    ordinary Groovy interpolation in Groovy code, not shell text in a script block.
    """
    offenders = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        for match in re.finditer(r"(?<!\\)\$\{([A-Za-z_][A-Za-z0-9_]*)", line):
            name = match.group(1)
            if name.split(".")[0] in NEXTFLOW_NAMES or name.startswith(("params", "task")):
                continue
            offenders.append(f"{path.name}:{number}: ${{{name}}} in {line.strip()[:70]}")
    assert not offenders, offenders


@pytest.mark.parametrize("path", NF_FILES, ids=lambda p: p.name)
def test_backslashes_are_valid_groovy_escapes(path: Path):
    """A trailing backslash is a line continuation; any other must be doubled."""
    offenders = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        body = line.rstrip()
        if body.endswith("\\"):
            body = body[:-1]
        for match in re.finditer(r"\\(.)", body):
            if match.group(1) not in VALID_GROOVY_ESCAPES:
                offenders.append(f"{path.name}:{number}: \\{match.group(1)} in {body.strip()[:70]}")
    assert not offenders, offenders


def test_workflow_lock_covers_every_workflow_file():
    """Editing a module without re-freezing would leave the lock describing dead code."""
    import json

    lock = json.loads((ROOT / "configs/roihu_workflow_lock.json").read_text(encoding="utf-8"))
    locked = set(lock["files"])
    present = {
        str(p.relative_to(ROOT)).replace("\\", "/")
        for p in WORKFLOW.rglob("*")
        if p.is_file()
    }
    assert locked == present, {"only_in_lock": locked - present, "only_on_disk": present - locked}

def test_locked_files_have_no_carriage_returns():
    """The workflow lock hashes raw bytes, so a CRLF working copy produces a lock that
    fails on the Linux execution platform.

    This happened: the lock was regenerated on Windows, where Python writes CRLF by
    default, and every one of the 14 hashes mismatched on Roihu even though git stored
    the correct LF content. .gitattributes mandates eol=lf; this asserts the working
    copy actually honours it for the files the lock covers.
    """
    offenders = [
        str(path.relative_to(ROOT))
        for path in sorted(WORKFLOW.rglob("*"))
        if path.is_file() and b"\r\n" in path.read_bytes()
    ]
    assert not offenders, offenders

@pytest.mark.parametrize("path", MODULE_FILES, ids=lambda p: p.name)
def test_no_early_exit_readers_in_pipelines(path: Path):
    """Every script block runs under `set -euo pipefail`.

    A reader that exits early closes the pipe, the writer takes SIGPIPE, and pipefail
    turns that into a failed command. This killed EXTRACT_HLA_AND_CONVERT with exit 141
    and SPECHLA_BAM with exit 1 and an empty stderr, on the first run that got far enough
    to execute processes at all.

    The offending readers are `grep -q`, `head -n`, and an `awk` program containing
    `exit`, when any of them is on the right-hand side of a pipe. Reading a file directly
    is fine, which is why `awk 'NR==2{...; exit}' R1.fastq` is not flagged.
    """
    offenders = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if re.search(r"\|\s*grep\s+-[A-Za-z]*q", line):
            offenders.append(f"{path.name}:{number}: grep -q after a pipe")
        if re.search(r"\|\s*head\b", line):
            offenders.append(f"{path.name}:{number}: head after a pipe")
        if re.search(r"\|\s*awk[^|]*\bexit\b", line):
            offenders.append(f"{path.name}:{number}: awk with exit after a pipe")
    assert not offenders, offenders

@pytest.mark.parametrize("path", NF_FILES, ids=lambda p: p.name)
def test_nested_language_escapes_are_doubled(path: Path):
    """Inside a script block a *valid* Groovy escape is as wrong as an invalid one.

    Groovy interprets \\t and \\n and substitutes a real tab or newline before the shell,
    or a nested Python heredoc, ever sees the text. OptiType solved its ILP and then died
    parsing its own result because of this:

        out.write("Gene<TAB>Allele1<TAB>Allele2<NEWLINE>
        SyntaxError: EOL while scanning string literal

    test_backslashes_are_valid_groovy_escapes catches the opposite mistake. Both rules are
    needed: a backslash in a script block must be doubled whichever kind it is.
    """
    offenders = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        for match in re.finditer(r"(?<!\\)\\([ntr])", line):
            offenders.append(
                f"{path.name}:{number}: single-backslash escape in {line.strip()[:70]}"
            )
    assert not offenders, offenders


@pytest.mark.parametrize("path", MODULE_FILES, ids=lambda p: p.name)
def test_region_queries_are_preceded_by_an_index_guard(path: Path):
    """A region query needs a .bai, and CRAM_TO_BAM does not propagate one.

    CRAM_TO_BAM indexes its output but declares only path("${sample_id}.bam"), so Nextflow
    stages the BAM into a consumer's work directory without the index beside it:

        [E::idx_find_and_load] Could not retrieve index file for 'HG00096.bam'
        samtools view: Random alignment retrieval only works for indexed BAM

    EXTRACT_HLA_AND_CONVERT and KOURAMI_BAM already guarded; SPECHLA_BAM and HLAHD_BAM did
    not, and SpecHLA was the last caller failing because of it.
    """
    text = path.read_text(encoding="utf-8")
    does_region_query = re.search(r"samtools view[^\n]*\$\{chr\}:", text)
    if not does_region_query:
        pytest.skip("no region query in this module")
    assert ".bai" in text and "samtools index" in text, (
        f"{path.name} performs a region query without guaranteeing the BAM index"
    )
