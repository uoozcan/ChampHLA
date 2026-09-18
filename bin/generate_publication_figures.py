#!/usr/bin/env python3.11
"""Compatibility wrapper for the hub-owned publication figure generator."""

from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parent.parent
HUB_SCRIPTS = REPO_ROOT / "figure_hub" / "scripts"
sys.path.insert(0, str(HUB_SCRIPTS))

from fig_benchmark_main import (  # noqa: E402
    EXPECTED_ACTIVE_STEMS,
    generate_all,
    load_publication_context,
    main,
    parse_args,
    write_publication_captions,
    write_publication_manifest_update,
)


if __name__ == "__main__":
    raise SystemExit(main())
