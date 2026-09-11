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
