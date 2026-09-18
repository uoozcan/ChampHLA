#!/bin/bash
# PanelHLA runtime environment on CSC Roihu. Source this; do not execute it.
#
#   source "$(dirname "$0")/../bin/roihu_env.sh"
#
# WHY THIS FILE EXISTS
#
# `module load nextflow` is a Puhti idiom and cannot work on Roihu, for two reasons
# that are each invisible until you hit them:
#
#  1. /etc/profile.d/zz-csc-env.sh -- the script that defines `module` at all --
#     returns immediately when $PS1 is unset, i.e. in every non-interactive shell.
#     A SLURM batch script is non-interactive, so a job that says `module load ...`
#     dies with "module: command not found" one second in. CSC's escape hatch is
#     CSC_ENV_INIT_NON_INTERACTIVE=yes, which must be exported BEFORE sourcing.
#
#  2. Once `module` exists, the bare name `nextflow` has no default version and its
#     modulefile lives behind a prerequisite. Both `bio-apps/v202603` and an
#     explicit version are required.
#
# A third trap, for whoever edits this next: `module` is a shell function, so
# `module load x | head` runs it in a subshell and silently throws away everything
# it set. Never pipe it.
#
# Verified on roihu-cpu-login1, 2026-09-18: nextflow 25.10.2, openjdk 17.0.11,
# samtools 1.21.

export CSC_ENV_INIT_NON_INTERACTIVE=yes

# zz-csc-env.sh tests [ -z "$PS1" ] without guarding it, so under `set -u` -- which
# any careful launcher uses -- sourcing it aborts with "PS1: unbound variable".
# Drop nounset just across the source and put it back exactly as it was.
_panelhla_nounset=0
case "$-" in *u*) _panelhla_nounset=1; set +u ;; esac
# shellcheck disable=SC1091
source /etc/profile.d/zz-csc-env.sh
# An `x && y` here would return non-zero whenever x is false, which aborts a
# caller that uses `set -e` without `set -u`. Use a plain if.
if [ "${_panelhla_nounset}" = 1 ]; then set -u; fi
unset _panelhla_nounset

module load bio-apps/v202603
module load nextflow/25.10.2-standalone
module load samtools/1.21

# Fail loudly here rather than 20 minutes into a run.
_panelhla_missing=""
for _t in nextflow java samtools; do
    command -v "$_t" >/dev/null 2>&1 || _panelhla_missing="${_panelhla_missing} ${_t}"
done
if [ -n "${_panelhla_missing}" ]; then
    echo "roihu_env.sh: not on PATH after module load:${_panelhla_missing}" >&2
    echo "  Check 'module spider nextflow' -- CSC repoints versions between releases." >&2
    return 1 2>/dev/null || exit 1
fi
unset _panelhla_missing _t

# Nextflow writes its framework and plugins here; the default is the home quota,
# which is 15 GB shared with everything else.
export NXF_HOME="${NXF_HOME:-/scratch/project_2008084/.nextflow}"
mkdir -p "${NXF_HOME}"
