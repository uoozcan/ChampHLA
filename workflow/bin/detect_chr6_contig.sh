#!/bin/bash
# Determine the primary chromosome-6 name from a SAM header without relying on
# container-specific awk character-class behaviour.  In particular, mawk 1.3.3
# in the pinned HLA-HD image does not treat a literal tab as [[:space:]].
set -euo pipefail

header=${1:?usage: detect_chr6_contig.sh SAM_HEADER}
test -s "${header}"

has_6=0
has_chr6=0
while IFS=$'\t' read -r -a fields; do
  [[ ${#fields[@]} -gt 0 && "${fields[0]}" == "@SQ" ]] || continue
  for field in "${fields[@]:1}"; do
    case "${field}" in
      SN:6) has_6=1 ;;
      SN:chr6) has_chr6=1 ;;
    esac
  done
done < "${header}"

if [[ $((has_6 + has_chr6)) -ne 1 ]]; then
  printf 'expected exactly one primary chromosome-6 contig (6 or chr6); found 6=%s chr6=%s\n' \
    "${has_6}" "${has_chr6}" >&2
  exit 2
fi

if [[ "${has_chr6}" -eq 1 ]]; then
  printf 'chr6\n'
else
  printf '6\n'
fi
