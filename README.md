# ChampHLA publication candidate

This is the isolated, canonical release candidate for the ChampHLA manuscript
recovery and confirmation work. All source projects and the supplied DOCX are
read-only inputs recorded in `ORIGINALS_MANIFEST.json`.

The default publication route is a benchmark/software paper. The conditional
method manuscript is activated only by `decisions/publication_route.json` after
the frozen confirmation gates pass. Historical WGS performance from the
mate-discarding chromosome-6 slice is invalid and cannot enter an accuracy
claim, abstract, or primary figure.

## Scientific boundary

- Research-only HLA-A, HLA-B, and HLA-C typing.
- Exact unordered two-field correctness is primary.
- Only loci with uniquely resolvable laboratory truth enter the locked
  denominator; a subject may contribute one, two, or three A/B/C loci.
- Missing and `no_consensus` calls count as incorrect at the fixed denominator.
- `SimpleTwoThirdsConsensus` is the abstaining primary baseline for the narrow
  claim that Guarded CC resolves consensus failures.
- Always-call plurality, raw CC, MV-floor, best-single-tool, MetaConsensus, and
  individual callers remain visible comparators.
- Subject-unseen 1000G is same-resource confirmation, not independent external
  validation. HPRC is the independent WGS validation arm.

## Canonical checks

```bash
python -m pip install -e '.[test]'
pytest -q
python -m champhla_recovery.cli audit-manuscript-claims \
  --manuscript manuscripts/benchmark/manuscript.md \
  --registry result_registry.tsv \
  --claims manuscripts/claim_audit.tsv \
  --output artifacts/manuscript_audit.json
```

Truth-bearing files must remain outside runtime prediction inputs. The sealed
truth directory is ignored by Git; only its checksum manifest is releasable.
