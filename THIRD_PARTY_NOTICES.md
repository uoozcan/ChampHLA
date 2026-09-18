# Third-party notices

PanelHLA is licensed under the MIT License (see `LICENSE`). That licence covers this
repository's own code and nothing else.

PanelHLA does not redistribute any HLA caller. Each is obtained and executed
separately — as a Singularity container pulled at run time, or in SpecHLA's case
from a local checkout — so each remains under its own terms. Because the callers run
as separate processes rather than being linked, a copyleft licence on a caller does
not extend to this repository.

**This matters for what you may do with results, not only with code.** Several
callers restrict commercial use. An assembled pipeline is constrained by the most
restrictive component it invokes, so in practice a full PanelHLA panel is for
research use even though PanelHLA itself is MIT. Check the terms of every caller you
enable before using it commercially.

This file is a description of licence terms, not legal advice.

## Verified from the deployed installation

Each entry below was read from the licence file shipped with the deployed tool on
CSC Roihu, at the path given.

| Caller | Licence | Evidence |
|---|---|---|
| OptiType | BSD 3-Clause | `/projappl/project_2008084/OptiType/LICENSE` |
| Kourami | BSD 3-Clause | `hla_tools/kourami/kourami-0.9.6/LICENSE` |
| arcasHLA | GPL-3.0 | `/projappl/project_2008084/arcas/.../LICENSE` |
| SpecHLA | MIT | `/projappl/project_2008084/SpecHLAx/License` |

SpecHLA bundles two components with their own terms, both permissive:

| Component | Licence | Evidence |
|---|---|---|
| HapCUT2 | BSD 2-Clause | `hla_tools/HapCUT2-1.3.4/LICENSE` |
| SpecHap | MIT | `hla_tools/SpecHap-1.0.1/LICENSE.md` |

## TO CONFIRM before submission

These callers are deployed only as container images and ship no licence file that
could be read from the installation, so their terms could not be verified here. They
must be confirmed from upstream and recorded before release. Two of them are widely
understood to restrict commercial use, which is precisely why they must be checked
rather than assumed.

| Caller | Where to confirm | Note |
|---|---|---|
| HLA-HD | Kyoto University distribution terms | Built from a separately obtained source copy (`HLA_collections/hlahd.dockerfile`); a licence agreement is understood to be required, and academic use free. Confirm. |
| POLYSOLVER | Broad Institute distribution terms | Understood to restrict use to non-commercial research. Confirm. |
| T1K | upstream repository | No licence file in the deployment. |
| seq2HLA | upstream repository | No licence file in the deployment. Not a member of any voting panel in this study. |

## Reference data

IPD-IMGT/HLA allele databases are redistributed inside several caller images under
their own terms. The IPD-IMGT/HLA database is released under a Creative Commons
Attribution–NoDerivs licence. Callers bundle releases spanning 3.11.0 to 3.63.0; the
per-caller resolution is recorded in
`provenance/external_evidence/db_release_resolution.json` and tabulated in the
manuscript Methods.

Human sequence data used here — 1000 Genomes, Geuvadis, NCI-60, HPRC — are public and
consented, and are subject to the data-use terms of their respective consortia.
