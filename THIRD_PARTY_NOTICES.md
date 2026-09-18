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
| POLYSOLVER | **See below — code and image differ** | `polysolver.sif:/home/polysolver/license/POLYSOLVER_License.docx` |

### POLYSOLVER: the code and the image are not under the same terms

POLYSOLVER's own source is distributed under a BSD-style licence (Broad Institute and
Dana-Farber, 2013–2017), which does not restrict commercial use. Its licence then
adds:

> The docker is available for academic/non-profit users only. Commercial use requires
> a Novoalign license or commercial licenses for other software as applicable.

The image bundles GATK, MuTect, Novoalign, eigentools (GPL-3), VCFtools, samtools,
PLINK, PLINK/SEQ, Strelka and BioPerl. **PanelHLA invokes the image, not the source**,
so the academic and non-profit restriction is the one that applies in practice. This
is the single most restrictive term across the panel and is what makes an assembled
PanelHLA run research-use.

Quoting only the BSD header would misdescribe it; so would calling POLYSOLVER
non-commercial software. Both halves have to be stated together.

SpecHLA bundles two components with their own terms, both permissive:

| Component | Licence | Evidence |
|---|---|---|
| HapCUT2 | BSD 2-Clause | `hla_tools/HapCUT2-1.3.4/LICENSE` |
| SpecHap | MIT | `hla_tools/SpecHap-1.0.1/LICENSE.md` |

## TO CONFIRM before submission

These callers ship no licence of their own inside their container image. That was
checked rather than assumed: each image was searched with `singularity exec
--no-mount tmp`, excluding conda package metadata and the host `/tmp`. The result is
recorded so the search is not repeated.

| Caller | In-image search result | Where to confirm |
|---|---|---|
| HLA-HD 1.4.0 | Absent. The two `LICENSE` files in the image (`/usr/local/bin/LICENSE`, `/home/biodocker/bin/LICENSE`) both belong to **bowtie2**, which sits in those directories; `/app/hlahd.1.4.0/Readme.txt` is installation instructions only. | Kyoto University distribution terms. A licence agreement is understood to be required, with academic use free. |
| T1K | Absent. Only conda dependency licences under `/usr/local/conda-meta/`. | Upstream repository. |
| seq2HLA | Absent. Only conda dependency licences. Not a member of any voting panel in this study. | Upstream repository. |

**A licence file inside a container is not automatically the caller's licence.** The
HLA-HD image demonstrates it: a GPL-3 file sits two directories from the tool and
belongs to something else entirely. Attribution here is by reading, not by
proximity.

## Reference data

IPD-IMGT/HLA allele databases are redistributed inside several caller images under
their own terms. The IPD-IMGT/HLA database is released under a Creative Commons
Attribution–NoDerivs licence. Callers bundle releases spanning 3.11.0 to 3.63.0; the
per-caller resolution is recorded in
`provenance/external_evidence/db_release_resolution.json` and tabulated in the
manuscript Methods.

Human sequence data used here — 1000 Genomes, Geuvadis, NCI-60, HPRC — are public and
consented, and are subject to the data-use terms of their respective consortia.
