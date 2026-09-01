from pathlib import Path

from champhla_recovery.assays import build_assay_manifest
from champhla_recovery.io import read_tsv


def test_build_assay_manifest(tmp_path: Path):
    wgs = tmp_path / "wgs.index"
    wgs.write_text(
        "##x=y\n#ENA_FILE_PATH\tMD5SUM\tRUN_ID\tSAMPLE_NAME\tPOPULATION\n"
        "ftp://host/HG00001.final.cram\tabc\tERR1\tHG00001\tGBR\n", encoding="utf-8",
    )
    wes = tmp_path / "wes.index"
    wes.write_text(
        "#CRAM\tCRAM_MD5\tCRAI\tCRAI_MD5\n"
        "ftp:/host/data/GBR/HG00001/exome_alignment/a.cram\tdef\tftp:/host/a.crai\tghi\n",
        encoding="utf-8",
    )
    rna = tmp_path / "rna.tsv"
    rna.write_text(
        "Characteristics[individual]\tAssay Name\tFactor Value[ancestry category]\tComment[ENA_RUN]\n"
        "HG00001\tHG00001.lib\tBritish\tERR2\nHG00001\tHG00001.lib\tBritish\tERR2\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.tsv"
    result = build_assay_manifest(str(wgs), str(wes), str(rna), str(output), str(tmp_path / "manifest.json"))
    assert result["counts"] == {"wgs": 1, "wes": 1, "rnaseq": 1}
    rows = read_tsv(output)
    assert len(rows) == 3
    assert all(row["input_url"].startswith("https://") for row in rows)

