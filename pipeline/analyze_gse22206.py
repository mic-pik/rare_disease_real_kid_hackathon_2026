#!/usr/bin/env python3
"""Descriptive, donor-aware analysis of the frozen GSE22206 microarray files.

Requires numpy and pandas. No case-control significance tests: one case donor.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd


URLS = {
    "series.soft": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22206&targ=self&form=text&view=full",
    "matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22206/matrix/GSE22206_series_matrix.txt.gz",
    "platform.soft": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6104&targ=self&form=text&view=full",
    "non_normalized.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22206/suppl/GSE22206_non-normalized.txt.gz",
}


def table(text, begin, end):
    return pd.read_csv(io.StringIO(text.split(begin, 1)[1].split(end, 1)[0].strip()), sep="\t")


def contrasts(x):
    """Rows are probes/genes, columns follow verified GEO sample order."""
    result = pd.DataFrame(index=x.index)
    values = x.to_numpy(float)
    for name, cases, controls in [
        ("fibro", [2, 3], [6, 7]),
        ("lcl_control1", [0, 1], [4, 5]),
        ("lcl_control3", [0, 1], [8, 9]),
        ("lcl_relative1", [0, 1], [10, 11]),
        ("lcl_relative2", [0, 1], [12, 13]),
    ]:
        d = values[:, cases, None] - values[:, None, controls]
        result[name + "_delta"] = d.mean(axis=(1, 2))
        result[name + "_pair_min"] = d.min(axis=(1, 2))
        result[name + "_pair_max"] = d.max(axis=(1, 2))
    result["lcl_delta"] = (result.lcl_control1_delta + result.lcl_control3_delta) / 2
    result["lcl_control_difference"] = values[:, [4, 5]].mean(axis=1) - values[:, [8, 9]].mean(axis=1)
    return result


def robust_direction(c):
    """All culture comparisons agree across both tissues and both LCL controls."""
    names = ["fibro", "lcl_control1", "lcl_control3"]
    up = (c[[n + "_pair_min" for n in names]] > 0).all(axis=1)
    down = (c[[n + "_pair_max" for n in names]] < 0).all(axis=1)
    return pd.Series(np.where(up, 1, np.where(down, -1, 0)), index=c.index)


def analyze(source, output, lincs_gene_info=None):
    output.mkdir(parents=True, exist_ok=True)
    text = gzip.open(source / "matrix.txt.gz", "rt").read()
    metadata = {}
    for row in csv.reader(text.splitlines(), delimiter="\t"):
        if row and row[0].startswith("!Sample_"):
            metadata.setdefault(row[0], []).append(row[1:])
    titles = metadata["!Sample_title"][0]
    expected = ["LCL-proband", "Fibro-proband", "LCL-control1", "Fibro-control2", "LCL-control3", "LCL-relative1", "LCL-relative2"]
    assert titles == [f"{s}-{r}" for s in expected for r in [1, 2]], "Unexpected sample layout"
    processing = metadata["!Sample_data_processing"][0]
    assert all("variance-stabilizing" in s and "robust spline" in s for s in processing)
    x = table(text, "!series_matrix_table_begin", "!series_matrix_table_end").set_index("ID_REF")
    assert x.columns.tolist() == metadata["!Sample_geo_accession"][0]
    assert x.shape == (22184, 14) and x.index.is_unique
    assert np.isfinite(x.to_numpy()).all()
    raw = gzip.open(source / "non_normalized.txt.gz", "rt").read()
    raw = pd.read_csv(io.StringIO(raw[raw.index("ID_REF\t"):]), sep="\t").set_index("ID_REF")
    assert raw.index.is_unique and raw.shape[1] == 28
    assert [s.strip() for s in raw.columns[::2]] == metadata["!Sample_description"][0]
    detection = raw.iloc[:, 1::2].reindex(x.index)
    detection.columns = x.columns
    assert np.isfinite(detection.to_numpy()).all()
    assert ((detection >= 0) & (detection <= 1)).all().all()
    ann = table((source / "platform.soft").read_text(), "!platform_table_begin", "!platform_table_end").set_index("ID")
    assert ann.index.is_unique and x.index.isin(ann.index).all()
    ann = ann.reindex(x.index)
    samples = pd.DataFrame({"accession": x.columns, "title": titles,
                            "donor": [s.split("-")[1] for s in titles],
                            "cell_type": [s.split("-")[0] for s in titles],
                            "culture_replicate": [s.split("-")[2] for s in titles],
                            "description": metadata["!Sample_description"][0]})
    samples.to_csv(output / "samples.csv", index=False)
    x.corr().to_csv(output / "sample_correlations.csv")
    qc = x.describe().T
    qc["detected_probes_p_lt_001"] = (detection < .01).sum(axis=0)
    qc.to_csv(output / "sample_qc.csv")
    probes = contrasts(x)
    probes.insert(0, "symbol", ann.Symbol)
    probes.insert(1, "entrez_id", ann.Entrez_Gene_ID)
    # Retain expression in either condition, requiring both cultures in at least
    # one donor for each cell type. This preserves decreased/induced transcripts.
    detected = detection < .01
    fibro_ok = detected.iloc[:, [2, 3]].all(axis=1) | detected.iloc[:, [6, 7]].all(axis=1)
    lcl_ok = pd.concat([detected.iloc[:, pair].all(axis=1) for pair in [[0, 1], [4, 5], [8, 9]]], axis=1).any(axis=1)
    probes["detected_fibro"] = fibro_ok
    probes["detected_lcl"] = lcl_ok
    # Historical platform symbols are retained; no speculative modern remapping.
    unambiguous = ann.Symbol.fillna("").str.fullmatch(r"[A-Za-z0-9_.-]+")
    probes["eligible"] = fibro_ok & lcl_ok & unambiguous
    probes["robust_direction"] = robust_direction(probes)
    probes.to_csv(output / "probe_contrasts.csv")
    kept = probes[probes.eligible]
    genes_x = x.loc[kept.index].groupby(kept.symbol).median()
    genes = contrasts(genes_x)
    genes["n_probes"] = kept.groupby("symbol").size()
    genes["robust_direction"] = robust_direction(genes)
    # Require every eligible probe to support the gene-level direction.
    directions = kept.groupby("symbol").robust_direction.agg(["min", "max"])
    genes["all_probes_agree"] = (directions["min"] == genes.robust_direction) & (directions["max"] == genes.robust_direction) & (genes.robust_direction != 0)
    genes["lcl_exceeds_control_difference"] = genes.lcl_delta.abs() > genes.lcl_control_difference.abs()
    genes["minimum_abs_delta"] = genes[["fibro_delta", "lcl_control1_delta", "lcl_control3_delta"]].abs().min(axis=1)
    genes["exploratory_signature"] = genes.all_probes_agree & genes.lcl_exceeds_control_difference
    genes["relative_direction_supported"] = False
    for sign in [-1, 1]:
        suffix = "pair_min" if sign == 1 else "pair_max"
        supported = pd.Series(True, index=genes.index)
        for relative in ["lcl_relative1", "lcl_relative2"]:
            # Both relatives are sensitivity comparators, not healthy controls.
            probe_supported = kept[relative + "_" + suffix] * sign > 0
            supported &= probe_supported.groupby(kept.symbol).all()
        genes.loc[(genes.robust_direction == sign) & supported, "relative_direction_supported"] = True
    genes.sort_values("minimum_abs_delta", ascending=False).to_csv(output / "gene_contrasts.csv")
    signature = genes[genes.exploratory_signature].sort_values("minimum_abs_delta", ascending=False)
    signature.to_csv(output / "exploratory_signature.csv")
    signature[signature.relative_direction_supported].to_csv(output / "relative_supported_signature.csv")
    for direction, label in [(1, "up"), (-1, "down")]:
        signature[signature.robust_direction == direction].head(100).to_csv(output / f"top100_{label}.csv")
    # Query gene exported even if it fails cross-tissue consistency filters.
    probes[probes.symbol == "BUB1B"].to_csv(output / "bub1b_probe_contrasts.csv")
    x.loc[ann.Symbol == "BUB1B"].to_csv(output / "bub1b_expression.csv")
    detection.loc[ann.Symbol == "BUB1B"].to_csv(output / "bub1b_detection.csv")
    summary = {
        "study": "GSE22206", "arrays": len(samples), "case_donors": 1,
        "control_donors": 3, "relative_donors": 2,
        "probes": len(x), "eligible_probes": len(kept), "eligible_genes": len(genes),
        "signature_genes": len(signature), "signature_up": int((signature.robust_direction == 1).sum()),
        "signature_down": int((signature.robust_direction == -1).sum()),
        "relative_supported_signature_genes": int(signature.relative_direction_supported.sum()),
        "effect_size_sensitivity": {str(t): int((signature.minimum_abs_delta >= t).sum()) for t in [.25, .5, 1]},
        "replicate_correlations": {s: float(x.iloc[:, 2*i].corr(x.iloc[:, 2*i+1])) for i, s in enumerate(expected)},
        "scale": "submitted lumi variance-stabilized / robust-spline-normalized units; differences are NOT log2 fold changes",
        "inference": "descriptive single-case contrasts; no differential-expression p-values or FDR; no drug rescue claim",
        "annotation": "archived GPL6104 symbols; ambiguous symbols excluded; no current-genome probe remapping",
        "signature_rule": "detected in both cultures of >=1 donor per cell type; all eligible probes and all case/control culture pairs agree across both tissues; LCL case difference exceeds control1/control3 difference; relatives excluded from primary contrasts",
    }
    if lincs_gene_info is not None:
        info = pd.read_csv(lincs_gene_info, sep="\t", dtype={"pr_gene_id": str})
        # Match exact stable gene IDs rather than guessing aliases from old symbols.
        ids = kept.groupby("symbol").entrez_id.agg(lambda s: sorted(set(str(int(v)) for v in s.dropna())))
        unique_ids = ids[ids.map(len) == 1].map(lambda s: s[0])
        query = signature.join(unique_ids.rename("entrez_id"), how="inner").reset_index()
        query = query.merge(info, left_on="entrez_id", right_on="pr_gene_id", how="inner", validate="many_to_one")
        # Exclude gene-ID collisions rather than double-weighting one gene.
        query = query[~query.entrez_id.duplicated(keep=False)]
        query.to_csv(output / "lincs_overlap_all.csv", index=False)
        landmarks = query[query.pr_is_lm == 1]
        landmarks.to_csv(output / "lincs_measured_landmark_query.csv", index=False)
        summary["lincs_overlap_genes"] = len(query)
        summary["lincs_measured_landmarks"] = len(landmarks)
        summary["lincs_measured_landmarks_up"] = int((landmarks.robust_direction == 1).sum())
        summary["lincs_measured_landmarks_down"] = int((landmarks.robust_direction == -1).sum())
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    hashes = lambda paths: {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    manifest = {"source_urls": URLS, "inputs_sha256": hashes(source / name for name in URLS),
                "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "outputs_sha256": hashes(p for p in output.iterdir() if p.name != "manifest.json" and p.is_file()),
                "numpy": np.__version__, "pandas": pd.__version__}
    if lincs_gene_info is not None:
        manifest["lincs_gene_info_sha256"] = hashlib.sha256(lincs_gene_info.read_bytes()).hexdigest()
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--lincs-gene-info", type=Path)
    args = parser.parse_args()
    analyze(args.sources, args.output, args.lincs_gene_info)
