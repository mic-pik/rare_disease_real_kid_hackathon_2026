"""Normalise the downloaded sources into one knowledge graph of verified edges.

Every training example is generated from these edges, so this is the only place
where facts enter the pipeline. Nothing downstream may assert anything that is
not traceable to a row here.

Outputs (parquet, in --outdir):
    genes.parquet         gene -> symbol, name, function, location, class, constraint
    gene_pathway.parquet  gene -> Reactome pathway
    gene_disease.parquet  gene -> disease, association score, evidence datatypes
    diseases.parquet      disease id -> name, description, therapeutic area
    drugs.parquet         drug -> name, type, max clinical phase, approval
    drug_target.parquet   drug -> target gene, mechanism of action, action type
    drug_disease.parquet  drug -> disease indication, phase
    gene_phenotype.parquet gene -> HPO phenotype (from HPO, not Open Targets)

Usage:
    python -I src/build_kg.py --raw data/raw --outdir data/kg
"""

import argparse
import glob
import os

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


def as_seq(v):
    """Open Targets list columns arrive as numpy arrays, not Python lists;
    an isinstance(v, list) check silently drops every one of them."""
    if isinstance(v, (list, tuple, np.ndarray)):
        return list(v)
    return []


def read_table(d: str) -> pd.DataFrame:
    files = sorted(glob.glob(os.path.join(d, "*.parquet")))
    if not files:
        raise FileNotFoundError(f"no parquet in {d}")
    return pq.ParquetDataset(files).read().to_pandas()


def first_or_none(x):
    """Open Targets stores many fields as lists; take the primary entry."""
    s = as_seq(x)
    return s[0] if s else None


def build_genes(raw: str, out: str) -> pd.DataFrame:
    t = read_table(f"{raw}/opentargets/target")
    g = pd.DataFrame({
        "gene_id": t["id"],
        "symbol": t["approvedSymbol"],
        "name": t["approvedName"],
        "biotype": t["biotype"],
        "function": t["functionDescriptions"].map(
            lambda v: " ".join(as_seq(v)) or None),
        "uniprot": t["proteinIds"].map(
            lambda v: next((p["id"] for p in as_seq(v)
                            if p.get("source") == "uniprot_swissprot"), None)),
        "subcellular": t["subcellularLocations"].map(
            lambda v: "; ".join(sorted({p["location"] for p in as_seq(v)
                                        if p.get("location")})) or None),
        "target_class": t["targetClass"].map(
            lambda v: (first_or_none(v) or {}).get("label")),
        "chromosome": t["genomicLocation"].map(
            lambda v: v.get("chromosome") if isinstance(v, dict) else None),
    })
    # pLI / LOEUF: how intolerant the gene is to loss of function - directly
    # relevant to arguing whether a truncating allele is plausibly causal.
    def constraint(v, kind):
        for c in as_seq(v):
            if c.get("constraintType") == kind:
                return c.get("score")
        return None
    g["pli"] = t["constraint"].map(lambda v: constraint(v, "lof"))

    # Reactome pathways come bundled with the target table.
    rows = []
    for gid, paths in zip(t["id"], t["pathways"]):
        for p in as_seq(paths):
            rows.append((gid, p.get("pathwayId"), p.get("pathway"),
                         p.get("topLevelTerm")))
    gp = pd.DataFrame(rows, columns=["gene_id", "pathway_id", "pathway", "top_level"])

    g.to_parquet(f"{out}/genes.parquet", index=False)
    gp.drop_duplicates().to_parquet(f"{out}/gene_pathway.parquet", index=False)
    print(f"genes {len(g):,}  gene-pathway {len(gp):,}", flush=True)
    return g


def build_diseases(raw: str, out: str) -> pd.DataFrame:
    d = read_table(f"{raw}/opentargets/disease")
    keep = pd.DataFrame({
        "disease_id": d["id"],
        "name": d["name"],
        "description": d.get("description"),
        "therapeutic_areas": (d["therapeuticAreas"].map(lambda v: "; ".join(as_seq(v)) or None)
                              if "therapeuticAreas" in d else None),
    })
    keep.to_parquet(f"{out}/diseases.parquet", index=False)
    print(f"diseases {len(keep):,}", flush=True)
    return keep


def build_associations(raw: str, out: str) -> None:
    a = read_table(f"{raw}/opentargets/association_overall_direct")
    a = a.rename(columns={"targetId": "gene_id", "diseaseId": "disease_id",
                          "score": "assoc_score", "associationScore": "assoc_score",
                          "evidenceCount": "n_evidence"})
    cols = [c for c in ["gene_id", "disease_id", "assoc_score", "n_evidence"] if c in a]
    a = a[cols]

    try:
        dt = read_table(f"{raw}/opentargets/association_by_datatype_direct")
        # Real columns are aggregationType/aggregationValue/associationScore,
        # not datatypeId/score.
        dt = dt.rename(columns={"targetId": "gene_id", "diseaseId": "disease_id",
                                "aggregationValue": "datatype",
                                "associationScore": "datatype_score"})
        dt = dt[["gene_id", "disease_id", "datatype", "datatype_score"]]
        dt.to_parquet(f"{out}/gene_disease_evidence.parquet", index=False)
        print(f"gene-disease evidence rows {len(dt):,}", flush=True)
    except FileNotFoundError:
        print("no datatype table yet - skipping evidence breakdown", flush=True)

    a.to_parquet(f"{out}/gene_disease.parquet", index=False)
    print(f"gene-disease {len(a):,}", flush=True)


def build_drugs(raw: str, out: str) -> None:
    m = read_table(f"{raw}/opentargets/drug_molecule")
    # Real schema: maximumClinicalStage is a string enum ('APPROVAL', 'PHASE_3',
    # ...). There is no isApproved / hasBeenWithdrawn / yearOfFirstApproval
    # column - asking for them with .get() silently produced an all-None table.
    # canonicalSmiles is kept so approved drugs can be scored by the naive Bayes
    # target model.
    drugs = pd.DataFrame({
        "drug_id": m["id"],
        "drug_name": m["name"],
        "drug_type": m.get("drugType"),
        "max_stage": m["maximumClinicalStage"],
        "is_approved": m["maximumClinicalStage"].eq("APPROVAL"),
        "smiles": m.get("canonicalSmiles"),
        "description": m.get("description"),
    })
    drugs.to_parquet(f"{out}/drugs.parquet", index=False)
    print(f"drugs {len(drugs):,} ({int(drugs.is_approved.sum()):,} approved, "
          f"{int(drugs.smiles.notna().sum()):,} with SMILES)", flush=True)

    moa = read_table(f"{raw}/opentargets/drug_mechanism_of_action")
    rows = []
    for _, r in moa.iterrows():
        targets = as_seq(r.get("targets"))
        chembls = as_seq(r.get("chemblIds"))
        for drug in chembls:
            for gene in targets:
                rows.append((drug, gene, r.get("mechanismOfAction"),
                             r.get("actionType"), r.get("targetName")))
    dt = pd.DataFrame(rows, columns=["drug_id", "gene_id", "mechanism",
                                     "action_type", "target_name"]).drop_duplicates()
    dt.to_parquet(f"{out}/drug_target.parquet", index=False)
    print(f"drug-target {len(dt):,}", flush=True)

    try:
        ind = read_table(f"{raw}/opentargets/clinical_indication")
        # Flat table: drugId / diseaseId / maxClinicalStage - there is no nested
        # "indications" list (that was the old Open Targets layout).
        di = ind.rename(columns={"drugId": "drug_id", "diseaseId": "disease_id",
                                 "maxClinicalStage": "max_phase"})
        di = di[["drug_id", "disease_id", "max_phase"]].dropna(subset=["drug_id"])
        di.to_parquet(f"{out}/drug_disease.parquet", index=False)
        print(f"drug-disease {len(di):,}", flush=True)
    except FileNotFoundError:
        print("no clinical_indication table yet", flush=True)


def build_phenotypes(raw: str, out: str) -> None:
    f = f"{raw}/hpo/genes_to_phenotype.txt"
    if not os.path.exists(f):
        print("no HPO file yet", flush=True)
        return
    h = pd.read_csv(f, sep="\t")
    h.columns = [c.lstrip("#").strip() for c in h.columns]
    h.to_parquet(f"{out}/gene_phenotype.parquet", index=False)
    print(f"gene-phenotype {len(h):,}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--only", nargs="*", default=None,
                   help="subset of: genes diseases assoc drugs phenotypes")
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    steps = {
        "genes": lambda: build_genes(args.raw, args.outdir),
        "diseases": lambda: build_diseases(args.raw, args.outdir),
        "assoc": lambda: build_associations(args.raw, args.outdir),
        "drugs": lambda: build_drugs(args.raw, args.outdir),
        "phenotypes": lambda: build_phenotypes(args.raw, args.outdir),
    }
    for name, fn in steps.items():
        if args.only and name not in args.only:
            continue
        print(f"--- {name} ---", flush=True)
        try:
            fn()
        except FileNotFoundError as e:
            print(f"  skipped: {e}", flush=True)


if __name__ == "__main__":
    main()
