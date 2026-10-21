"""Agents trialled in the *other* genome-instability syndromes - mechanism-class transfer.

Why this exists
---------------
For an ultra-rare constitutional chromosomal-instability syndrome (the index
case: mosaic variegated aneuploidy from biallelic BUB1B loss) the usual
repurposing move - enumerate targets reachable from the causal gene's pathway
and look for drugs against them - returns nothing usable, because the mitotic
checkpoint and its neighbours are pan-essential: inhibiting them is
chemotherapy, not therapy. That search exhausts the *target* axis and
concludes "no approved drug".

There is a second axis. Fanconi anaemia, ataxia-telangiectasia, Bloom, Werner,
Nijmegen breakage, dyskeratosis congenita, Rothmund-Thomson and the rest share
the *therapeutic problem* - constitutional genome instability with growth
failure, marrow failure and cancer predisposition - while differing in the
molecular lesion. What is actually given to those children is a real
hypothesis for this child, and it arrives with paediatric safety and dosing
data already attached, which de-novo target hypotheses never do.

Not phenotype similarity
------------------------
This is deliberately *not* an HPO phenotype-overlap search. That was tested on
the index patient and returned noise (pemphigoid gestationis,
Legg-Calve-Perthes) because the patient's HPO annotation is generic
(microcephaly, growth delay, intellectual disability). The sibling set here is
defined by *mechanism class*: the cell-biological failure mode, curated and
listed explicitly below.

How the sibling set was chosen (provenance)
-------------------------------------------
A hand-curated, explicitly enumerated seed set (`SEED`), one row per disorder,
each carrying its Orphanet ID, the mechanism subclass, and the MeSH heading
that ChEMBL's `drug_indication.mesh_heading` uses (None where ChEMBL has no
such heading - those disorders are kept, with zero agents, because "nothing
has been tried" is itself an answer). Inclusion rule, applied by hand:

  core      a constitutional, monogenic disorder whose primary cell-biological
            defect is maintenance of genome or chromosome integrity - DNA
            repair, replication-fork/helicase function, DNA-damage signalling,
            telomere maintenance, sister-chromatid cohesion, or mitotic
            chromosome segregation (the index disease's own subclass).
  adjacent  shares the clinical problem (bone-marrow failure plus cancer
            predisposition, or a progeroid course) but the primary lesion is a
            different mechanism class - ribosome biogenesis, nuclear lamina.
            Kept, tiered, and excluded by default: useful as a contrast, not
            as a mechanism-class transfer.

Deliberately excluded: 'Anemia, Aplastic' (41 agents) and 'Bone Marrow Failure
Disorders' (5). They are not constitutional genome-instability syndromes, they
are largely the acquired/consequence pool, and including them would swamp the
table with transplant-conditioning and marrow-support agents - exactly the
class this tool exists to separate out.

Orphanet widening was attempted and rejected. The local Orphanet parse
(`hierarchy.tsv`) gives parent nodes as bare IDs of type `Category` with no
names available in the parsed TSVs, so a widened set would not be auditable.
The nearest category shared by MVA, Fanconi and Nijmegen is 183422 (47
children), whose membership - Noonan, Peutz-Jeghers, Costello, Denys-Drash -
shows it is "genetic cancer-predisposing syndrome", not genome instability.
Orphanet is therefore used only for what it does cleanly here: attaching
causal genes to each curated disorder from `genes.tsv`, so the tool can be
entered by gene as well as by disease name.

The distinction that matters in the output
------------------------------------------
A phase-ranked list of what was given to Fanconi patients is dominated by
busulfan, fludarabine, cyclophosphamide, alemtuzumab, methotrexate,
ciclosporin - transplant conditioning and GvHD prophylaxis. Those treat the
*consequence* (marrow failure) and transfer nothing except the observation
that transplant is possible. The transferable agents are the handful aimed at
the instability itself: metformin and quercetin (aldehyde/oxidative DNA damage
in FA), niacinamide (NAD+ for PARP/sirtuin-dependent repair in A-T),
amifostine (radical scavenging). Every row is therefore labelled with a
curated `intent`, keyed on (disease, drug) because the same molecule can be
given for different reasons in different syndromes - sirolimus is GvHD
prophylaxis in Fanconi and progerin clearance in progeria.

Usage:
    python -I src/build_sibling_diseases.py \
        --out ../mva_llm/data/kg/sibling_diseases.parquet
"""
from pathlib import Path
import argparse
import os
import sqlite3

import pandas as pd

CHEMBL = str(Path(__file__).resolve().parents[1] / "data/chembl37_db/chembl_37/chembl_37_sqlite/chembl_37.db")
ORPHANET = str(Path(__file__).resolve().parents[1] / "data/orphanet")
OUT = str(Path(__file__).resolve().parents[2] / "mva_llm/data/kg/sibling_diseases.parquet")

# disease, orpha, mesh_heading (None = ChEMBL has no such indication heading),
# tier, mechanism subclass, curated key genes
SEED = [
    # --- the index disease itself: recorded, flagged, never returned as a sibling
    ("Mosaic variegated aneuploidy syndrome", 1052, None, "index",
     "mitotic chromosome segregation (spindle assembly checkpoint)",
     "BUB1B BUB1 BUB3 CEP57 TRIP13"),

    # --- core: constitutional genome/chromosome instability
    ("Fanconi anemia", 84, "Fanconi Anemia", "core",
     "interstrand-crosslink repair (FA-BRCA pathway)", "FANCA FANCC FANCD2 BRCA2"),
    ("Ataxia-telangiectasia", 100, "Ataxia Telangiectasia", "core",
     "double-strand-break signalling (ATM)", "ATM"),
    ("Ataxia-telangiectasia-like disorder", 251347, None, "core",
     "double-strand-break sensing (MRN complex)", "MRE11"),
    ("Nijmegen breakage syndrome", 647, "Nijmegen Breakage Syndrome", "core",
     "double-strand-break sensing (MRN complex)", "NBN"),
    ("Nijmegen breakage syndrome-like disorder", 240760, None, "core",
     "double-strand-break sensing (MRN complex)", "MRE11 RAD50"),
    ("Bloom syndrome", 125, "Bloom Syndrome", "core",
     "RecQ helicase, replication-fork resolution", "BLM"),
    ("Werner syndrome", 902, "Werner Syndrome", "core",
     "RecQ helicase, replication/telomere maintenance (progeroid)", "WRN"),
    ("Rothmund-Thomson syndrome", 221016, "Rothmund-Thomson Syndrome", "core",
     "RecQ helicase (RECQL4); type 1 is ANAPC1, an APC/C subunit", "RECQL4 ANAPC1"),
    ("Dyskeratosis congenita", 1775, "Dyskeratosis Congenita", "core",
     "telomere maintenance", "DKC1 TERT TERC TINF2 RTEL1"),
    ("Xeroderma pigmentosum", 910, "Xeroderma Pigmentosum", "core",
     "nucleotide-excision repair", "XPA XPC ERCC2 ERCC3 ERCC4 ERCC5 DDB2"),
    ("Cockayne syndrome", 90321, None, "core",
     "transcription-coupled nucleotide-excision repair", "ERCC6 ERCC8"),
    ("Trichothiodystrophy", 33364, None, "core",
     "nucleotide-excision repair / TFIIH", "ERCC2 ERCC3 GTF2H5 MPLKIP"),
    ("LIG4 syndrome", 99812, None, "core",
     "non-homologous end joining", "LIG4"),
    ("Cernunnos-XLF deficiency", 169079, None, "core",
     "non-homologous end joining", "NHEJ1"),
    ("Seckel syndrome", 808, None, "core",
     "replication-stress response (ATR) / centrosome integrity",
     "ATR ATRIP RBBP8 CEP152 PLK4"),
    ("Warsaw breakage syndrome", 280558, None, "core",
     "sister-chromatid cohesion (cohesinopathy)", "DDX11"),
    ("Roberts syndrome", 3103, None, "core",
     "sister-chromatid cohesion; premature centromere separation", "ESCO2"),
    ("ICF syndrome", 2268, None, "core",
     "centromeric heterochromatin instability (DNA methylation)",
     "DNMT3B ZBTB24 CDCA7 HELLS"),
    ("Li-Fraumeni syndrome", 524, "Li-Fraumeni Syndrome", "core",
     "DNA-damage response / apoptosis (TP53)", "TP53"),

    # --- adjacent: shared clinical problem, different mechanism class
    ("Hutchinson-Gilford progeria syndrome", 740, "Progeria", "adjacent",
     "nuclear lamina (progerin); genome instability is secondary", "LMNA ZMPSTE24"),
    ("Shwachman-Diamond syndrome", 811, "Shwachman-Diamond Syndrome", "adjacent",
     "ribosome maturation (ribosomopathy); inherited marrow failure", "SBDS"),
    ("Diamond-Blackfan anemia", 124, "Anemia, Diamond-Blackfan", "adjacent",
     "ribosome biogenesis (ribosomopathy); inherited marrow failure",
     "RPS19 RPL5 RPL11"),
]

# Curated intent. Keyed on drug name, or on ("MeSH heading", "DRUG") where the
# same molecule is given for a different reason in a different syndrome.
# (intent, rationale)
_CONDITIONING = "consequence_directed"
INTENT = {
    # --- aimed at the genome instability / DNA-damage load itself
    "METFORMIN": ("instability_directed",
        "Complex-I inhibitor trialled to lower endogenous aldehyde and "
        "oxidative DNA-damage load and slow clonal evolution; the only "
        "repurposing-style candidate with paediatric Fanconi dosing data."),
    "METFORMIN HYDROCHLORIDE": ("instability_directed",
        "Same hypothesis as metformin: reduce the endogenous DNA-damage load "
        "on haematopoietic stem cells rather than treat marrow failure after "
        "the fact. Phase 2 in Fanconi anaemia, in children."),
    "QUERCETIN": ("instability_directed",
        "Oral flavonoid antioxidant, trialled in Fanconi anaemia to reduce "
        "reactive-oxygen-driven DNA damage in marrow progenitors."),
    "NIACINAMIDE": ("instability_directed",
        "NAD+ precursor; the rationale is substrate supply for PARP- and "
        "sirtuin-dependent DNA repair, i.e. the damage response itself."),
    "AMIFOSTINE": ("instability_directed",
        "Free-radical scavenger / cytoprotectant, given to lower the DNA "
        "lesion burden rather than to treat its consequences."),
    "AFAMELANOTIDE": ("instability_directed",
        "MC1R agonist increasing eumelanin: reduces the UV-induced lesion "
        "load upstream of the repair defect. Note this only transfers where "
        "the damage source is exogenous UV - it has no bearing on a mitotic "
        "mis-segregation defect."),
    "LONAFARNIB": ("other_mechanism_directed",
        "Farnesyltransferase inhibitor blocking progerin farnesylation - "
        "aimed squarely at the progeria lesion. Approved, so the paediatric "
        "safety data are real, but the mechanism is lamin processing, not "
        "genome instability."),
    "PROGERININ": ("other_mechanism_directed",
        "Blocks the progerin-lamin A interaction; lesion-directed for "
        "progeroid laminopathy (and trialled in Werner syndrome), not a "
        "general genome-instability agent."),
    "PRAVASTATIN": ("other_mechanism_directed",
        "Given with a bisphosphonate in progeria to block protein "
        "prenylation upstream of progerin; lamina-directed, not DNA-repair "
        "directed."),
    "ZOLEDRONIC ACID": ("other_mechanism_directed",
        "Farnesyl-diphosphate-synthase inhibitor, used in the progeria "
        "prenylation-blockade combination (and for the bone disease)."),
    "LEUCINE": ("other_mechanism_directed",
        "Translation substrate trialled in Diamond-Blackfan anaemia; aimed "
        "at the ribosomopathy lesion, a different mechanism class."),
    "TRIFLUOPERAZINE": ("other_mechanism_directed",
        "Repurposed in Diamond-Blackfan anaemia for an erythropoiesis effect "
        "unrelated to genome maintenance."),
    "BITOPERTIN": ("other_mechanism_directed",
        "GlyT1 inhibitor modulating erythroid haem synthesis; "
        "ribosomopathy/erythropoiesis directed."),
    "SOTATERCEPT": ("other_mechanism_directed",
        "Activin-ligand trap acting on late erythroid maturation."),

    # --- transplant conditioning and stem-cell handling
    "BUSULFAN": (_CONDITIONING, "Myeloablative conditioning for haematopoietic "
        "stem-cell transplant. Treats the consequence (marrow failure) and is "
        "itself an alkylator - in a chromosome-breakage syndrome it is given "
        "at reduced intensity precisely because these patients cannot "
        "tolerate DNA-crosslinking agents."),
    "TREOSULFAN": (_CONDITIONING, "Reduced-toxicity alkylator conditioning "
        "before transplant."),
    "MELPHALAN": (_CONDITIONING, "Alkylator conditioning before transplant."),
    "FLUDARABINE": (_CONDITIONING, "Lymphodepleting conditioning before "
        "transplant."),
    "FLUDARABINE PHOSPHATE": (_CONDITIONING, "Lymphodepleting conditioning "
        "before transplant."),
    "CYCLOPHOSPHAMIDE": (_CONDITIONING, "Conditioning / post-transplant "
        "GvHD prophylaxis; a DNA-crosslinking agent, markedly more toxic in "
        "crosslink-repair-deficient patients."),
    "ALEMTUZUMAB": (_CONDITIONING, "Anti-CD52 serotherapy for transplant "
        "conditioning and graft tolerance."),
    "BRIQUILIMAB": (_CONDITIONING, "Anti-cKIT antibody for antibody-based, "
        "genotoxic-free conditioning - interesting for this patient group but "
        "still transplant conditioning."),
    "PLERIXAFOR": (_CONDITIONING, "CXCR4 mobiliser used to harvest "
        "haematopoietic stem cells for transplant or gene therapy."),
    "ANCESTIM": (_CONDITIONING, "Stem-cell-factor analogue used for "
        "mobilisation."),

    # --- GvHD prophylaxis / immunosuppression, i.e. transplant aftercare
    "METHOTREXATE": ("consequence_directed", "GvHD prophylaxis after "
        "transplant (an antifolate; not a disease-modifying agent here)."),
    "MYCOPHENOLATE MOFETIL": ("consequence_directed", "GvHD prophylaxis."),
    "CYCLOSPORINE": ("consequence_directed", "GvHD prophylaxis / "
        "immunosuppression for marrow failure."),
    "TACROLIMUS ANHYDROUS": ("consequence_directed", "GvHD prophylaxis."),
    "TACROLIMUS": ("consequence_directed", "GvHD prophylaxis."),
    "METHYLPREDNISOLONE": ("consequence_directed", "Steroid for GvHD / "
        "immune-mediated marrow failure."),
    "PREDNISONE": ("consequence_directed", "Steroid for GvHD or (in "
        "Diamond-Blackfan) to drive erythropoiesis; symptomatic."),
    "DACLIZUMAB": ("consequence_directed", "Anti-IL2R immunosuppression."),
    "ETANERCEPT": ("consequence_directed", "TNF blockade for GvHD."),
    "RITUXIMAB": ("consequence_directed", "Anti-CD20, for post-transplant "
        "lymphoproliferative disease and autoimmune cytopenias."),

    # --- marrow support and other downstream consequences
    "ELTROMBOPAG": ("consequence_directed", "Thrombopoietin-receptor agonist "
        "stimulating residual stem cells; supports the failing marrow, does "
        "not touch the instability."),
    "FILGRASTIM": ("consequence_directed", "G-CSF for neutropenia; note the "
        "standing concern about clonal evolution in marrow-failure "
        "syndromes."),
    "DANAZOL": ("consequence_directed", "Attenuated androgen; raises blood "
        "counts, and in telomere disease lengthens telomeres via oestrogen-"
        "receptor-mediated TERT induction - arguably the one marrow-support "
        "agent with a mechanistic claim, but still aimed at the cytopenia."),
    "OXANDROLONE": ("consequence_directed", "Anabolic androgen for marrow "
        "support and growth."),
    "DEFERASIROX": ("consequence_directed", "Iron chelation for transfusional "
        "overload - a consequence of a consequence."),

    # --- symptomatic
    "DEXAMETHASONE SODIUM PHOSPHATE": ("symptomatic", "Low-dose intra-"
        "erythrocyte dexamethasone trialled in ataxia-telangiectasia for "
        "neurological function; the mechanism is not DNA repair and the "
        "effect is symptomatic."),
    "TRIHEPTANOIN": ("symptomatic", "Anaplerotic medium-chain triglyceride; "
        "metabolic support in ataxia-telangiectasia."),
    "BACLOFEN": ("symptomatic", "GABA-B agonist for spasticity."),

    # --- gene-specific
    "MOZAFANCOGENE AUTOTEMCEL": ("gene_specific_correction",
        "Lentiviral FANCA gene therapy of autologous stem cells. Corrects one "
        "specific gene - it is the strongest Fanconi result in this table and "
        "the least transferable, since nothing about it carries to BUB1B."),

    # --- indication mappings that look spurious
    "TRENONACOG ALFA": ("unclassified",
        "Recombinant factor IX. There is no plausible ataxia-telangiectasia "
        "rationale; treat this as a ChEMBL indication-mapping artefact rather "
        "than evidence."),
    # --- context-specific overrides
    ("Li-Fraumeni Syndrome", "METFORMIN"): ("instability_directed",
        "Trialled as cancer chemoprevention in germline TP53 carriers, on a "
        "metabolic / DNA-damage-load rationale rather than to treat an "
        "established tumour."),
    ("Progeria", "EVEROLIMUS"): ("other_mechanism_directed",
        "mTOR inhibition to clear progerin by autophagy - lesion-directed in "
        "progeria, unlike its use as GvHD prophylaxis elsewhere."),
    ("Fanconi Anemia", "SIROLIMUS"): ("consequence_directed",
        "mTOR inhibition used here as GvHD prophylaxis after transplant, not "
        "for the mTOR/genome-instability rationale it carries in progeria."),
}

AGENT_SQL = """
SELECT md.pref_name, md.chembl_id, MAX(di.max_phase_for_ind) AS max_phase_for_ind,
       md.max_phase AS molecule_max_phase, md.first_approval, md.oral,
       md.parenteral, md.black_box_warning, md.withdrawn_flag, md.molecule_type,
       (SELECT GROUP_CONCAT(DISTINCT dm.mechanism_of_action)
          FROM drug_mechanism dm WHERE dm.molregno = md.molregno) AS moa,
       (SELECT GROUP_CONCAT(DISTINCT ac.level2_description)
          FROM molecule_atc_classification mac
          JOIN atc_classification ac ON ac.level5 = mac.level5
         WHERE mac.molregno = md.molregno) AS atc_class,
       (SELECT GROUP_CONCAT(DISTINCT mac.level5)
          FROM molecule_atc_classification mac
         WHERE mac.molregno = md.molregno) AS atc_codes,
       (SELECT COUNT(*) FROM drug_warning dw WHERE dw.molregno = md.molregno) AS n_warnings
  FROM drug_indication di
  JOIN molecule_dictionary md ON md.molregno = di.molregno
 WHERE di.mesh_heading = ?
 GROUP BY md.molregno
 ORDER BY max_phase_for_ind DESC, md.pref_name
"""


def orphanet_genes(directory=ORPHANET):
    """orpha -> set of causal gene symbols, from the local Orphanet parse."""
    g = pd.read_csv(os.path.join(directory, "genes.tsv"), sep="\t", dtype=str)
    return g.groupby("orpha")["symbol"].apply(lambda s: set(s.dropna())).to_dict()


def intent_for(mesh, drug):
    hit = INTENT.get((mesh, drug)) or INTENT.get(drug)
    if hit is None:
        return "unclassified", ("Not reviewed by hand. Unclassified means no "
                                "curated call was made, not that the mechanism "
                                "is unknown.")
    return hit


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chembl", default=CHEMBL)
    p.add_argument("--orphanet", default=ORPHANET)
    p.add_argument("--out", default=OUT)
    args = p.parse_args()

    og = orphanet_genes(args.orphanet)
    con = sqlite3.connect(f"file:{args.chembl}?mode=ro", uri=True)
    rows = []
    for disease, orpha, mesh, tier, mech, curated_genes in SEED:
        genes = sorted(set(curated_genes.split()) | og.get(str(orpha), set()))
        agents = con.execute(AGENT_SQL, (mesh,)).fetchall() if mesh else []
        base = {"disease": disease, "orpha": f"ORPHA:{orpha}",
                "mesh_heading": mesh, "tier": tier, "mechanism_subclass": mech,
                "genes": " ".join(genes), "n_agents": len(agents),
                "in_chembl": mesh is not None}
        if not agents:
            rows.append({**base, "drug": None, "chembl_id": None,
                         "max_phase_for_ind": None, "molecule_max_phase": None,
                         "first_approval": None, "oral": None,
                         "parenteral": None, "black_box_warning": None,
                         "withdrawn": None, "molecule_type": None, "moa": None,
                         "atc_class": None, "atc_codes": None, "n_warnings": None,
                         "intent": None, "intent_rationale": None})
            continue
        for (name, cid, ph, mph, appr, oral, par, bbw, wd, mtype, moa, atc,
             atcc, nw) in agents:
            intent, why = intent_for(mesh, (name or "").upper())
            rows.append({**base, "drug": name, "chembl_id": cid,
                         "max_phase_for_ind": None if ph is None else float(ph),
                         "molecule_max_phase": None if mph is None else float(mph),
                         "first_approval": None if appr is None else int(appr),
                         "oral": bool(oral), "parenteral": bool(par),
                         "black_box_warning": bool(bbw),
                         "withdrawn": bool(wd), "molecule_type": mtype,
                         "moa": moa, "atc_class": atc, "atc_codes": atcc,
                         "n_warnings": int(nw or 0),
                         "intent": intent, "intent_rationale": why})
    con.close()

    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_parquet(args.out, index=False)

    d = out.drop_duplicates("disease")
    print(f"{len(d)} disorders ({(d.tier == 'core').sum()} core, "
          f"{(d.tier == 'adjacent').sum()} adjacent, "
          f"{(d.tier == 'index').sum()} index)")
    print(f"{out.drug.notna().sum()} disease-agent rows, "
          f"{out.drug.dropna().nunique()} distinct agents")
    print("\nagents per disorder:")
    print(d.set_index("disease")[["tier", "n_agents"]].to_string())
    print("\nintent breakdown (core tier):")
    print(out[(out.tier == "core") & out.drug.notna()]
          .intent.value_counts().to_string())
    unc = out[out.intent == "unclassified"].drug.dropna().unique()
    if len(unc):
        print(f"\nunclassified agents ({len(unc)}): {', '.join(sorted(unc))}")
    print(f"\nwrote {args.out} ({len(out):,} rows)")


if __name__ == "__main__":
    main()
