"""BM25 index over the local abstract corpus, so a mechanism can be *cited*.

Why this exists. The pipeline's other tools read structured databases - ChEMBL,
Open Targets, DepMap, Orphanet. Those record what a drug was developed for and
what a gene does, but not the mechanistic link that makes a candidate from one
disease relevant to another. The v2 run failed on exactly that: the
sibling-disease tool transferred nicotinamide from ataxia-telangiectasia with
the rationale ChEMBL records for it ("NAD+ precursor, substrate for PARP- and
sirtuin-dependent DNA repair"), the direction gate judged that rationale
against a chromosome-missegregation endpoint, and - correctly, on what it was
given - rejected it. The link that would have changed the verdict (SIRT2
deacetylates BubR1 at K668 and stabilises it, so an NAD+ precursor raises BubR1
in a BubR1-hypomorphic animal) is in the literature and in no local database.

So: retrieval over abstracts we already hold, returning real text with a PMID,
a year and a licence. Not a summary and not the model's recall - the point is
that a claim entering the pipeline carries a citation someone can check.

Deliberately BM25 rather than dense embeddings. A local PubMedBERT is on this
machine, but the queries here are keyword-precise ("SIRT2 BubR1 acetylation",
"KIF18A inhibitor"), lexical matching is competitive on those, and - the real
reason - a reviewer can see exactly why a record was returned. That matches how
the rest of this project treats its models.

Known bias, recorded here because it limits every answer the tool gives: the
corpus was assembled from a fixed set of Europe PMC queries about MVA, the
spindle assembly checkpoint and chromosomal instability. It is not a sample of
the literature; absence from it is weak evidence of absence in general.

Usage:
    python -I src/build_literature_index.py --out data/kg/literature.parquet
"""

import argparse
import glob
import json
import os
import re

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Both the original corpus and any later extension; a new directory is picked
# up without touching this file.
CORPUS_GLOBS = ["data/raw/literature/*.jsonl", "data/raw/literature_v2/*.jsonl"]

TAG_RE = re.compile(r"<[^>]+>")


def clean_text(s: str) -> str:
    """Abstracts carry HTML (<h4>Introduction</h4>) and escaped entities."""
    if not s:
        return ""
    s = s.replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&")
    return TAG_RE.sub(" ", s).strip()


def load(root: str) -> pd.DataFrame:
    rows, seen = [], set()
    for pat in CORPUS_GLOBS:
        for path in sorted(glob.glob(os.path.join(root, pat))):
            source = os.path.splitext(os.path.basename(path))[0]
            with open(path) as fh:
                for line in fh:
                    try:
                        d = json.loads(line)
                    except ValueError:
                        continue
                    pmid = d.get("pmid")
                    if not pmid or pmid in seen:
                        continue
                    seen.add(pmid)
                    rows.append({
                        "pmid": str(pmid),
                        "doi": d.get("doi"),
                        "title": clean_text(d.get("title") or ""),
                        "abstract": clean_text(d.get("abstract") or ""),
                        "year": d.get("year"),
                        "journal": d.get("journal"),
                        "licence": d.get("licence"),
                        "corpus_query": source,
                    })
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", default=ROOT)
    p.add_argument("--out", default=os.path.join(ROOT, "data/kg/literature.parquet"))
    args = p.parse_args()

    df = load(args.root)
    if df.empty:
        raise SystemExit("no abstracts found; is the corpus present?")
    df = df[df.abstract.str.len() > 0].reset_index(drop=True)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_parquet(args.out, index=False)

    by = df.corpus_query.value_counts().to_string()
    print(f"{len(df):,} unique abstracts with text -> {args.out}\n")
    print(by, flush=True)
    yrs = pd.to_numeric(df.year, errors="coerce")
    print(f"\nyears {int(yrs.min())}-{int(yrs.max())}, "
          f"{int((yrs >= 2020).sum()):,} from 2020 or later", flush=True)


if __name__ == "__main__":
    main()
