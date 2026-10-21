"""Retrieve real abstracts, with citations, for a mechanistic claim.

See src/build_literature_index.py for why this tool exists and what the corpus
is biased towards. The short version: the structured databases record what a
drug was developed for, not the mechanism that would justify moving it to a new
disease, and the v2 run rejected a good candidate for exactly that reason.

The tool returns retrieved text, a PMID and a licence. It deliberately does not
summarise: a claim that enters the pipeline through this tool should be
quotable and checkable, and the stage prompts ask the model to cite rather than
paraphrase. Retrieval is not endorsement - the corpus contains contradictory
findings on the same question, and surfacing both is the point.
"""

from __future__ import annotations

import functools
import os
import re

import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer

KG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                  "data/kg")
INDEX = os.path.join(KG, "literature.parquet")

# Keeps gene symbols (BUB1B, KIF18A) and hyphenated agent codes (AMG-650,
# VLS-1488) as single tokens; the sklearn default splits them.
TOKEN_RE = r"(?u)\b[\w][\w-]+\b"

K1, B = 1.5, 0.75          # standard BM25 constants
TITLE_WEIGHT = 3           # a term in the title is worth more than one in body


@functools.lru_cache(maxsize=1)
def _index():
    df = pd.read_parquet(INDEX)
    text = ((df.title + " ") * TITLE_WEIGHT + df.abstract).str.lower()
    vec = CountVectorizer(token_pattern=TOKEN_RE, min_df=2)
    tf = vec.fit_transform(text)                      # sparse counts
    dl = np.asarray(tf.sum(axis=1)).ravel().astype(np.float64)
    avgdl = dl.mean()
    n = tf.shape[0]
    dfreq = np.asarray((tf > 0).sum(axis=0)).ravel()
    idf = np.log(1.0 + (n - dfreq + 0.5) / (dfreq + 0.5))
    return df, vec, tf.tocsc(), dl, avgdl, idf


def _score(terms: list[str]) -> np.ndarray:
    df, vec, tf, dl, avgdl, idf = _index()
    vocab = vec.vocabulary_
    scores = np.zeros(tf.shape[0], dtype=np.float64)
    norm = K1 * (1 - B + B * dl / avgdl)
    for t in terms:
        j = vocab.get(t)
        if j is None:
            continue
        col = tf[:, j].tocoo()
        f = col.data.astype(np.float64)
        rows = col.row
        scores[rows] += idf[j] * (f * (K1 + 1)) / (f + norm[rows])
    return scores


def _snippet(abstract: str, terms: list[str], width: int = 420) -> str:
    """The span of the abstract that actually matched, not the first N chars."""
    low = abstract.lower()
    hits = [low.find(t) for t in terms if t in low]
    hits = [h for h in hits if h >= 0]
    if not hits:
        return abstract[:width].strip()
    centre = int(np.median(hits))
    start = max(0, centre - width // 2)
    out = abstract[start:start + width].strip()
    return ("..." if start else "") + out + ("..." if start + width < len(abstract) else "")


def literature_evidence(query: str, gene: str | None = None, top: int = 12,
                        min_year: int | None = None) -> dict:
    """Find published abstracts supporting *or contradicting* a claim.

    Use when a candidate's plausibility turns on a mechanism that the
    structured databases do not record - how one protein regulates another's
    stability, what an agent did in a related disease, whether a target's
    vulnerability depends on a pathway this patient has lost. Quote what comes
    back and cite the PMID; do not restate it from memory.
    """
    try:
        df, vec, *_ = _index()
    except FileNotFoundError:
        return {"ok": False, "error": f"no literature index at {INDEX}; "
                                      "run src/build_literature_index.py"}

    q = f"{query} {gene}" if gene else query
    terms = [t.lower() for t in re.findall(TOKEN_RE, q)]
    if not terms:
        return {"ok": False, "error": "empty query"}

    scores = _score(terms)
    known = [t for t in terms if t in vec.vocabulary_]
    if not known or not scores.any():
        return {"ok": True, "query": q, "n_hits": 0, "hits": [],
                "terms_not_in_corpus": [t for t in terms
                                        if t not in vec.vocabulary_],
                "interpretation": (
                    "Nothing in the local corpus matches. That is weak evidence: "
                    "the corpus is a fixed set of queries about MVA, the spindle "
                    "assembly checkpoint and chromosomal instability, not a "
                    "sample of the literature. Report it as 'not found locally', "
                    "never as 'no such evidence exists', and do not substitute "
                    "your own recollection for a retrieved citation."),
                "source": "Europe PMC abstracts held locally"}

    mask = np.ones(len(df), dtype=bool)
    if min_year is not None:
        yr = pd.to_numeric(df.year, errors="coerce")
        mask &= (yr >= min_year).fillna(False).to_numpy()
    scores = np.where(mask, scores, 0.0)

    order = np.argsort(-scores)[:max(1, int(top))]
    hits = []
    for i in order:
        if scores[i] <= 0:
            continue
        r = df.iloc[i]
        hits.append({"pmid": r.pmid, "doi": r.doi, "year": r.year,
                     "journal": r.journal, "licence": r.licence,
                     "title": r.title,
                     "matched_text": _snippet(r.abstract, known),
                     "bm25_score": round(float(scores[i]), 2)})

    return {"ok": True, "query": q, "n_hits": len(hits), "hits": hits,
            "terms_not_in_corpus": [t for t in terms
                                    if t not in vec.vocabulary_],
            "corpus_size": int(len(df)),
            "showing_top_n_of_corpus": int(top),
            "interpretation": (
                "These are retrieved abstracts, ranked by lexical match (BM25). "
                "Ranking is relevance, not quality and not agreement: a highly "
                "ranked abstract may report a negative result or contradict "
                "another hit, and where the corpus disagrees with itself you "
                "should say so rather than take the top hit. Quote the "
                "`matched_text` and cite the PMID for any claim you carry "
                "forward; an abstract is also not a full paper, so treat an "
                "effect size or a dose read from one as approximate. The corpus "
                "is biased by construction towards MVA, the spindle assembly "
                "checkpoint and chromosomal instability. Finally, this is a "
                "ranked slice, not the corpus: a relevant paper whose title "
                "and abstract happen to use different words than your query "
                "can sit just below the cut. Before concluding that something "
                "is not reported, re-query with different phrasing - the "
                "mechanism, the residue, the protein, the phenotype - and "
                "raise `top`."),
            "source": "Europe PMC abstracts held locally; "
                      "see src/build_literature_index.py"}


SCHEMA = {"type": "function", "function": {
    "name": "literature_evidence",
    "description": (
        "Search locally held published abstracts for evidence about a "
        "mechanism, and get back real text with a PMID. Use it when a "
        "candidate's plausibility depends on a link the structured databases "
        "do not record - how one protein regulates another's stability, what "
        "an agent actually did in a related disease, or whether a target's "
        "vulnerability depends on a pathway this patient has lost. Returns "
        "contradicting findings too, which you should report rather than hide. "
        "Cite the PMID; do not paraphrase from memory."),
    "parameters": {"type": "object", "properties": {
        "query": {"type": "string",
                  "description": "Free text, e.g. 'SIRT2 BubR1 acetylation stability'"},
        "gene": {"type": "string", "description": "Optional gene symbol to add"},
        "top": {"type": "integer", "description": "How many abstracts (default 12). Raise it before concluding something is absent."},
        "min_year": {"type": "integer", "description": "Optional earliest year"}},
        "required": ["query"]}}}


if __name__ == "__main__":
    import json
    for q in ["SIRT2 BubR1 acetylation stability lifespan",
              "KIF18A inhibitor chromosomally unstable aneuploid selective",
              "spindle assembly checkpoint competency KIF18A sensitivity"]:
        r = literature_evidence(q, top=3)
        print(f"\n=== {q}  ({r['n_hits']} hits)")
        for h in r["hits"]:
            print(f"  {h['pmid']} {h['year']} [{h['bm25_score']}] {h['title'][:88]}")
