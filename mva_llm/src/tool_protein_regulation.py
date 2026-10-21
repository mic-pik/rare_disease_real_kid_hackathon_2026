"""What is known to control this protein's level — asked systematically.

Why this exists. For a loss-of-function dosage disease the decisive question is
"what raises this protein", and the answer lives in the post-translational
literature. Across three otherwise-identical pipeline runs the model reached
three different conclusions about the same candidate, and the cause was not
reasoning: it was which words it happened to put in a free-text search.

    run 5  "BUB1B BUBR1 acetylation SIRT2 stability"   -> found SIRT2, misread it
    run 6  "SIRT2 BubR1 deacetylation stability"       -> found the K250 paper only
    run 7  "BUB1B BUBR1 protein stability ubiquitination degradation"
                                                        -> dropped "acetylation",
                                                           never reached SIRT2 at all

One term decided whether the mechanism was discovered. That is not a judgement
the model should be making ad hoc, so this tool takes it away: given a gene, it
runs a fixed battery of regulation queries and returns the union. Coverage
becomes a property of the tool, identical on every run, and the model's job
goes back to reading the evidence rather than guessing how to find it.

It is deliberately recall-oriented. Some hits will be irrelevant; that is the
right trade when missing one paper silently flips the conclusion.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import re

import pandas as pd

from tool_literature import INDEX, literature_evidence

# Each axis is a way a protein's abundance or activity gets controlled. The
# list is fixed so that coverage does not vary between runs.
AXES = {
    "acetylation": "{g} {p} acetylation lysine acetyltransferase deacetylase",
    # "abundance lifespan levels" is load-bearing: without it this axis is won
    # by general NAD-metabolism reviews and by BUB1B cancer-expression papers,
    # and the one study reporting an NAD precursor raising this protein in vivo
    # sits at rank 44.
    "sirtuin_nad": "{g} {p} sirtuin SIRT deacetylase NAD abundance lifespan levels",
    "ubiquitination": "{g} {p} ubiquitination E3 ligase proteasome degradation",
    "phosphorylation": "{g} {p} phosphorylation kinase phosphatase stability",
    "protein_level": "{g} {p} protein level abundance half-life stability turnover",
    "transcription": "{g} {p} transcription promoter expression upregulation",
    # The therapeutic question stated directly. Asking what *modifies* a protein
    # is not the same as asking what *raises* it, and for a dosage disease the
    # second is the one that matters; papers reporting an intervention that
    # restored a protein in vivo often describe it in exactly these words.
    "raising_the_protein":
        "{g} {p} increase raise restore abundance levels in vivo treatment",
}



def _alias_from_corpus(gene: str, min_count: int = 5, margin: float = 2.0) -> str:
    """Learn what this protein is also called, from how papers write it.

    The gene symbol and the protein name are frequently different words, and
    the mechanistic literature prefers the protein: 258 abstracts here say
    "BubR1" and never "BUB1B". Querying the symbol alone therefore fails to
    retrieve them at all - no amount of re-ranking helps, because they were
    never in the candidate set.

    Two earlier attempts at this were worse than useless. Requiring the caller
    to pass the alias meant it was simply omitted, and the affected axes
    returned nothing while *looking* like genuine negative results. Deriving it
    statistically, by term enrichment, returned co-expressed genes and GEO
    accession numbers.

    This reads the convention instead: authors write "BUB1B/BUBR1" or
    "BUB1B (BubR1)" the first time they use it. Counting those constructions is
    precise rather than clever, and the result is reported in the output so a
    wrong alias is visible. Requires a clear margin over the runner-up, since
    "FOXM1/BUB1B" is the same construction between two different genes.
    """
    corpus = pd.read_parquet(INDEX)
    text = corpus.title.fillna("") + " " + corpus.abstract.fillna("")
    g = re.escape(gene)
    w = r"([A-Za-z][A-Za-z0-9-]{2,11})"
    pats = [rf"{g}\s*/\s*{w}", rf"{g}\s*\(\s*{w}\s*\)",
            rf"{w}\s*/\s*{g}", rf"{w}\s*\(\s*{g}\s*\)",
            rf"{g}[, ]+also known as\s+{w}"]
    counts: dict[str, int] = {}
    for s_ in text:
        for pat in pats:
            for m in re.findall(pat, s_, flags=re.I):
                k = m.lower()
                if k != gene.lower():
                    counts[k] = counts.get(k, 0) + 1
    if not counts:
        return ""
    ranked = sorted(counts.items(), key=lambda kv: -kv[1])
    top, n = ranked[0]
    runner = ranked[1][1] if len(ranked) > 1 else 0
    if n >= min_count and n >= margin * max(runner, 1):
        return top
    return ""


def protein_regulation(gene: str, protein_alias: str = "", top_per_axis: int = 6) -> dict:
    """Everything the corpus says about what sets this protein's level.

    Use on the causal gene of a loss-of-function disease *before* proposing any
    intervention that aims to raise the product. Returns hits grouped by the
    mechanism axis that found them, so a route can be chosen from what the
    literature actually reports rather than from which search term came to mind.
    """
    g = (gene or "").strip()
    if not g:
        return {"ok": False, "error": "no gene given"}
    p = protein_alias.strip()
    auto_alias = ""
    if not p:
        auto_alias = _alias_from_corpus(g)
        p = auto_alias

    # Anchor every axis on the gene itself. Without this the broad axes drown:
    # "sirtuin SIRT deacetylase NAD" returns general NAD-metabolism reviews that
    # never mention this protein, and they outrank the one paper that does -
    # which is exactly how the SIRT2/BubR1 link was missed. We are asking what
    # regulates *this* protein, so a paper that never names it cannot answer.
    # Prefer papers that name the protein, but do not *require* it. A hard
    # filter looked right and was wrong: the gene symbol and the protein name
    # are different words, 258 abstracts here say "BubR1" and never "BUB1B",
    # and filtering on the symbol alone silently deleted the whole acetylation
    # literature - the axis came back empty and the model duly reported that
    # nothing is known. Deriving the alias automatically was worse: it returned
    # co-expressed genes and GEO accessions. So: rank, don't exclude. Papers
    # naming the protein rise to the top; a general NAD-metabolism review that
    # never mentions it gets no boost and falls on its own.
    names = {g.lower()} | ({p.lower()} if p else set())
    corpus = pd.read_parquet(INDEX)
    # Refuse rather than return noise. BM25 always returns *something*, so a
    # gene the corpus has never heard of still yields a full page of confident-
    # looking hits about other proteins - the same shape of failure as a tool
    # that answers when it should error.
    _txt = (corpus.title.fillna("") + " " + corpus.abstract.fillna("")).str.lower()
    _n = int(_txt.str.contains(g.lower(), regex=False).sum())
    if _n < 3:
        return {"ok": False, "gene": g, "n_papers_mentioning_gene": _n,
                "error": (f"the local corpus contains {_n} abstracts mentioning "
                          f"{g!r}, too few to say anything about how this protein "
                          "is regulated. This is a statement about the corpus, "
                          "which was built around MVA and the spindle assembly "
                          "checkpoint, not about the gene."),
                "source": "local Europe PMC abstract corpus"}
    text = (corpus.title.fillna("") + " " + corpus.abstract.fillna("")).str.lower()
    mentions = set(corpus.pmid[text.apply(lambda t: any(n in t for n in names))])

    axes, seen, total, dropped = {}, set(), 0, 0
    for axis, tmpl in AXES.items():
        q = tmpl.format(g=g, p=p).replace("  ", " ").strip()
        # Over-retrieve, then keep only papers that actually name the protein.
        r = literature_evidence(q, top=top_per_axis * 8)
        ranked = sorted(r.get("hits", []),
                        key=lambda h: (h["pmid"] in mentions, h["bm25_score"]),
                        reverse=True)
        hits = []
        for h in ranked[:top_per_axis]:
            total += 1
            hits.append({**h,
                         "names_the_protein": h["pmid"] in mentions,
                         "first_seen_on_axis": h["pmid"] not in seen})
            seen.add(h["pmid"])
        dropped += max(0, len(ranked) - top_per_axis)
        axes[axis] = {"query": q, "n": len(hits), "hits": hits}

    return {"ok": True, "gene": g, "protein_alias": p or None,
            "protein_alias_source": ("supplied by caller" if protein_alias.strip()
                                     else f"derived from corpus usage: {auto_alias!r}"
                                     if auto_alias else "none found"),
            "axes_searched": list(AXES),
            "n_unique_papers": len(seen), "n_hits_total": total,
            "lower_ranked_not_shown": dropped,
            "ranking": "papers naming the protein are ranked first, not filtered",
            "results_by_axis": axes,
            "interpretation": (
                "A fixed battery, run identically every time, so that whether a "
                "mechanism is found does not depend on the wording you chose. "
                "Read every axis before concluding a route is unsupported - the "
                "paper that matters is often on an axis you would not have "
                "searched, and papers reporting *opposite* effects on the same "
                "protein frequently sit on different axes (one describing a "
                "modification that stabilises, another one that promotes "
                "degradation). Where they conflict, check whether they concern "
                "different residues, compartments or timescales before deciding "
                "either is wrong. Recall is favoured over precision here, so "
                "expect irrelevant hits and discard them explicitly."),
            "source": "local Europe PMC abstract corpus via literature_evidence"}


SCHEMA = {"type": "function", "function": {
    "name": "protein_regulation",
    "description": (
        "What is known to control a protein's abundance: acetylation, sirtuins "
        "and NAD, ubiquitination and proteasomal degradation, phosphorylation, "
        "measured protein level and half-life, and transcription. Runs a fixed "
        "battery of literature searches and returns the union, so coverage does "
        "not depend on the wording you chose. Use this before proposing any "
        "intervention meant to RAISE a deficient protein, and read every axis "
        "before saying a route is unsupported."),
    "parameters": {"type": "object", "properties": {
        "gene": {"type": "string", "description": "HGNC symbol, e.g. BUB1B"},
        "protein_alias": {"type": "string",
                          "description": "Protein name if it differs, e.g. BubR1"},
        "top_per_axis": {"type": "integer", "description": "Default 6"}},
        "required": ["gene"]}}}


if __name__ == "__main__":
    r = protein_regulation("BUB1B", "BubR1")
    print(f"{r['n_unique_papers']} unique papers across {len(r['axes_searched'])} axes")
    for axis, d in r["results_by_axis"].items():
        print(f"\n--- {axis}")
        for h in d["hits"][:4]:
            print(f"   {h['pmid']} {h['year']} {h['title'][:78]}")
