"""Evo 2 7B features for ClinVar-subset variants: ref/alt log-likelihood delta + embeddings.

Per variant (window of --window bp, variant start at the centre):
  delta_ll  sum log P(alt window) - sum log P(ref window)
  emb       (4*4096,) layer LAYER: [ref_site, alt_site, ref_mean, alt_mean] (fp16)
Indels: alt window = ref left context + ALT + ref right context, trimmed to --window bp.

Output (resumable; finished chunks are skipped):
  <out>/shard{i}of{n}_chunk{k}.npz   var_key, delta_ll, emb

Usage: .venv-evo2/bin/python scripts/evo2_features.py --subset <clinvar_subset.tsv.gz> --fasta <genome.fa> \
         --ckpt <evo2_7b.pt> --out <dir> [--shard 0/4] [--window 4096] [--batch 4] [--kernels] [--limit N]
"""
import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import torch

LAYER = "blocks.28.mlp.l3"
SPLIT_ORDER = {"val": 0, "test": 1, "train": 2}  # val/test first so evaluation can start early


class Evo2Scorer:
    def __init__(self, ckpt, fasta, window, kernels=False):
        from evo2 import Evo2
        self.model = Evo2("evo2_7b", local_path=ckpt, use_kernels=kernels)
        self.tok = self.model.tokenizer
        self.fa = pysam.FastaFile(fasta)
        self.window = window

    def windows(self, chrom, pos, ref, alt):
        w, c = self.window, self.window // 2
        start = pos - 1 - c
        # Extra right context so the alt window can still be filled to w bp after a deletion;
        # pad with N where the window runs past a chromosome end
        end = start + w + len(ref)
        seq = self.fa.fetch(chrom, max(start, 0), end).upper()
        seq = "N" * max(-start, 0) + seq
        seq += "N" * (end - start - len(seq))
        assert seq[c:c + len(ref)] == ref, f"REF mismatch {chrom}:{pos}"
        alt_seq = (seq[:c] + alt + seq[c + len(ref):])[:w]
        return seq[:w], alt_seq

    @torch.no_grad()
    def score(self, variants):
        """variants: list of (chrom, pos, ref, alt). Returns delta_ll (n,), emb (n, 4*4096) fp16."""
        seqs = [s for v in variants for s in self.windows(*v)]  # ref0, alt0, ref1, alt1, ...
        ids = torch.tensor([self.tok.tokenize(s) for s in seqs], dtype=torch.long, device="cuda:0")
        (logits, _), emb = self.model(ids, return_embeddings=True, layer_names=[LAYER])
        logp = logits[:, :-1].float().log_softmax(-1)
        ll = logp.gather(-1, ids[:, 1:, None]).squeeze(-1).sum(-1)
        delta = (ll[1::2] - ll[0::2]).cpu().numpy()
        e = emb[LAYER].float()
        c = self.window // 2
        feats = torch.cat([e[0::2, c], e[1::2, c], e[0::2].mean(1), e[1::2].mean(1)], 1)
        return delta, feats.cpu().numpy().astype(np.float16)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subset", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--shard", default="0/1", help="i/n: this process handles every n-th chunk starting at i")
    ap.add_argument("--window", type=int, default=4096)  # 4 kb: 2x faster than 8 kb, same val AUROC (benchmark)
    ap.add_argument("--batch", type=int, default=4, help="variants per forward pass (2x sequences)")
    ap.add_argument("--chunk", type=int, default=1000)
    ap.add_argument("--kernels", action="store_true", help="Vortex Triton kernels for Hyena convolutions")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    shard_i, shard_n = map(int, args.shard.split("/"))

    df = pd.read_csv(args.subset, sep="\t", dtype={"chrom": str}, usecols=["var_key", "chrom", "pos", "ref", "alt", "split"])
    df = df.assign(o=df.split.map(SPLIT_ORDER)).sort_values(["o", "chrom", "pos"], kind="stable").reset_index(drop=True)
    if args.limit:
        df = df.head(args.limit)
    chunks = [df.iloc[i:i + args.chunk] for i in range(0, len(df), args.chunk)]
    mine = [(k, c) for k, c in enumerate(chunks) if k % shard_n == shard_i
            and not (out / f"shard{shard_i}of{shard_n}_chunk{k:05d}.npz").exists()]
    print(f"{len(df):,} variants, {len(chunks)} chunks, {len(mine)} to do for shard {args.shard}", flush=True)
    if not mine:
        return

    scorer = Evo2Scorer(args.ckpt, args.fasta, args.window, args.kernels)
    t0, done = time.time(), 0
    for k, c in mine:
        vs = list(c[["chrom", "pos", "ref", "alt"]].itertuples(index=False, name=None))
        deltas, embs = [], []
        for i in range(0, len(vs), args.batch):
            d, e = scorer.score(vs[i:i + args.batch])
            deltas.append(d)
            embs.append(e)
        np.savez(out / f"shard{shard_i}of{shard_n}_chunk{k:05d}.npz", var_key=c.var_key.to_numpy(),
                 delta_ll=np.concatenate(deltas), emb=np.concatenate(embs))
        done += len(c)
        rate = done / (time.time() - t0)
        left = sum(len(cc) for kk, cc in mine if kk > k)
        print(f"chunk {k}: {rate:.2f} variants/s, {left:,} left (~{left / rate / 3600:.1f} h)", flush=True)


if __name__ == "__main__":
    main()
