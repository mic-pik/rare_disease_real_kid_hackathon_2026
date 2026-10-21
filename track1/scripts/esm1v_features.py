"""ESM-1v features for missense variants: masked-marginal scores (5 models) + embeddings (model 1).

Per variant:
  scores   (5,)      log P(mt) - log P(wt) at the masked site, one per ESM-1v model (ensemble = mean)
  emb      (4*1280,) model 1 last layer: [wt_site, mt_site, wt_mean, mt_mean] (fp16)
Proteins longer than 1022 aa are cut to a 1022 aa window around the site.

QC runs over all missense variants first (no model): ENSP in peptide FASTA, position in range,
WT residue matches VEP. Only QC-passing variants are embedded.

Output (resumable; finished chunks are skipped):
  <out>/qc.tsv.gz                    one row per VEP missense variant with qc status
  <out>/shard{i}of{n}_chunk{k}.npz   var_key, scores, emb

Usage: .venv/bin/python scripts/esm1v_features.py --vep <vep.tsv> --subset <clinvar_subset.tsv.gz> \
         --pep <pep.all.fa.gz> --models <models_dir> --out <dir> [--device cuda] [--shard 0/4] [--limit N]
"""
import argparse
import gzip
import io
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmForMaskedLM

MAX_AA = 1022  # 1024 positions incl. <cls>/<eos>
AA20 = set("ACDEFGHIKLMNPQRSTVWY")
SPLIT_ORDER = {"val": 0, "test": 1, "train": 2}  # val/test first so evaluation can start early

ap = argparse.ArgumentParser()
ap.add_argument("--vep", required=True)
ap.add_argument("--subset", required=True)
ap.add_argument("--pep", required=True)
ap.add_argument("--models", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--device", default="cuda")
ap.add_argument("--shard", default="0/1", help="i/n: this process handles every n-th chunk starting at i")
ap.add_argument("--chunk", type=int, default=2000)
ap.add_argument("--batch", type=int, default=16)
ap.add_argument("--limit", type=int, default=None, help="only the first N QC-passing variants (testing)")
args = ap.parse_args()
out = Path(args.out)
out.mkdir(parents=True, exist_ok=True)
shard_i, shard_n = map(int, args.shard.split("/"))


def read_vep(path):
    txt = "".join(line for line in Path(path).read_text().splitlines(keepends=True) if not line.startswith("##")).lstrip("#")
    return pd.read_csv(io.StringIO(txt), sep="\t", dtype=str, na_values=["-"],
                       usecols=["Uploaded_variation", "Consequence", "ENSP", "Protein_position", "Amino_acids"])


def read_pep(path):
    seqs, name, buf = {}, None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0].split(".")[0], []
            else:
                buf.append(line.strip())
    seqs[name] = "".join(buf)
    # Terminal '*' is a stop marker; internal ones are not in the ESM vocab
    return {k: v.rstrip("*").replace("*", "X") for k, v in seqs.items()}


# ---- QC (all missense variants, no model) ----
vep = read_vep(args.vep)
vep = vep[vep.Consequence.str.split(",").str[0] == "missense_variant"].rename(columns={"Uploaded_variation": "var_key"})
aa = vep.Amino_acids.str.split("/", expand=True)
vep["wt"], vep["mt"] = aa[0], aa[1]
vep["pos"] = pd.to_numeric(vep.Protein_position, errors="coerce")
vep["ensp"] = vep.ENSP
pep = read_pep(args.pep)


def qc(r):
    if r.wt not in AA20 or r.mt not in AA20:
        return "multi_or_nonstandard_aa"  # e.g. ST/SA: adjacent-residue MNVs VEP calls missense
    if pd.isna(r.ensp) or r.ensp not in pep:
        return "ensp_missing"
    if pd.isna(r.pos) or not 1 <= r.pos <= len(pep[r.ensp]):
        return "pos_out_of_range"
    return "ok" if pep[r.ensp][int(r.pos) - 1] == r.wt else "wt_mismatch"


vep["qc"] = [qc(r) for r in vep.itertuples()]
vep["prot_len"] = vep.ensp.map(lambda e: len(pep.get(e, "")) if isinstance(e, str) else 0)
subset = pd.read_csv(args.subset, sep="\t", dtype={"chrom": str}, usecols=["var_key", "split", "label"])
vep = vep.merge(subset, on="var_key", how="left")
vep[["var_key", "ensp", "pos", "wt", "mt", "prot_len", "qc", "split"]].to_csv(out / "qc.tsv.gz", sep="\t", index=False)
print(f"Missense variants: {len(vep):,}")
print(vep.qc.value_counts().to_string())
ok = vep[vep.qc == "ok"]
print(f"Proteins: {ok.ensp.nunique():,}; variants in proteins > {MAX_AA} aa (windowed): {(ok.prot_len > MAX_AA).sum():,}")

# ---- Feature extraction ----
todo = ok.assign(o=ok.split.map(SPLIT_ORDER).fillna(3)).sort_values(["o", "ensp", "pos"]).reset_index(drop=True)
if args.limit:
    todo = todo.head(args.limit)
todo["pos"] = todo.pos.astype(int)
chunks = [todo.iloc[i:i + args.chunk] for i in range(0, len(todo), args.chunk)]
mine = [(k, c) for k, c in enumerate(chunks) if k % shard_n == shard_i]
mine = [(k, c) for k, c in mine if not (out / f"shard{shard_i}of{shard_n}_chunk{k:05d}.npz").exists()]
print(f"Chunks: {len(chunks)} total, {len(mine)} to do for shard {args.shard}")
if not mine:
    raise SystemExit

model_dirs = sorted(Path(args.models).glob("esm1v_t33_650M_UR90S_*"))
tok = AutoTokenizer.from_pretrained(model_dirs[0])
dtype = torch.float16 if args.device.startswith("cuda") else torch.float32
models = [EsmForMaskedLM.from_pretrained(d, dtype=dtype).to(args.device).eval() for d in model_dirs]
wt_ids = {a: tok.convert_tokens_to_ids(a) for a in AA20}


def window(ensp, pos):
    seq = pep[ensp]
    start = 0 if len(seq) <= MAX_AA else min(max(pos - 1 - MAX_AA // 2, 0), len(seq) - MAX_AA)
    return seq[start:start + MAX_AA], pos - 1 - start


def batches(items):
    for i in range(0, len(items), args.batch):
        yield items[i:i + args.batch]


@torch.no_grad()
def run(model, seqs, sites, mask=False, hidden=False):
    enc = tok(seqs, return_tensors="pt", padding=True).to(args.device)
    ids = enc["input_ids"].clone()
    rows = torch.arange(len(seqs), device=args.device)
    tok_sites = torch.tensor(sites, device=args.device) + 1  # +1 for <cls>
    if mask:
        ids[rows, tok_sites] = tok.mask_token_id
    o = model(input_ids=ids, attention_mask=enc["attention_mask"], output_hidden_states=hidden)
    res = {"logp": o.logits[rows, tok_sites].float().log_softmax(-1)}
    if hidden:
        h = o.hidden_states[-1].float()
        m = enc["attention_mask"].clone()
        m[:, 0] = 0
        m[rows, enc["attention_mask"].sum(1) - 1] = 0  # drop <cls>/<eos> from the mean
        res["site"] = h[rows, tok_sites]
        res["mean"] = (h * m[..., None]).sum(1) / m.sum(1, keepdim=True)
    return res


@torch.no_grad()
def wt_states(seq):
    """Model-1 last-layer states for an unmodified window: (per-residue [L, 1280], mean [1280])."""
    enc = tok([seq], return_tensors="pt").to(args.device)
    h = models[0](**enc, output_hidden_states=True).hidden_states[-1][0, 1:-1].float()
    return h, h.mean(0)


t0 = time.time()
for k, c in mine:
    wins = [window(r.ensp, r.pos) for r in c.itertuples()]
    seqs = [w[0] for w in wins]
    sites = [w[1] for w in wins]
    mts = [s[:i] + m + s[i + 1:] for (s, i), m in zip(wins, c.mt)]
    wt_tok = torch.tensor([wt_ids[a] for a in c.wt])
    mt_tok = torch.tensor([wt_ids[a] for a in c.mt])
    idx = list(range(len(c)))

    scores = np.zeros((len(c), len(models)), dtype=np.float32)
    for mi, model in enumerate(models):
        for b in batches(idx):
            lp = run(model, [seqs[i] for i in b], [sites[i] for i in b], mask=True)["logp"].cpu()
            scores[b, mi] = (lp[range(len(b)), mt_tok[b]] - lp[range(len(b)), wt_tok[b]]).numpy()

    emb = np.zeros((len(c), 4 * 1280), dtype=np.float16)
    # WT once per protein window (chunk is sorted by protein), mutant once per variant
    wt_cache = {}
    for i in idx:
        if seqs[i] not in wt_cache:
            wt_cache[seqs[i]] = wt_states(seqs[i])
    for b in batches(idx):
        m = run(models[0], [mts[i] for i in b], [sites[i] for i in b], hidden=True)
        w_site = torch.stack([wt_cache[seqs[i]][0][sites[i]] for i in b])
        w_mean = torch.stack([wt_cache[seqs[i]][1] for i in b])
        emb[b] = torch.cat([w_site, m["site"], w_mean, m["mean"]], 1).cpu().numpy().astype(np.float16)
    del wt_cache

    np.savez(out / f"shard{shard_i}of{shard_n}_chunk{k:05d}.npz", var_key=c.var_key.to_numpy(), scores=scores, emb=emb)
    done = sum(len(cc) for kk, cc in mine if kk <= k)
    rate = done / (time.time() - t0)
    print(f"chunk {k}: {len(c)} variants, {rate:.2f} variants/s, ensemble mean score {scores.mean():.3f}", flush=True)
