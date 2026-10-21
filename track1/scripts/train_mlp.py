"""MLP on scalar features + full embedding blocks (GPU). Early stopping on val AUROC; 3-seed ensemble.

Selection on val only. Test predictions are saved but not printed (see evaluate_test.py).

Usage: .venv-evo2/bin/python scripts/train_mlp.py <cache_dir> <out_dir> <feature_set>
  feature_set: all | evo | esm   (scalars are always included)
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn

sys.path.insert(0, str(Path(__file__).parent))
from model_common import consequence_vocab, emb_blocks, evaluate, load_cache, scalar_features  # noqa: E402

BLOCKS = {
    "all": ["evo_site_diff", "evo_mean_diff", "evo_ref_site", "esm_site_diff", "esm_mean_diff", "esm_wt_site"],
    "evo": ["evo_site_diff", "evo_mean_diff", "evo_ref_site"],
    "esm": ["esm_site_diff", "esm_mean_diff", "esm_wt_site"],
}
SEEDS = [0, 1, 2]
EPOCHS, PATIENCE, BATCH = 40, 5, 1024
LR, WD, DROPOUT, HIDDEN = 3e-4, 0.05, 0.3, (512, 128)

cache, out, fs = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
out.mkdir(parents=True, exist_ok=True)
dev = "cuda"
t0 = time.time()
meta, evo, esm = load_cache(cache)
vocab = consequence_vocab(meta)
split, y = meta.split.to_numpy(), meta.label.to_numpy().astype(np.float32)
tr, va = split == "train", split == "val"

blocks = emb_blocks(evo, esm, np.arange(len(meta)))
X = np.hstack([scalar_features(meta, vocab).to_numpy()] + [blocks[b] for b in BLOCKS[fs]]).astype(np.float32)
del blocks
mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-6
X = torch.from_numpy((X - mu) / sd).to(dev, torch.float16)  # ~9 GB for "all"
Y = torch.from_numpy(y).to(dev)
tr_idx = torch.from_numpy(np.where(tr)[0]).to(dev)
print(f"[{fs}] X {tuple(X.shape)} on GPU, loaded in {time.time() - t0:.0f}s", flush=True)


def make_model(d_in):
    layers, d = [], d_in
    for h in HIDDEN:
        layers += [nn.Linear(d, h), nn.GELU(), nn.Dropout(DROPOUT)]
        d = h
    return nn.Sequential(*layers, nn.Linear(d, 1)).to(dev)


@torch.no_grad()
def predict(model):
    model.eval()
    out_ = [torch.sigmoid(model(X[i:i + 8192].float()).squeeze(1)) for i in range(0, len(X), 8192)]
    return torch.cat(out_).cpu().numpy()


val_meta = meta[va].reset_index(drop=True)
ens = np.zeros(len(meta))
for seed in SEEDS:
    torch.manual_seed(seed)
    model = make_model(X.shape[1])
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
    best, best_state, bad = -1, None, 0
    for ep in range(EPOCHS):
        model.train()
        perm = tr_idx[torch.randperm(len(tr_idx), device=dev)]
        for i in range(0, len(perm), BATCH):
            b = perm[i:i + BATCH]
            loss = nn.functional.binary_cross_entropy_with_logits(model(X[b].float()).squeeze(1), Y[b])
            opt.zero_grad()
            loss.backward()
            opt.step()
        auc = evaluate(val_meta, predict(model)[va], "")[0]["auroc"]
        if auc > best:
            best, bad = auc, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
        print(f"[{fs}] seed {seed} epoch {ep}: val AUROC {auc:.4f}", flush=True)
        if bad >= PATIENCE:
            break
    model.load_state_dict(best_state)
    ens += predict(model) / len(SEEDS)
    torch.save({"state": best_state, "mu": mu, "sd": sd, "blocks": BLOCKS[fs], "vocab": vocab,
                "hidden": HIDDEN, "dropout": DROPOUT}, out / f"mlp_{fs}_seed{seed}.pt")
    print(f"[{fs}] seed {seed}: best val AUROC {best:.4f}  [{time.time() - t0:.0f}s]", flush=True)

name = f"MLP[scalars+{fs}]"
pd.DataFrame({"var_key": meta.var_key, "split": split, "label": y, name: ens}).to_parquet(
    out / f"preds_mlp_{fs}.parquet", index=False)
res = pd.DataFrame(evaluate(val_meta, ens[va], name))
res.to_csv(out / f"val_results_mlp_{fs}.tsv", sep="\t", index=False)
print(res.round(4).to_string(index=False))
