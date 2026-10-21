"""Compact storage and fast scoring for fitted LaplacianNB models.

A fitted laplaciannb.LaplacianNB keeps `feature_names_` as a Python dict with
one entry per distinct fingerprint bit. Pickling ~700 of those is both slow and
large, and its `_joint_log_likelihood` densifies to (n_samples, n_features),
which does not fit for wide feature spaces.

`save` strips a fitted model down to two arrays; `LaplacianModel` reproduces the
library's scoring arithmetic exactly (see `check_matches_library`) but with
sparse lookups, so a molecule can be scored against hundreds of targets.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


BG_PCTS = (50.0, 90.0, 95.0, 99.0, 99.9)


def save(path: str, model, meta: dict | None = None,
         background: tuple[float, float] = (0.0, 1.0),
         background_scores: np.ndarray | None = None) -> None:
    """Serialise a fitted laplaciannb.LaplacianNB to an .npz.

    `background` is (mean, std) of this model's decision function over a random
    molecule sample; it puts otherwise incomparable per-target scores on a
    common scale, which is what makes cross-target ranking meaningful.

    `background_scores` is the raw sample those statistics came from. Its
    percentiles are stored so a verdict can be stated as "above this target's
    own 99th background percentile" rather than assuming the scores are normal.
    That matters once low-support targets are included: with 30 actives the
    score distribution is heavier-tailed, so a flat z cut means something
    different for every target.
    """
    ids = np.fromiter(model.feature_names_.keys(), dtype=np.int64,
                      count=len(model.feature_names_))
    cols = np.fromiter(model.feature_names_.values(), dtype=np.int64,
                       count=len(model.feature_names_))
    order = np.argsort(ids, kind="stable")       # numeric order, for searchsorted
    np.savez_compressed(
        path,
        feature_ids=ids[order].astype(np.uint32),
        log_prob=np.ascontiguousarray(model.feature_log_prob_[:, cols[order]],
                                      dtype=np.float32),
        classes=np.asarray(model.classes_),
        class_log_prior=np.asarray(model.class_log_prior_, dtype=np.float64),
        background=np.asarray(background, dtype=np.float64),
        background_pcts=np.asarray(BG_PCTS, dtype=np.float64),
        background_quantiles=np.asarray(
            np.percentile(background_scores, BG_PCTS)
            if background_scores is not None and len(background_scores)
            else [np.nan] * len(BG_PCTS), dtype=np.float64),
        meta=np.array(repr(meta or {})),
    )


class LaplacianModel:
    """Scorer for a saved LaplacianNB."""

    def __init__(self, feature_ids, log_prob, classes, class_log_prior, meta=None,
                 background=(0.0, 1.0), background_pcts=None,
                 background_quantiles=None):
        self.feature_ids = feature_ids
        self.log_prob = log_prob
        self.classes = classes
        self.class_log_prior = class_log_prior
        self.meta = meta or {}
        self.background = background
        self.background_pcts = background_pcts
        self.background_quantiles = background_quantiles

    def background_percentile(self, scores):
        """Where a score sits in this target's own background distribution.

        Returns the highest stored percentile the score exceeds, or 0.0. This is
        the honest way to express confidence for a low-support target: it makes
        no normality assumption and is specific to that model's own scale.
        """
        import numpy as _np
        s = _np.atleast_1d(scores)
        if self.background_quantiles is None or _np.isnan(
                self.background_quantiles).all():
            return _np.zeros(len(s))
        out = _np.zeros(len(s))
        for pct, q in zip(self.background_pcts, self.background_quantiles):
            out = _np.where(s >= q, pct, out)
        return out

    @classmethod
    def load(cls, path: str) -> "LaplacianModel":
        z = np.load(path, allow_pickle=False)
        meta = eval(str(z["meta"]), {"__builtins__": {}}, {})  # noqa: S307 - our own repr
        bg = tuple(z["background"]) if "background" in z.files else (0.0, 1.0)
        pcts = z["background_pcts"] if "background_pcts" in z.files else None
        quants = (z["background_quantiles"]
                  if "background_quantiles" in z.files else None)
        return cls(z["feature_ids"], z["log_prob"], z["classes"],
                   z["class_log_prior"], meta, bg, pcts, quants)

    def joint_log_likelihood(self, X: sp.csr_matrix) -> np.ndarray:
        """Sum of per-bit log weights, matching LaplacianNB._joint_log_likelihood."""
        X = X.tocsr()
        n = X.shape[0]
        idx = X.indices.astype(np.int64, copy=False)
        pos = np.searchsorted(self.feature_ids, idx)
        np.clip(pos, 0, len(self.feature_ids) - 1, out=pos)
        known = self.feature_ids[pos] == idx          # bits unseen in training
        rows = np.repeat(np.arange(n), np.diff(X.indptr))
        r, c = rows[known], pos[known]
        jll = np.empty((n, self.log_prob.shape[0]))
        for k in range(self.log_prob.shape[0]):
            jll[:, k] = np.bincount(r, weights=self.log_prob[k][c], minlength=n)
        return jll

    def decision_function(self, X: sp.csr_matrix) -> np.ndarray:
        """Score for the positive class; higher means more likely active."""
        jll = self.joint_log_likelihood(X)
        pos = int(np.argmax(self.classes))
        return jll[:, pos] - jll[:, 1 - pos] if jll.shape[1] == 2 else jll[:, pos]

    def zscore(self, X: sp.csr_matrix) -> np.ndarray:
        """Decision function standardised against the random-molecule background,
        so scores from different targets can be ranked against each other."""
        mean, std = self.background
        return (self.decision_function(X) - mean) / (std if std > 0 else 1.0)

    def predict_proba(self, X: sp.csr_matrix) -> np.ndarray:
        jll = self.joint_log_likelihood(X)
        jll = jll - jll.max(axis=1, keepdims=True)
        p = np.exp(jll)
        return p / p.sum(axis=1, keepdims=True)


def check_matches_library(model, saved: "LaplacianModel", X: sp.csr_matrix,
                          rtol: float = 1e-4) -> None:
    """Assert the compact scorer agrees with the fitted library model."""
    ref = model.predict_proba(X)
    got = saved.predict_proba(X)
    if not np.allclose(ref, got, rtol=rtol, atol=1e-6):
        raise AssertionError(
            f"scorer disagrees with laplaciannb: max |diff| = "
            f"{np.abs(ref - got).max():.3g}")
