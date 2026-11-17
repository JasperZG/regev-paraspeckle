"""Is the reproducible response organised along the genome?

If paraspeckle proximity is a property of where a gene sits in the nucleus rather
than of its RNA sequence, neighbouring genes should share it. Correlation of the
response between gene pairs as a function of genomic separation, against a
within-chromosome permutation baseline.
"""
import pandas as pd, numpy as np
from scipy import stats

M = pd.read_pickle("master.pkl")
res = pd.read_pickle("cqn_response.pkl")
K = M.loc[res.index].copy()
K["y1"] = res[["y_20old", "y_40odd"]].mean(axis=1)
K["y2"] = res[["y_40old", "y_40even", "y_80"]].mean(axis=1)
K["y"] = res.y_cqn
K["mid"] = (K.start + K.end) / 2
K = K.sort_values(["seq", "mid"])

EDGES = [0, 0.1, 0.25, 0.5, 1, 2, 5, 10, 25, 50, 1e4]
def profile(col, perm=False, rng=None):
    acc = {i: [[], []] for i in range(len(EDGES) - 1)}
    for ch, g in K.groupby("seq"):
        v = g[col].values.copy()
        p = g.mid.values
        if perm:
            v = rng.permutation(v)
        n = len(v)
        if n < 10: continue
        for i in range(n):
            j = i + 1
            while j < n and (p[j] - p[i]) < 50e6:
                d = (p[j] - p[i]) / 1e6
                b = np.searchsorted(EDGES, d, "right") - 1
                acc[b][0].append(v[i]); acc[b][1].append(v[j])
                j += 1
    return {b: (np.corrcoef(a, c)[0, 1] if len(a) > 500 else np.nan, len(a)) for b, (a, c) in acc.items()}

rng = np.random.default_rng(0)
real = profile("y")
perms = [profile("y", True, rng) for _ in range(5)]
print("separation (Mb)      n pairs     r(real)   r(permuted, mean+-sd)")
for b in range(len(EDGES) - 1):
    r, n = real[b]
    pv = [p[b][0] for p in perms]
    if np.isnan(r): continue
    print(f"  {EDGES[b]:6.2f} - {EDGES[b+1]:<7.2f} {n:9d}   {r:+.4f}     {np.nanmean(pv):+.4f} +- {np.nanstd(pv):.4f}")

print("\n== same, for the two batches separately (shared spatial signal?) ==")
for col in ["y1", "y2"]:
    pr = profile(col)
    print(f"  {col}: " + "  ".join(f"{EDGES[b]:.2g}-{EDGES[b+1]:.2g}Mb r={pr[b][0]:+.3f}" for b in range(4) if not np.isnan(pr[b][0])))

print("\n== cross-batch: does batch1 at gene i predict batch2 at its neighbour? ==")
a1, a2 = [], []
for ch, g in K.groupby("seq"):
    v1 = g.y1.values; v2 = g.y2.values; p = g.mid.values
    for i in range(len(v1) - 1):
        if p[i + 1] - p[i] < 1e6:
            a1.append(v1[i]); a2.append(v2[i + 1])
print(f"  neighbour pairs <1 Mb: n={len(a1)}  r(batch1[i], batch2[i+1]) = {np.corrcoef(a1, a2)[0,1]:+.4f}")

print("\n== published HyPro-seq NEAT1 targets (Yap et al. 2022) ==")
neat = M[M.gene_name == "NEAT1"].iloc[0]; NPOS = (neat.start + neat.end) / 2
for g in ["NAA40", "CCDC57", "PC", "EML3", "SF3B2"]:
    r = K[K.gene_name == g]
    if len(r):
        r = r.iloc[0]
        d = abs(r.mid - NPOS) / 1e6 if r.seq == "chr11" else np.nan
        print(f"  {g:8s} {r.seq:6s} dNEAT1={d if r.seq=='chr11' else float('nan'):6.2f} Mb  "
              f"baseMean={r.baseMean:8.1f}  y={r.y:+.3f}  y1={r.y1:+.3f}  y2={r.y2:+.3f}  "
              f"pct={stats.percentileofscore(K.y, r.y):.0f}")
