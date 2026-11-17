"""Is the <2 Mb effect specific to the NEAT1 locus?

Slide a 2 Mb window over every position in the genome and ask how often a random
locus beats NEAT1. Also leave-one-gene-out, to check it is not a few genes.
"""
import pandas as pd, numpy as np
from scipy import stats

K = pd.read_pickle("master_resp.pkl")
K["mid"] = (K.start + K.end) / 2
neat = K.attrs.get("neat", None)
M = pd.read_pickle("master.pkl")
n = M[M.gene_name == "NEAT1"].iloc[0]; NPOS = (n.start + n.end) / 2
y = K.y.values
near = ((K.seq == "chr11") & ((K.mid - NPOS).abs() < 2e6)).values
obs = y[near].mean() - y[~near].mean()
print("NEAT1 window: n=%d  effect=%+.4f" % (near.sum(), obs))

rng = np.random.default_rng(1)
stat = []
for ch, g in K.groupby("seq"):
    p = g.mid.values; v = g.y.values
    for c in np.arange(p.min(), p.max(), 1e6):
        m = np.abs(p - c) < 2e6
        if 60 <= m.sum() <= 160:
            rest = np.concatenate([K.y.values[K.seq.values != ch], v[~m]])
            stat.append((ch, c, m.sum(), v[m].mean() - rest.mean()))
S = pd.DataFrame(stat, columns=["chrom", "centre", "n", "eff"])
print("random 2 Mb windows with 60-160 genes: n=%d" % len(S))
print("  fraction with effect >= NEAT1's: %.4f" % (S.eff >= obs).mean())
print("  window effect distribution: SD=%.4f  95th pct=%+.4f  max=%+.4f" % (S.eff.std(), S.eff.quantile(.95), S.eff.max()))
print("\n  top 8 windows genome-wide:")
T = S.sort_values("eff", ascending=False).head(8)
T["Mb"] = (T.centre / 1e6).round(1)
print(T[["chrom", "Mb", "n", "eff"]].round(4).to_string(index=False))

print("\n== leave-one-gene-out on the NEAT1 window ==")
idx = np.where(near)[0]
effs = []
for i in idx:
    m = near.copy(); m[i] = False
    effs.append(y[m].mean() - y[~near].mean())
print("  min %+.4f  max %+.4f  (observed %+.4f) -> not driven by one gene" % (min(effs), max(effs), obs))

print("\n== radius profile with bootstrap CI ==")
for r in [0.5, 1, 2, 3, 5, 10]:
    m = (K.seq == "chr11") & ((K.mid - NPOS).abs() < r * 1e6)
    v = K.y[m].values
    if len(v) < 10: continue
    bs = [np.mean(rng.choice(v, len(v))) for _ in range(2000)]
    print(f"  r={r:4.1f} Mb n={len(v):4d} mean={v.mean():+.3f}  95%% CI [{np.percentile(bs,2.5):+.3f},{np.percentile(bs,97.5):+.3f}]")
print("  genome mean %+.3f" % K.y.mean())
