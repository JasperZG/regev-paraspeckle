"""Is the NEAT1-window effect just residual GC over-correction?

chr11q13 has median gene GC 0.559 against 0.436 genome-wide, and the next best
windows in the genome-wide scan were on chr16, also GC-rich. So the window has to
be compared against GC-matched windows, not against the whole genome.
"""
import pandas as pd, numpy as np
from scipy import stats

M = pd.read_pickle("master.pkl")
res = pd.read_pickle("cqn_response.pkl")
K = M.loc[res.index].copy()
K["y"] = res.y_cqn
K["mid"] = (K.start + K.end) / 2
n1 = M[M.gene_name == "NEAT1"].iloc[0]; NPOS = (n1.start + n1.end) / 2
near = (K.seq == "chr11") & ((K.mid - NPOS).abs() < 2e6)
obs = K.y[near].mean() - K.y[~near].mean()
gcN = K.gc_gene[near].median()
print("NEAT1 window: n=%d  effect=%+.4f  median GC=%.3f (genome %.3f)" % (near.sum(), obs, gcN, K.gc_gene.median()))

rows = []
for ch, g in K.groupby("seq"):
    p = g.mid.values
    for cc in np.arange(p.min(), p.max(), 1e6):
        m = np.abs(p - cc) < 2e6
        if 60 <= m.sum() <= 160:
            rows.append((ch, cc, m.sum(), g.y.values[m].mean() - K.y[~((K.seq == ch) & (np.abs(K.mid - cc) < 2e6))].mean(),
                         np.median(g.gc_gene.values[m])))
S = pd.DataFrame(rows, columns=["chrom", "centre", "n", "eff", "gc"])
S["is_neat"] = (S.chrom == "chr11") & (np.abs(S.centre - NPOS) < 1.1e6)
print("\nwindows: n=%d ; correlation between window GC and window effect: rho=%+.3f" %
      (len(S), stats.spearmanr(S.gc, S.eff)[0]))
hi = S[S.gc > 0.50]
print("GC-matched comparison (window median GC > 0.50): n=%d windows" % len(hi))
print(hi.assign(Mb=(hi.centre / 1e6).round(1)).sort_values("eff", ascending=False)[["chrom", "Mb", "n", "gc", "eff", "is_neat"]].round(3).to_string(index=False))
nb = S[S.is_neat].eff.max()
print("\nNEAT1 window effect %+.4f ; rank among GC-matched windows: %d/%d" %
      (nb, int((hi.eff >= nb).sum()), len(hi)))

print("\n== within-window test: does GC still predict y inside the NEAT1 neighbourhood? ==")
nw = K[near]
print("  rho(GC, y) inside the window = %+.3f (n=%d)" % (stats.spearmanr(nw.gc_gene, nw.y)[0], len(nw)))
print("  rho(GC, y) genome-wide after CQN = %+.3f" % stats.spearmanr(K.gc_gene, K.y)[0])

print("\n== gene-by-gene GC-matched control for the 93 window genes ==")
rng = np.random.default_rng(0)
pool = K[~near]
diffs = []
for _ in range(2000):
    picks = []
    for _, r in nw.iterrows():
        cand = pool[(pool.gc_gene - r.gc_gene).abs() < 0.015]
        if len(cand):
            picks.append(cand.y.values[rng.integers(len(cand))])
    diffs.append(nw.y.mean() - np.mean(picks))
diffs = np.array(diffs)
print("  mean difference vs GC-matched genes: %+.4f  95%% CI [%+.4f, %+.4f]  p=%.3f"
      % (diffs.mean(), np.percentile(diffs, 2.5), np.percentile(diffs, 97.5), (diffs <= 0).mean()))

print("\n== and matched on GC + length + expression ==")
pool2 = pool.copy()
diffs2 = []
for _ in range(1000):
    picks = []
    for _, r in nw.iterrows():
        cand = pool2[((pool2.gc_gene - r.gc_gene).abs() < 0.02) &
                     ((pool2.log10_len_gene - r.log10_len_gene).abs() < 0.3) &
                     ((pool2.log10_baseMean - r.log10_baseMean).abs() < 0.4)]
        if len(cand):
            picks.append(cand.y.values[rng.integers(len(cand))])
    diffs2.append(nw.y.mean() - np.mean(picks))
diffs2 = np.array(diffs2)
print("  mean difference: %+.4f  95%% CI [%+.4f, %+.4f]  p=%.3f  (matched n per draw ~%d)"
      % (diffs2.mean(), np.percentile(diffs2, 2.5), np.percentile(diffs2, 97.5), (diffs2 <= 0).mean(), len(nw)))
