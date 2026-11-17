"""Conditional-quantile-style normalisation on GC and length.

Each library gets its own smooth log-count ~ s(GC) + s(log length) curve removed,
which is what CQN/EDASeq do. Response afterwards:
    y = mean(residual over the 10 probe libraries) - mean(residual over the 2 neg)
Any feature that still predicts y is not riding the library-level GC/length axis.
"""
import pandas as pd, numpy as np
from scipy import stats
from numpy.polynomial import polynomial as P

M = pd.read_pickle("master.pkl")
cnt = pd.read_csv("data/counts.tsv", sep="\t", index_col=0)
K = M[M.keep].copy()
c = cnt.loc[K.index]
ok = c.min(axis=1) >= 5
c = c[ok]; Ksub = K[ok]
lib = cnt.sum()
y = np.log2(c.div(lib, axis=1) * 1e6 + 1)

def bs(x, df=5):
    """natural-ish basis: quantile knots, truncated cubic + linear terms"""
    x = np.asarray(x, float)
    q = np.quantile(x, np.linspace(0, 1, df + 1)[1:-1])
    cols = [np.ones_like(x), x, x**2, x**3] + [np.clip(x - k, 0, None) ** 3 for k in q]
    return np.column_stack(cols)

X = np.column_stack([bs(Ksub.gc_gene.values), bs(Ksub.log10_len_gene.values)[:, 1:],
                     bs(Ksub.log10_len_exon.values)[:, 1:]])
R = pd.DataFrame(index=y.index, columns=y.columns, dtype=float)
for s in y.columns:
    b, *_ = np.linalg.lstsq(X, y[s].values, rcond=None)
    R[s] = y[s].values - X @ b
negs = ["NEAT1_neg_Rep1", "NEAT1_neg_Rep2"]
probe = [s for s in y.columns if s not in negs]
res = pd.DataFrame(index=R.index)
res["y_cqn"] = R[probe].mean(axis=1) - R[negs].mean(axis=1)
for lab, cols in [("20old", ["NEAT1_20old_Rep1","NEAT1_20old_Rep2"]), ("40old", ["NEAT1_40old_Rep1","NEAT1_40old_Rep2"]),
                  ("40odd", ["NEAT1_40odd_Rep1","NEAT1_40odd_Rep2"]), ("40even", ["NEAT1_40even_Rep1","NEAT1_40even_Rep2"]),
                  ("80", ["NEAT1_80_Rep1","NEAT1_80_Rep2"])]:
    res["y_" + lab] = R[cols].mean(axis=1) - R[negs].mean(axis=1)
# per-replicate pairing: each probe rep against one neg rep (independent noise)
res["y_repA"] = R[["NEAT1_20old_Rep1","NEAT1_40old_Rep1","NEAT1_40odd_Rep1","NEAT1_40even_Rep1","NEAT1_80_Rep1"]].mean(axis=1) - R["NEAT1_neg_Rep1"]
res["y_repB"] = R[["NEAT1_20old_Rep2","NEAT1_40old_Rep2","NEAT1_40odd_Rep2","NEAT1_40even_Rep2","NEAT1_80_Rep2"]].mean(axis=1) - R["NEAT1_neg_Rep2"]
res.to_pickle("cqn_response.pkl")
R.to_pickle("cqn_resid.pkl")

print("genes:", len(res))
print("residual GC/length association per library (should be ~0):")
for s in [negs[0], probe[0], probe[-1]]:
    print(f"  {s:20s} rho(GC)={stats.spearmanr(Ksub.gc_gene, R[s])[0]:+.3f} rho(len)={stats.spearmanr(Ksub.log10_len_gene, R[s])[0]:+.3f}")
print("\nreproducibility of the CQN response:")
print("  rep A vs rep B (fully independent libraries): r = %.3f  rho = %.3f"
      % (np.corrcoef(res.y_repA, res.y_repB)[0, 1], stats.spearmanr(res.y_repA, res.y_repB)[0]))
print("  40even vs 40odd: r = %.3f" % np.corrcoef(res.y_40even, res.y_40odd)[0, 1])
print("  across probe sets:")
print(res[["y_20old", "y_40old", "y_40odd", "y_40even", "y_80"]].corr().round(3).to_string())
print("\nresponse vs the old uncorrected LFC: r = %.3f" % np.corrcoef(res.y_cqn, Ksub.lfc)[0, 1])
print("residual GC/len association of the response: rho(GC)=%.3f rho(len)=%.3f"
      % (stats.spearmanr(Ksub.gc_gene, res.y_cqn)[0], stats.spearmanr(Ksub.log10_len_gene, res.y_cqn)[0]))
top = res.y_cqn.sort_values(ascending=False).head(30)
print("\ntop 30 by CQN response:")
print(pd.DataFrame({"gene": M.gene_name.reindex(top.index), "type": M.gene_type.reindex(top.index),
                    "y": top.round(2), "bm": M.baseMean.reindex(top.index).round(0),
                    "repA": res.y_repA.reindex(top.index).round(2), "repB": res.y_repB.reindex(top.index).round(2)}).to_string())
