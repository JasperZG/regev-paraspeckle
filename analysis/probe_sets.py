"""Structure across probe sets.

80 = all oligos; 40even/40odd = disjoint halves of it; 20old/40old = an earlier design.
A gene enriched because one particular oligo cross-hybridises with it should appear
in the sets containing that oligo and not the others. A gene enriched because it is
near paraspeckles should appear in all of them.
"""
import pandas as pd, numpy as np
from scipy import stats

M = pd.read_pickle("master.pkl")
res = pd.read_pickle("cqn_response.pkl")
sets = ["y_20old", "y_40old", "y_40odd", "y_40even", "y_80"]
Y = res[sets]
print("correlation between probe sets (CQN response):")
print(Y.corr().round(3).to_string())
print("\nhierarchical structure: 40old/40even/80 vs 20old/40odd")
print("  mean r within {40old,40even,80}: %.3f" % np.mean([Y.y_40old.corr(Y.y_40even), Y.y_40old.corr(Y.y_80), Y.y_40even.corr(Y.y_80)]))
print("  r(20old, 40odd): %.3f" % Y.y_20old.corr(Y.y_40odd))
print("  mean r across the two groups: %.3f" %
      np.mean([Y[a].corr(Y[b]) for a in ["y_40old", "y_40even", "y_80"] for b in ["y_20old", "y_40odd"]]))

Z = (Y - Y.mean()) / Y.std()
res2 = pd.DataFrame(index=Y.index)
res2["consistent"] = Z.min(axis=1)          # enriched in the weakest set too
res2["spread"] = Z.max(axis=1) - Z.min(axis=1)
res2["mean"] = Z.mean(axis=1)
nm = M.gene_name.reindex(Y.index)
print("\n== genes enriched in EVERY probe set (top 25 by the weakest set) ==")
t = res2.sort_values("consistent", ascending=False).head(25)
print(pd.DataFrame({"gene": nm[t.index], "type": M.gene_type.reindex(t.index)}).join(Z.round(2)).join(t[["consistent"]].round(2)).to_string())
print("\n== most set-specific genes (largest spread, top 15) ==")
t = res2[res2["mean"] > 0].sort_values("spread", ascending=False).head(15)
print(pd.DataFrame({"gene": nm[t.index]}).join(Z.round(2)).join(t[["spread"]].round(2)).to_string())
res2.to_pickle("probeset_structure.pkl")

print("\n== how much of the response is shared across all five sets? ==")
u, s, vt = np.linalg.svd(Z.fillna(0).values - Z.fillna(0).values.mean(0), full_matrices=False)
print("  variance explained by PC1..PC3:", np.round(s**2 / (s**2).sum(), 3)[:3])
print("  PC1 loadings:", dict(zip(sets, np.round(vt[0], 2))))
print("  PC2 loadings:", dict(zip(sets, np.round(vt[1], 2))))
