"""Does RBP binding predict the paraspeckle response?

Peak counts scale with transcript length and expression, so everything is tested
as a density (peaks per kb of the region) and partialled on length + abundance,
and every candidate must replicate in both batches. Non-paraspeckle RBPs act as
the specificity control.
"""
import pandas as pd, numpy as np
from scipy import stats

M = pd.read_pickle("master.pkl")
res = pd.read_pickle("cqn_response.pkl")
E = pd.read_pickle("eclip_features.pkl")
K = M.loc[res.index].copy()
K["y"] = res.y_cqn
K["y1"] = res[["y_20old", "y_40odd"]].mean(axis=1)
K["y2"] = res[["y_40old", "y_40even", "y_80"]].mean(axis=1)

def rr(x, covs):
    rx = stats.rankdata(x)
    X = np.column_stack([np.ones(len(x))] + [stats.rankdata(c) for c in covs])
    b, *_ = np.linalg.lstsq(X, rx, rcond=None)
    return rx - X @ b

PARA = {"NONO", "SFPQ", "FUS", "TARDBP", "HNRNPK", "MATR3", "HNRNPU", "EWSR1", "TAF15", "AKAP8L", "ZC3H11A"}
rows = []
for col in E.columns:
    lab = col.rsplit("_", 2)[0]
    rbp = lab.split("_")[0]
    reg = col.rsplit("_", 2)[1]
    L = (K.len_gene if reg == "gene" else K.len_exon) / 1000
    dens = E[col].reindex(K.index) / L.clip(lower=0.1)
    if (dens > 0).sum() < 300:
        continue
    cov = [K.log10_len_gene.values, K.log10_len_exon.values, K.log10_baseMean.values]
    rows.append(dict(rbp=rbp, panel=lab, region=reg, ngene=int((E[col] > 0).sum()),
                     rho=stats.spearmanr(dens, K.y)[0],
                     rho1=stats.spearmanr(dens, K.y1)[0], rho2=stats.spearmanr(dens, K.y2)[0],
                     partial=np.corrcoef(rr(dens.values, cov), rr(K.y.values, cov))[0, 1],
                     paraspeckle=rbp in PARA))
T = pd.DataFrame(rows)
T["replicates"] = np.where(np.sign(T.rho1) == np.sign(T.rho2), np.minimum(T.rho1.abs(), T.rho2.abs()), 0) * np.sign(T.rho1)
T = T.sort_values("partial", ascending=False)
pd.set_option("display.width", 220)
print(T.round(3).to_string(index=False))
T.to_csv("res_eclip.tsv", sep="\t", index=False)

print("\n== paraspeckle proteins vs control RBPs (exon density, partial rho) ==")
for grp, s in T[T.region == "exon"].groupby("paraspeckle"):
    print(f"  {'paraspeckle' if grp else 'control   '}: n={len(s)} median partial={s.partial.median():+.3f} "
          f"max={s.partial.max():+.3f} ({s.loc[s.partial.idxmax(), 'panel'] if len(s) else '-'})")

best = T[(T.region == "exon") & T.paraspeckle].head(4).panel.tolist()
print("\n== binary: bound vs unbound, for the strongest panels ==")
for lab in best:
    col = f"{lab}_exon_npk"
    b = E[col].reindex(K.index) > 0
    print(f"  {lab:16s} bound n={b.sum():5d} mean y={K.y[b].mean():+.3f} vs unbound {K.y[~b].mean():+.3f}  "
          f"p={stats.mannwhitneyu(K.y[b], K.y[~b]).pvalue:.1e}")
