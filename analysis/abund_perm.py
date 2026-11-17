"""Same 2-vs-10 permutation test, now for the abundance trend and for the
lncRNA / antisense enrichment that survived CQN.

Also: is the gene-level counting exon-based or gene-body-based?
"""
import pandas as pd, numpy as np, itertools
from scipy import stats

M = pd.read_pickle("master.pkl")
cnt = pd.read_csv("data/counts.tsv", sep="\t", index_col=0)
K = M[M.keep]
c = cnt.loc[K.index]
ok = c.min(axis=1) >= 5
cp = (c / c.sum() * 1e6)[ok]
Ks = K[ok]
samples = list(cnt.columns)
negs = ["NEAT1_neg_Rep1", "NEAT1_neg_Rep2"]

print("== counting mode ==")
print("  rho(log10 cpm_neg, log10 exon length) = %+.3f" % stats.spearmanr(Ks.log10_len_exon, np.log10(cp[negs].mean(axis=1) + .01))[0])
print("  rho(log10 cpm_neg, log10 intron length) = %+.3f" % stats.spearmanr(Ks.log10_len_intron, np.log10(cp[negs].mean(axis=1) + .01))[0])
print("  rho(log10 cpm_neg, log10 gene length) = %+.3f" % stats.spearmanr(Ks.log10_len_gene, np.log10(cp[negs].mean(axis=1) + .01))[0])

# CQN residual basis identical to cqn.py so the test matches the response used
def bs(x, df=5):
    x = np.asarray(x, float)
    q = np.quantile(x, np.linspace(0, 1, df + 1)[1:-1])
    return np.column_stack([np.ones_like(x), x, x**2, x**3] + [np.clip(x - k, 0, None)**3 for k in q])
X = np.column_stack([bs(Ks.gc_gene.values), bs(Ks.log10_len_gene.values)[:, 1:], bs(Ks.log10_len_exon.values)[:, 1:]])
y = np.log2(cp + 1)
R = pd.DataFrame({s: y[s].values - X @ np.linalg.lstsq(X, y[s].values, rcond=None)[0] for s in samples}, index=cp.index)

lab_bm = np.log10(Ks.baseMean.values)
is_lnc = (Ks.gene_type == "lncRNA").astype(float).values
is_as = Ks.gene_name.astype(str).str.contains("-AS|-DT|^LINC").astype(float).values

rows = []
for pair in itertools.combinations(samples, 2):
    rest = [s for s in samples if s not in pair]
    d = R[rest].mean(axis=1).values - R[list(pair)].mean(axis=1).values
    rows.append(dict(pair=" + ".join(p.replace("NEAT1_", "") for p in pair),
                     is_true=set(pair) == set(negs),
                     same_cond=pair[0].rsplit("_", 1)[0] == pair[1].rsplit("_", 1)[0],
                     rho_bm=stats.spearmanr(lab_bm, d)[0],
                     d_lnc=d[is_lnc == 1].mean() - d[is_lnc == 0].mean(),
                     d_as=d[is_as == 1].mean() - d[is_as == 0].mean()))
P = pd.DataFrame(rows)
true = P[P.is_true].iloc[0]
print("\n== permutation over 2-vs-10 sample splits, on CQN residuals ==")
for col, nicer in [("rho_bm", "rho(log baseMean, response)"), ("d_lnc", "lncRNA - other mean response"),
                   ("d_as", "antisense/divergent - other")]:
    v = P[col].values; t = true[col]
    side = (v <= t).mean() if t < 0 else (v >= t).mean()
    print(f"  {nicer:32s} true={t:+.4f}  rank {int((np.abs(v)>=abs(t)).sum())}/66 by |.|  one-sided p={side:.3f}  "
          f"null SD={v.std():.4f} range[{v.min():+.3f},{v.max():+.3f}]")
print("\n  replicate-pair splits only:")
print(P[P.same_cond][["pair", "is_true", "rho_bm", "d_lnc", "d_as"]].round(4).to_string(index=False))

print("\n== chr11<2Mb effect under the same permutation ==")
neat = M[M.gene_name == "NEAT1"].iloc[0]; NPOS = (neat.start + neat.end) / 2
near = ((Ks.seq == "chr11") & (((Ks.start + Ks.end) / 2 - NPOS).abs() < 2e6)).values
print("  n near =", near.sum())
vals = []
for pair in itertools.combinations(samples, 2):
    rest = [s for s in samples if s not in pair]
    d = R[rest].mean(axis=1).values - R[list(pair)].mean(axis=1).values
    vals.append((set(pair) == set(negs), d[near].mean() - d[~near].mean()))
V = pd.DataFrame(vals, columns=["is_true", "diff"])
t = V[V.is_true]["diff"].iloc[0]
print(f"  true diff={t:+.4f}  rank {int((V['diff'].abs()>=abs(t)).sum())}/66  one-sided p={(V['diff']>=t).mean():.3f}  null SD={V['diff'].std():.4f}")
