"""What separates {40old,40even,80} from {20old,40odd}?

If it is a capture-oligo effect it should be gene-intrinsic and scattered.
If it is a batch / handling effect it should look like a global sample axis.
"""
import pandas as pd, numpy as np
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram

M = pd.read_pickle("master.pkl")
R = pd.read_pickle("cqn_resid.pkl")
groups = pd.read_pickle("probe_groups.pkl")
short = [c.replace("NEAT1_", "") for c in R.columns]
R.columns = short

print("== sample-sample correlation of CQN residuals ==")
C = R.corr(method="spearman")
print(C.round(2).to_string())
Zc = linkage(1 - C.values, "average")
d = dendrogram(Zc, labels=list(C.columns), no_plot=True)
print("\nclustering order:", d["ivl"])

X = (R - R.mean(axis=1).values[:, None]).values
u, s, vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
print("\nPCA of libraries (var explained):", np.round(s**2 / (s**2).sum(), 3)[:4])
for i in range(3):
    print(f"  PC{i+1}:", dict(zip(short, np.round(vt[i], 2))))

print("\n== are the even-only genes clustered in the genome? ==")
for lab in ["even_only", "odd_only", "both_halves", "flat_control"]:
    idx = groups[lab]
    sub = M.loc[idx]
    per_chr = sub.seq.value_counts(normalize=True)
    # nearest-neighbour distance within the group
    dmin = []
    for ch, g in sub.groupby("seq"):
        p = np.sort(((g.start + g.end) / 2).values)
        if len(p) > 1:
            dd = np.diff(p); dmin += list(np.minimum(np.r_[dd, np.inf], np.r_[np.inf, dd]))
    print(f"  {lab:14s} n={len(sub):4d} top chrom={per_chr.index[0]} ({per_chr.iloc[0]:.2f})  "
          f"median nn-dist={np.median(dmin)/1e6:.1f} Mb  pairs <1Mb: {np.mean(np.array(dmin) < 1e6):.3f}")

print("\n== expression level and biotype of each group ==")
for lab in ["even_only", "odd_only", "both_halves", "flat_control"]:
    sub = M.loc[groups[lab]]
    print(f"  {lab:14s} median baseMean={sub.baseMean.median():8.1f}  lncRNA frac={np.mean(sub.gene_type=='lncRNA'):.2f} "
          f" median GC={sub.gc_gene.median():.3f}  median log10 len={sub.log10_len_gene.median():.2f}")

print("\n== do even-only genes behave like this in the ORIGINAL (uncorrected) LFCs too? ==")
res = pd.read_pickle("cqn_response.pkl")
for lab in ["even_only", "odd_only", "both_halves"]:
    idx = groups[lab]
    print(f"  {lab:14s} median pooled apeglm LFC = {M.lfc.reindex(idx).median():+.3f}")
