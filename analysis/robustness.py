"""Independent re-derivation of the three headline results with different choices,
to check they are not artifacts of the normalisation I happened to pick."""
import pandas as pd, numpy as np, itertools
from scipy import stats

M = pd.read_pickle("master.pkl")
cnt = pd.read_csv("data/counts.tsv", sep="\t", index_col=0)
K = M[M.keep]
c = cnt.loc[K.index]
samples = list(cnt.columns); negs = ["NEAT1_neg_Rep1", "NEAT1_neg_Rep2"]
n1 = M[M.gene_name == "NEAT1"].iloc[0]; NPOS = (n1.start + n1.end) / 2

def deciles_norm(cp, Ks, nq):
    """alternative to splines: median-centre within GC x length bins (no smoothing)"""
    g = pd.qcut(Ks.gc_gene, nq, labels=False, duplicates="drop")
    l = pd.qcut(Ks.log10_len_gene, nq, labels=False, duplicates="drop")
    key = g.astype(str) + "_" + l.astype(str)
    y = np.log2(cp + 1)
    return y.groupby(key).transform(lambda v: v - v.median())

def run_variant(minc, nq, label):
    ok = c.min(axis=1) >= minc
    cp = (c / c.sum() * 1e6)[ok]; Ks = K[ok]
    R = deciles_norm(cp, Ks, nq)
    near = ((Ks.seq == "chr11") & (((Ks.start + Ks.end) / 2 - NPOS).abs() < 2e6)).values
    lnc = (Ks.gene_type == "lncRNA").values
    out = []
    for pair in itertools.combinations(samples, 2):
        rest = [s for s in samples if s not in pair]
        d = R[rest].mean(axis=1).values - R[list(pair)].mean(axis=1).values
        out.append((set(pair) == set(negs), d[near].mean() - d[~near].mean(),
                    d[lnc].mean() - d[~lnc].mean(),
                    stats.spearmanr(Ks.gc_gene.values, d)[0]))
    P = pd.DataFrame(out, columns=["t", "cis", "lnc", "gc"])
    tr = P[P.t].iloc[0]
    print(f"  {label:34s} n={ok.sum():5d}  cis {tr.cis:+.3f} (rank {int((P.cis>=tr.cis).sum())}/66)   "
          f"lnc {tr.lnc:+.3f} (rank {int((P.lnc>=tr.lnc).sum())}/66)   "
          f"GC {tr.gc:+.3f} (rank {int((P.gc.abs()>=abs(tr.gc)).sum())}/66)")

print("== bin-median normalisation instead of splines, several filters ==")
run_variant(5, 8, "min count 5, 8x8 bins")
run_variant(5, 12, "min count 5, 12x12 bins")
run_variant(10, 10, "min count 10, 10x10 bins")
run_variant(2, 10, "min count 2, 10x10 bins")

print("\n== no GC/length correction at all (raw CPM), same permutation ==")
ok = c.min(axis=1) >= 5
cp = (c / c.sum() * 1e6)[ok]; Ks = K[ok]
y = np.log2(cp + 1)
near = ((Ks.seq == "chr11") & (((Ks.start + Ks.end) / 2 - NPOS).abs() < 2e6)).values
lnc = (Ks.gene_type == "lncRNA").values
out = []
for pair in itertools.combinations(samples, 2):
    rest = [s for s in samples if s not in pair]
    d = y[rest].mean(axis=1).values - y[list(pair)].mean(axis=1).values
    out.append((set(pair) == set(negs), d[near].mean() - d[~near].mean(), d[lnc].mean() - d[~lnc].mean()))
P = pd.DataFrame(out, columns=["t", "cis", "lnc"])
tr = P[P.t].iloc[0]
print(f"  raw: cis {tr.cis:+.3f} rank {int((P.cis>=tr.cis).sum())}/66 ; lnc {tr.lnc:+.3f} rank {int((P.lnc>=tr.lnc).sum())}/66")

print("\n== cis effect computed from the published DE tables directly (no reprocessing) ==")
for f in ["80", "40even", "40odd", "40old", "20old"]:
    d = pd.read_csv(f"data/DE_condition_{f}.csv").set_index("Geneid")
    j = M.join(d[["baseMean", "log2FoldChange"]], rsuffix="_t")
    j = j[(j.baseMean_t > 10) & j.log2FoldChange_t.notna() & j.on_primary_chrom & (j.gene_name != "NEAT1")]
    nr = (j.seq == "chr11") & (((j.start + j.end) / 2 - NPOS).abs() < 2e6)
    print(f"  {f:7s} n_near={nr.sum():3d}  median LFC near {j.log2FoldChange_t[nr].median():+.3f} vs "
          f"rest {j.log2FoldChange_t[~nr].median():+.3f}  MWU p={stats.mannwhitneyu(j.log2FoldChange_t[nr], j.log2FoldChange_t[~nr], alternative='greater').pvalue:.3f}")
