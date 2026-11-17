"""Are the three surviving effects independent, and how big are they together?"""
import pandas as pd, numpy as np, itertools
from scipy import stats
from sklearn.linear_model import LinearRegression

M = pd.read_pickle("master.pkl")
res = pd.read_pickle("cqn_response.pkl")
K = M.loc[res.index].copy()
K["y"] = res.y_cqn
K["y1"] = res[["y_20old", "y_40odd"]].mean(axis=1)
K["y2"] = res[["y_40old", "y_40even", "y_80"]].mean(axis=1)
K["mid"] = (K.start + K.end) / 2
neat = M[M.gene_name == "NEAT1"].iloc[0]; NPOS = (neat.start + neat.end) / 2
K["near_neat1"] = ((K.seq == "chr11") & ((K.mid - NPOS).abs() < 2e6)).astype(int)
K["is_lnc"] = (K.gene_type == "lncRNA").astype(int)
K["is_as"] = K.gene_name.astype(str).str.contains("-AS|-DT|^LINC").astype(int)
K["noncoding"] = (K.len_cds == 0).astype(int)

print("== overlap of the three ==")
print(pd.crosstab(K.is_lnc, K.is_as))
print("  near-NEAT1 genes that are lncRNA: %.2f (background %.2f)" % (K[K.near_neat1 == 1].is_lnc.mean(), K.is_lnc.mean()))

print("\n== each effect with the others held out ==")
def eff(mask_col, subset, lab):
    s = K[subset]
    a = s[s[mask_col] == 1].y; b = s[s[mask_col] == 0].y
    if len(a) < 10: return
    print(f"  {lab:52s} n={len(a):4d}  mean {a.mean():+.3f} vs {b.mean():+.3f}  "
          f"diff {a.mean()-b.mean():+.3f}  MWU p={stats.mannwhitneyu(a,b).pvalue:.1e}")
eff("is_lnc", K.index.notna(), "lncRNA (all genes)")
eff("is_lnc", K.near_neat1 == 0, "lncRNA (excluding the NEAT1 neighbourhood)")
eff("near_neat1", K.index.notna(), "<2 Mb of NEAT1 (all genes)")
eff("near_neat1", K.is_lnc == 0, "<2 Mb of NEAT1 (protein-coding only)")
eff("is_as", K.is_lnc == 0, "antisense/divergent name (non-lncRNA only)")
eff("is_as", K.is_lnc == 1, "antisense/divergent name (within lncRNA)")

print("\n== linear model, and leave-one-chromosome-out CV ==")
K["z_bm"] = (np.log10(K.baseMean) - np.log10(K.baseMean).mean()) / np.log10(K.baseMean).std()
blocks = {
    "lncRNA only": ["is_lnc"],
    "cis only": ["near_neat1"],
    "lncRNA + cis": ["is_lnc", "near_neat1"],
    "lncRNA + cis + antisense": ["is_lnc", "near_neat1", "is_as"],
    "+ abundance": ["is_lnc", "near_neat1", "is_as", "z_bm"],
}
yv = K.y.values
for lab, cols in blocks.items():
    X = K[cols].values.astype(float)
    pred = np.zeros(len(yv))
    for ch in K.seq.unique():
        te = (K.seq == ch).values
        m = LinearRegression().fit(X[~te], yv[~te]); pred[te] = m.predict(X[te])
    r2 = 1 - ((yv - pred) ** 2).sum() / ((yv - yv.mean()) ** 2).sum()
    co = LinearRegression().fit(X, yv)
    print(f"  {lab:26s} CV R2={r2:+.4f}   " + "  ".join(f"{c}={v:+.3f}" for c, v in zip(cols, co.coef_)))

print("\n== how much signal is there to explain? ==")
r12 = np.corrcoef(K.y1, K.y2)[0, 1]
print(f"  batch1 vs batch2 r = {r12:.3f} -> reliability of the pooled response (Spearman-Brown) = {2*r12/(1+r12):.3f}")
print(f"  so the ceiling on R2 for any gene feature is about {2*r12/(1+r12):.2f}; the model above reaches {max(0,r2):.4f}")

print("\n== divergent lncRNAs: does the response track the partner promoter? ==")
import bioframe as bf
lnc = K[(K.is_lnc == 1)].reset_index()[["seq", "start", "end", "Geneid", "strand", "y"]].rename(columns={"seq": "chrom"})
pc = M[(M.gene_type == "protein_coding") & M.baseMean.notna()].reset_index()
pc["tss"] = np.where(pc.strand == "+", pc.start, pc.end)
pcw = pc.assign(start=pc.tss - 2000, end=pc.tss + 2000)[["seq", "start", "end", "Geneid", "strand", "baseMean"]].rename(columns={"seq": "chrom"})
lncw = lnc.assign(tss=np.where(lnc.strand == "+", lnc.start, lnc.end))
lncw = lncw.assign(start=lncw.tss - 2000, end=lncw.tss + 2000)
ov = bf.overlap(lncw, pcw, how="left", suffixes=("", "_p"))
ov["div"] = ov.Geneid_p.notna() & (ov.strand != ov.strand_p)
agg = ov.groupby("Geneid").agg(div=("div", "max"), pbm=("baseMean_p", "max"))
print("  lncRNAs with a divergent coding partner (TSS within 2 kb, opposite strand): %d / %d" % (agg.div.sum(), len(agg)))
a = K.y.reindex(agg.index[agg.div.astype(bool)]); b = K.y.reindex(agg.index[~agg.div.astype(bool)])
print(f"  divergent lncRNA mean y={a.mean():+.3f} (n={len(a)}) vs other lncRNA {b.mean():+.3f} (n={len(b)})  p={stats.mannwhitneyu(a.dropna(),b.dropna()).pvalue:.3f}")
K.to_pickle("master_resp.pkl")
