"""Characterise the genes that replicate across both batches."""
import pandas as pd, numpy as np, bioframe as bf
from scipy import stats

M = pd.read_pickle("master.pkl")
res = pd.read_pickle("cqn_response.pkl")
K = M.loc[res.index].copy()
K["y1"] = res[["y_20old", "y_40odd"]].mean(axis=1)
K["y2"] = res[["y_40old", "y_40even", "y_80"]].mean(axis=1)
K["y"] = res.y_cqn
K["mid"] = (K.start + K.end) / 2
z1 = (K.y1 - K.y1.mean()) / K.y1.std(); z2 = (K.y2 - K.y2.mean()) / K.y2.std()
K["zmin"] = np.minimum(z1, z2)

print("response distribution: SD(y)=%.3f  1%%/99%% = %.2f / %.2f" % (K.y.std(), K.y.quantile(.01), K.y.quantile(.99)))
print("genes with zmin > 2: %d ; zmin < -2: %d" % ((K.zmin > 2).sum(), (np.maximum(z1, z2) < -2).sum()))

top = K[K.zmin > 1.5].sort_values("zmin", ascending=False)
bot = K[np.maximum(z1, z2) < -1.5]
print("\nTOP (replicating in both batches), n=%d" % len(top))
print(", ".join(str(x) for x in top.gene_name.head(60)))
print("\nBOTTOM (depleted in both), n=%d" % len(bot))
print(", ".join(str(x) for x in bot.sort_values("y").gene_name.head(45)))

print("\n== composition ==")
for lab, s in [("top", top), ("bottom", bot), ("all", K)]:
    print(f"  {lab:7s} n={len(s):5d}  lncRNA={np.mean(s.gene_type=='lncRNA'):.2f}  "
          f"AS/DT name={np.mean(s.gene_name.astype(str).str.contains('-AS|-DT|^LINC')):.2f}  "
          f"median baseMean={s.baseMean.median():7.1f}  median GC={s.gc_gene.median():.3f}  "
          f"median log10len={s.log10_len_gene.median():.2f}  chr11={np.mean(s.seq=='chr11'):.3f}")

# does a gene overlap another gene on the opposite strand? (unstranded-counting bleed)
g = K.reset_index()[["seq", "start", "end", "Geneid", "strand"]].rename(columns={"seq": "chrom"})
allg = M[M.on_primary_chrom & (M.baseMean.fillna(0) > 0)].reset_index()[["seq", "start", "end", "Geneid", "strand", "baseMean"]].rename(columns={"seq": "chrom"})
ov = bf.overlap(g, allg, how="inner", suffixes=("", "_o"))
ov = ov[ov.Geneid != ov.Geneid_o]
anti = ov[ov.strand != ov.strand_o].groupby("Geneid").baseMean_o.max()
K["anti_partner_bm"] = anti.reindex(K.index).fillna(0)
K["has_anti"] = (K.anti_partner_bm > 0).astype(int)
K["anti_ratio"] = np.log2((K.anti_partner_bm + 1) / (K.baseMean + 1))
print("\n== antisense-overlap (unstranded counting bleed) ==")
print("  genes with an opposite-strand overlapping gene: %.2f" % K.has_anti.mean())
for lab, s in [("top", K.loc[top.index]), ("bottom", K.loc[bot.index]), ("all", K)]:
    print(f"  {lab:7s} has_anti={s.has_anti.mean():.2f}  median anti_ratio={s.anti_ratio.median():+.2f}")
print("  rho(y, anti_ratio) = %+.3f ; rho(y, has_anti) = %+.3f"
      % (stats.spearmanr(K.anti_ratio, K.y)[0], stats.spearmanr(K.has_anti, K.y)[0]))

print("\n== is the response just 'low expression' ? ==")
K["bmbin"] = pd.qcut(K.baseMean, 8)
print(K.groupby("bmbin", observed=True).agg(n=("y", "size"), med_y=("y", "median"), sd_y=("y", "std"),
                                            frac_top=("zmin", lambda v: (v > 1.5).mean())).round(3).to_string())
K.to_pickle("master_resp.pkl")
