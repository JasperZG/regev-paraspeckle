"""Re-test every feature against the CQN-corrected response."""
import pandas as pd, numpy as np
from scipy import stats

M = pd.read_pickle("master.pkl")
res = pd.read_pickle("cqn_response.pkl")
K = M.loc[res.index].copy()
K["y"] = res.y_cqn
K["yA"] = res.y_repA; K["yB"] = res.y_repB
K["intron_frac"] = K.len_intron / K.len_gene
K["n_introns"] = (K.n_exons - 1).clip(lower=0)
K["is_lnc"] = (K.gene_type == "lncRNA").astype(int)

feats = [c for c in K.columns if c.startswith(("gc_", "log10_len", "gene_frac", "exon_frac", "utr3_frac", "utr5_frac",
                                               "cds_frac", "gene_nAlu", "exon_nAlu", "utr3_nAlu", "gene_IRAlu",
                                               "exon_IRAlu", "utr3_IRAlu", "gene_has", "exon_has", "utr3_has"))]
feats += ["n_exons", "n_transcripts", "intron_frac", "n_introns", "log10_baseMean", "is_lnc", "log2_nuc_cyt", "level"]
rows = []
for c in feats:
    s = K[[c, "y", "yA", "yB"]].dropna()
    if len(s) < 200 or s[c].std() == 0:
        continue
    rows.append(dict(feature=c, n=len(s), rho=stats.spearmanr(s[c], s.y)[0],
                     rhoA=stats.spearmanr(s[c], s.yA)[0], rhoB=stats.spearmanr(s[c], s.yB)[0]))
T = pd.DataFrame(rows).set_index("feature")
T["min_abs_rep"] = T[["rhoA", "rhoB"]].abs().min(axis=1)
T = T.reindex(T.rho.abs().sort_values(ascending=False).index)
print(T.round(3).head(28).to_string())
T.to_csv("res_cqn_scan.tsv", sep="\t")

print("\n== gene_type ==")
K["bt"] = K.gene_type.where(K.gene_type.isin(["protein_coding", "lncRNA"]), "other")
print(K.groupby("bt").agg(n=("y", "size"), med=("y", "median"), mean=("y", "mean")).round(3).to_string())
print("MWU lnc vs pc p =", stats.mannwhitneyu(K.y[K.bt == "lncRNA"], K.y[K.bt == "protein_coding"]).pvalue)

print("\n== chromosome ==")
cm = K.groupby("seq").agg(n=("y", "size"), med=("y", "median")).sort_values("med")
print(cm.round(3).T.to_string())

print("\n== distance to NEAT1 on chr11 ==")
neat = M[M.gene_name == "NEAT1"].iloc[0]; npos = (neat.start + neat.end) / 2
c11 = K[K.seq == "chr11"].copy(); c11["d"] = ((c11.start + c11.end) / 2 - npos).abs() / 1e6
for lo, hi in [(0, 0.5), (0.5, 2), (2, 5), (5, 20), (20, 200)]:
    s = c11[(c11.d >= lo) & (c11.d < hi)]
    print(f"  {lo}-{hi} Mb n={len(s):4d} median y={s.y.median():+.3f}")
print("  all genes median y = %+.3f" % K.y.median())
