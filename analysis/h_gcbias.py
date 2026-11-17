"""Is the negative GC relationship technical (library/PCR GC bias) or biological?

Tests
  1. per-sample GC bias profile: rho(GC, log ratio of that sample to the mean of all)
  2. null contrast neg_Rep1 vs neg_Rep2 (no biology, same protocol) -> does GC still predict?
  3. pulldown vs pulldown contrasts (same input class)
  4. dose-response across probe sets 20old/40old/40even/40odd/80
"""
import pandas as pd, numpy as np
from scipy import stats

M = pd.read_pickle("master.pkl")
cnt = pd.read_csv("data/counts.tsv", sep="\t", index_col=0)
K = M[M.keep]
c = cnt.loc[K.index]
cpm = c / c.sum() * 1e6
ok = (c.min(axis=1) >= 5)
L = np.log2(cpm[ok] + 1)
gc = K.gc_gene[ok]
ln = K.log10_len_gene[ok]
print("genes used for bias profile:", ok.sum())

print("\n== 1. per-sample GC bias (log2 sample / mean of all samples) ==")
ctr = L.sub(L.mean(axis=1), axis=0)
for s in L.columns:
    print(f"  {s:20s} rho(GC) = {stats.spearmanr(gc, ctr[s])[0]:+.3f}   rho(len) = {stats.spearmanr(ln, ctr[s])[0]:+.3f}")

print("\n== 2. null and within-class contrasts (log2 ratio of summed CPM) ==")
def contrast(a, b, lab):
    x = np.log2((cpm.loc[ok, a].mean(axis=1) + 1) / (cpm.loc[ok, b].mean(axis=1) + 1))
    print(f"  {lab:42s} rho(GC)={stats.spearmanr(gc, x)[0]:+.3f}  rho(len)={stats.spearmanr(ln, x)[0]:+.3f}  SD={x.std():.3f}")
    return x
contrast(["NEAT1_neg_Rep1"], ["NEAT1_neg_Rep2"], "neg Rep1 vs neg Rep2 (null)")
contrast(["NEAT1_80_Rep1"], ["NEAT1_80_Rep2"], "80 Rep1 vs 80 Rep2 (null)")
contrast(["NEAT1_40even_Rep1"], ["NEAT1_40even_Rep2"], "40even Rep1 vs Rep2 (null)")
contrast(["NEAT1_40odd_Rep1", "NEAT1_40odd_Rep2"], ["NEAT1_40even_Rep1", "NEAT1_40even_Rep2"], "40odd vs 40even (probe halves)")
contrast(["NEAT1_80_Rep1", "NEAT1_80_Rep2"], ["NEAT1_40even_Rep1", "NEAT1_40even_Rep2"], "80 vs 40even (both pulldowns)")
probe = ["NEAT1_20old_Rep1","NEAT1_20old_Rep2","NEAT1_40old_Rep1","NEAT1_40old_Rep2","NEAT1_40odd_Rep1",
         "NEAT1_40odd_Rep2","NEAT1_40even_Rep1","NEAT1_40even_Rep2","NEAT1_80_Rep1","NEAT1_80_Rep2"]
contrast(probe, ["NEAT1_neg_Rep1", "NEAT1_neg_Rep2"], "all probe vs neg (the real contrast)")

print("\n== 3. dose-response: GC slope per probe set (vs neg) ==")
for lab, cols in [("20old", ["NEAT1_20old_Rep1","NEAT1_20old_Rep2"]), ("40old", ["NEAT1_40old_Rep1","NEAT1_40old_Rep2"]),
                  ("40odd", ["NEAT1_40odd_Rep1","NEAT1_40odd_Rep2"]), ("40even", ["NEAT1_40even_Rep1","NEAT1_40even_Rep2"]),
                  ("80", ["NEAT1_80_Rep1","NEAT1_80_Rep2"])]:
    x = np.log2((cpm.loc[ok, cols].mean(axis=1) + 1) / (cpm.loc[ok, ["NEAT1_neg_Rep1","NEAT1_neg_Rep2"]].mean(axis=1) + 1))
    neat = cpm.loc[M.index[M.gene_name == "NEAT1"], cols].mean(axis=1).iloc[0]
    print(f"  {lab:7s} NEAT1 CPM={neat:8.0f}  rho(GC)={stats.spearmanr(gc, x)[0]:+.3f}  rho(len)={stats.spearmanr(ln, x)[0]:+.3f}")

print("\n== 4. is the GC slope proportional to NEAT1 recovery? ==")
rows = []
for lab, cols in [("20old", ["NEAT1_20old_Rep1","NEAT1_20old_Rep2"]), ("40old", ["NEAT1_40old_Rep1","NEAT1_40old_Rep2"]),
                  ("40odd", ["NEAT1_40odd_Rep1","NEAT1_40odd_Rep2"]), ("40even", ["NEAT1_40even_Rep1","NEAT1_40even_Rep2"]),
                  ("80", ["NEAT1_80_Rep1","NEAT1_80_Rep2"])]:
    x = np.log2((cpm.loc[ok, cols].mean(axis=1) + 1) / (cpm.loc[ok, ["NEAT1_neg_Rep1","NEAT1_neg_Rep2"]].mean(axis=1) + 1))
    neat = cpm.loc[M.index[M.gene_name == "NEAT1"], cols].mean(axis=1).iloc[0]
    rows.append((lab, np.log2(neat / cpm.loc[M.index[M.gene_name == "NEAT1"], ["NEAT1_neg_Rep1","NEAT1_neg_Rep2"]].mean(axis=1).iloc[0]),
                 stats.spearmanr(gc, x)[0]))
d = pd.DataFrame(rows, columns=["set", "log2_NEAT1_enrich", "rho_GC"])
print(d.round(3).to_string(index=False))
print("  r(NEAT1 enrichment, GC slope) =", round(np.corrcoef(d.log2_NEAT1_enrich, d.rho_GC)[0, 1], 3))
