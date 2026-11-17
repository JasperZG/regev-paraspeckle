"""Batch-aware response, and cis-proximity to the NEAT1 locus.

Library PCA puts {neg, 20old, 40odd} on one side and {40old, 40even, 80} on the
other. The neg controls sit in the first group, so 20old/40odd vs neg is a
within-batch contrast and 40old/40even/80 vs neg is confounded with that axis.
"""
import pandas as pd, numpy as np
from scipy import stats

M = pd.read_pickle("master.pkl")
res = pd.read_pickle("cqn_response.pkl")
K = M.loc[res.index].copy()
K["y_within"] = res[["y_20old", "y_40odd"]].mean(axis=1)      # within-batch
K["y_across"] = res[["y_40old", "y_40even", "y_80"]].mean(axis=1)
K["y_all"] = res.y_cqn
K["mid"] = (K.start + K.end) / 2
neat = M[M.gene_name == "NEAT1"].iloc[0]
NPOS = (neat.start + neat.end) / 2
print("NEAT1 locus chr11:%.2f Mb" % (NPOS / 1e6))

print("\n== chr11 enrichment of the reproducible set ==")
groups = pd.read_pickle("probe_groups.pkl")
bg = (K.seq == "chr11").mean()
for lab in ["both_halves", "even_only", "odd_only", "flat_control"]:
    s = K.reindex(groups[lab])
    f = (s.seq == "chr11").mean()
    k = (s.seq == "chr11").sum()
    p = stats.binomtest(int(k), len(s), bg, alternative="greater").pvalue
    print(f"  {lab:14s} chr11 {k:3d}/{len(s):4d} = {f:.3f}  (background {bg:.3f})  p={p:.2g}")

print("\n== distance to NEAT1, within-batch response ==")
c11 = K[K.seq == "chr11"].copy()
c11["d"] = (c11.mid - NPOS).abs() / 1e6
bins = [0, 0.5, 1, 2, 4, 8, 16, 32, 200]
for lo, hi in zip(bins[:-1], bins[1:]):
    s = c11[(c11.d >= lo) & (c11.d < hi)]
    if len(s) == 0: continue
    print(f"  {lo:5.1f}-{hi:<5.1f} Mb  n={len(s):4d}  within={s.y_within.median():+.3f}  "
          f"across={s.y_across.median():+.3f}  all={s.y_all.median():+.3f}")
print(f"  genome-wide median           within={K.y_within.median():+.3f}  "
      f"across={K.y_across.median():+.3f}  all={K.y_all.median():+.3f}")

near = c11[c11.d < 2]; far = K[~((K.seq == "chr11") & ((K.mid - NPOS).abs() < 2e6))]
for col in ["y_within", "y_across", "y_all"]:
    u = stats.mannwhitneyu(near[col], far[col], alternative="greater")
    print(f"  <2 Mb vs rest, {col}: n={len(near)} median {near[col].median():+.3f} vs {far[col].median():+.3f}  p={u.pvalue:.2g}")

print("\n== is it chr11 as a whole, or the NEAT1 neighbourhood? ==")
for col in ["y_within", "y_all"]:
    print(f"  {col}: chr11 median {c11[col].median():+.3f} vs other chroms {K[K.seq!='chr11'][col].median():+.3f}; "
          f"rho(distance, y) within chr11 = {stats.spearmanr(c11.d, c11[col])[0]:+.3f}")

print("\n== same test around MALAT1 (adjacent locus, different body) ==")
mal = M[M.gene_name == "MALAT1"].iloc[0]; MPOS = (mal.start + mal.end) / 2
c11["dm"] = (c11.mid - MPOS).abs() / 1e6
for lo, hi in [(0, 0.5), (0.5, 2), (2, 8), (8, 200)]:
    s = c11[(c11.dm >= lo) & (c11.dm < hi)]
    print(f"  {lo}-{hi} Mb from MALAT1 n={len(s):4d} within={s.y_within.median():+.3f}")

print("\n== top genes within 5 Mb of NEAT1 ==")
t = c11[c11.d < 5].sort_values("y_all", ascending=False).head(18)
print(t[["gene_name", "gene_type", "d", "baseMean", "y_within", "y_across", "y_all", "lfc"]].round(3).to_string(index=False))
K.to_pickle("master_resp.pkl")
