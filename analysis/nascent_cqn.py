"""Intron vs exon enrichment after per-library GC/length correction.

The delivered counts are whole-gene-span counts, so they mix intronic and exonic
reads. Here exon-union and intron counts (featureCounts -s 1 on the SAFs in work/)
are corrected separately: for each library and each region type, fit
log2 CPM ~ s(GC of that region) + s(log length of that region) and keep the
residual. The response for a region is mean(residual over probe libraries) -
mean(residual over neg). Everything is then put through the 66-split label
permutation, which is the same bar used for the gene-level effects.
"""
import pandas as pd, numpy as np, itertools
from scipy import stats

feat = pd.read_csv("work/gene_features_gencode_v39.tsv", sep="\t", index_col=0)
ex = pd.read_csv("work/counts_exon.tsv", sep="\t", index_col=0)
it = pd.read_csv("work/counts_intron.tsv", sep="\t", index_col=0)
gs = pd.read_csv("work/counts_genespan.tsv", sep="\t", index_col=0)
lib = gs.sum()
negs = ["neg_Rep1", "neg_Rep2"]; samples = list(gs.columns); probes = [s for s in samples if s not in negs]
NPOS = (65422774 + 65445540) / 2

keep = (feat.on_primary_chrom & (feat.seq != "chrM") & (feat.gene_name != "NEAT1")).reindex(ex.index).fillna(False)
keep &= (ex[negs].sum(axis=1) >= 20) & (it[negs].sum(axis=1) >= 20) & (ex[probes].sum(axis=1) >= 50) & (it[probes].sum(axis=1) >= 50)
keep &= feat.gc_exon.notna() & feat.gc_intron.notna() & (feat.len_intron > 200)
K = feat[keep].copy()
print("genes:", len(K))

def bs(x, df=5):
    x = np.asarray(x, float); q = np.quantile(x, np.linspace(0, 1, df + 1)[1:-1])
    return np.column_stack([np.ones_like(x), x, x ** 2, x ** 3] + [np.clip(x - k, 0, None) ** 3 for k in q])

def resid(cnt, gc, loglen):
    y = np.log2(cnt.loc[K.index].div(lib, axis=1) * 1e6 + 0.5)
    X = np.column_stack([bs(gc), bs(loglen)[:, 1:]])
    R = pd.DataFrame(index=y.index, columns=y.columns, dtype=float)
    for s in y.columns:
        b, *_ = np.linalg.lstsq(X, y[s].values, rcond=None); R[s] = y[s].values - X @ b
    return R

RE = resid(ex, K.gc_exon.values, K.log10_len_exon.values)
RI = resid(it, K.gc_intron.values, K.log10_len_intron.values)
for s in [negs[0], probes[0], probes[-1]]:
    print(f"  residual GC assoc {s}: exon rho={stats.spearmanr(K.gc_exon, RE[s])[0]:+.3f} intron rho={stats.spearmanr(K.gc_intron, RI[s])[0]:+.3f}")
K["yE"] = RE[probes].mean(axis=1) - RE[negs].mean(axis=1)
K["yI"] = RI[probes].mean(axis=1) - RI[negs].mean(axis=1)
K["dIE"] = K.yI - K.yE
K["mid"] = (K.start + K.end) / 2
K["near"] = (K.seq == "chr11") & ((K.mid - NPOS).abs() < 2e6)
K["lnc"] = K.gene_type == "lncRNA"
best = pd.read_pickle("work/probe_pattern_exon.pkl")
spec = set(best[(best.b > 1.0) & (best.r2 > 0.8) & (best["max"] > 1.0)].index)
K["pattern"] = K.index.isin(spec)
alu = pd.read_pickle("work/alu_counts.pkl").reindex(K.index).fillna(0)
K["sense_alu"] = alu.n_sense.values  # sense-oriented Alu elements in exons (JXJ 572 cross-hybridisation substrate)
K["noalu"] = (K.sense_alu == 0) & ~K.pattern
print("genes with no sense exonic Alu and not a probe-signature gene:", int(K.noalu.sum()))
print("\ncorrected responses (median): exon %+.3f  intron %+.3f  dIE %+.3f ; frac dIE>0 %.3f" % (K.yE.median(), K.yI.median(), K.dIE.median(), (K.dIE > 0).mean()))
print("residual association of dIE with GC_exon-GC_intron: rho=%+.3f ; with log len intron: %+.3f" % (stats.spearmanr(K.gc_exon - K.gc_intron, K.dIE)[0], stats.spearmanr(K.log10_len_intron, K.dIE)[0]))
for lab, m in [("near NEAT1 <2Mb", K.near), ("lncRNA", K.lnc), ("lncRNA, Alu-free", K.lnc & K.noalu), ("protein_coding", K.gene_type == "protein_coding"), ("protein_coding, Alu-free", (K.gene_type == "protein_coding") & K.noalu), ("probe-pattern genes", K.pattern), ("sense exonic Alu >=2", K.sense_alu >= 2)]:
    s = K[m]; print(f"  {lab:20s} n={len(s):5d} yE {s.yE.median():+.3f}  yI {s.yI.median():+.3f}  dIE {s.dIE.median():+.3f}")

# reproducibility across independent library sets
b1 = ["20old_Rep1", "20old_Rep2", "40odd_Rep1", "40odd_Rep2"]; b2 = ["40old_Rep1", "40old_Rep2", "40even_Rep1", "40even_Rep2", "80_Rep1", "80_Rep2"]
d1 = (RI[b1].mean(axis=1) - RE[b1].mean(axis=1)) - (RI[negs].mean(axis=1) - RE[negs].mean(axis=1))
d2 = (RI[b2].mean(axis=1) - RE[b2].mean(axis=1)) - (RI[negs].mean(axis=1) - RE[negs].mean(axis=1))
print("dIE reproducibility batch1 vs batch2: r=%.3f" % np.corrcoef(d1, d2)[0, 1])
print("yI reproducibility: r=%.3f ; yE: r=%.3f" % (np.corrcoef(RI[b1].mean(axis=1) - RI[negs].mean(axis=1), RI[b2].mean(axis=1) - RI[negs].mean(axis=1))[0, 1],
                                                    np.corrcoef(RE[b1].mean(axis=1) - RE[negs].mean(axis=1), RE[b2].mean(axis=1) - RE[negs].mean(axis=1))[0, 1]))

# 66-split permutation
near = K.near.values; lnc = K.lnc.values; pat = K.pattern.values; noalu = K.noalu.values; pc = (K.gene_type == "protein_coding").values
rows = []
for pair in itertools.combinations(samples, 2):
    rest = [s for s in samples if s not in pair]
    yE = (RE[rest].mean(axis=1) - RE[list(pair)].mean(axis=1)).values
    yI = (RI[rest].mean(axis=1) - RI[list(pair)].mean(axis=1)).values
    d = yI - yE
    rows.append(dict(true=set(pair) == set(negs), med_dIE=np.median(d), frac_pos=(d > 0).mean(),
                     near_I=yI[near].mean() - yI[~near].mean(), near_E=yE[near].mean() - yE[~near].mean(), near_d=d[near].mean() - d[~near].mean(),
                     lnc_I=yI[lnc & ~pat].mean() - yI[~lnc & ~pat].mean(), lnc_E=yE[lnc & ~pat].mean() - yE[~lnc & ~pat].mean(),
                     lncA_E=yE[lnc & noalu].mean() - yE[pc & noalu].mean(), lncA_I=yI[lnc & noalu].mean() - yI[pc & noalu].mean(),
                     nearA_I=yI[near & noalu].mean() - yI[~near & noalu].mean(), nearA_E=yE[near & noalu].mean() - yE[~near & noalu].mean(),
                     alu_E=yE[K.sense_alu.values >= 2].mean() - yE[K.sense_alu.values == 0].mean()))
P = pd.DataFrame(rows); t = P[P["true"]].iloc[0]; P.to_pickle("work/nascent_perm66.pkl")
print("\n66-split permutation (rank 1 = most extreme of 66; floor p = 0.015):")
for c in ["med_dIE", "frac_pos", "near_I", "near_E", "near_d", "lnc_I", "lnc_E", "lncA_E", "lncA_I", "nearA_I", "nearA_E", "alu_E"]:
    print(f"  {c:9s} true {t[c]:+.3f}  rank {(P[c] >= t[c]).sum()}/66  null [{P[c].min():+.3f}, {P[c].max():+.3f}]")

# cis decay with the intron response vs the exon response
print("\ncis decay (mean response, genes within r Mb of NEAT1 on chr11):")
for r in [0.5, 1, 2, 3, 5, 10]:
    m = (K.seq == "chr11") & ((K.mid - NPOS).abs() < r * 1e6)
    print(f"  r={r:4.1f} Mb n={m.sum():4d}  intron {K.yI[m].mean():+.3f}  exon {K.yE[m].mean():+.3f}  (genome: intron {K.yI.mean():+.3f} exon {K.yE.mean():+.3f})")

# 2 Mb sliding-window null for the intron response
stat = []
for ch, g in K.groupby("seq"):
    p = g.mid.values; v = g.yI.values
    for cc in np.arange(p.min(), p.max(), 1e6):
        m = np.abs(p - cc) < 2e6
        if 40 <= m.sum() <= 160:
            rest = np.concatenate([K.yI.values[K.seq.values != ch], v[~m]])
            stat.append((ch, cc, v[m].mean() - rest.mean()))
S = pd.DataFrame(stat, columns=["chrom", "centre", "eff"])
obs = K.yI[near].mean() - K.yI[~near].mean()
print(f"\n2 Mb windows (40-160 genes): n={len(S)}; NEAT1 window intron effect {obs:+.3f}; windows beating it: {(S.eff >= obs).sum()}")
print(S.sort_values("eff", ascending=False).head(6).to_string(index=False))

# Yap et al. 2022 HyPro-seq targets
Y = pd.read_csv("ext/yap2022/neat1_hypro_targets.tsv", sep="\t")
yens = set(Y.gene_id.astype(str).str.split(".").str[0])
K["hypro"] = K.index.str.split(".").str[0].isin(yens)
for col in ["yI", "yE", "dIE"]:
    a = K.loc[K.hypro, col]; b = K.loc[~K.hypro, col]
    print(f"HyPro targets (n={len(a)}) {col}: median {a.median():+.3f} vs others {b.median():+.3f}  MWU p={stats.mannwhitneyu(a, b).pvalue:.2e}")
a = K[K.hypro & (K.seq != "chr11")]; b = K[~K.hypro & (K.seq != "chr11")]
print(f"HyPro targets off chr11 (n={len(a)}): yI {a.yI.median():+.3f} vs {b.yI.median():+.3f} p={stats.mannwhitneyu(a.yI, b.yI).pvalue:.2e}; yE {a.yE.median():+.3f} vs {b.yE.median():+.3f}")
K.to_pickle("work/nascent_cqn.pkl")
RE.to_pickle("work/resid_exon.pkl"); RI.to_pickle("work/resid_intron.pkl")
