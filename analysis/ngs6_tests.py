"""NGS6 as a prospective test of the second-pass conclusions.

NGS6 (data/rtRT/rtRT-seq/Feature count/NGS6_Featurecount.xlsx) has no JXJ 572, a
new 3' probe design (JXJ 1056-1091), five positive sets and six negative samples
(12 negative libraries). Only gene-span counts are available for it, no BAMs, so
there is no intron/exon split here.

Tests, all written down before this batch was read:
  1. the 86 probe-signature genes go flat
  2. the sense-Alu dose curve disappears
  3. TUB is up only in the sets carrying its JXJ 806 Region 6 probe
  4. GC fails once there are more than two negatives
  5. the cis intronic halo within 1-2 Mb
  6. what, if anything, is enriched in all five sets of both batches
"""
import pandas as pd, numpy as np
from scipy import stats

rng = np.random.default_rng(0)
NPOS = (65422774 + 65445540) / 2


def load(path):
    C = pd.read_csv(path, sep="\t", index_col=0)
    C.index = C.index.str.split(".").str[0]
    return C[~C.index.duplicated()]


def bs(x, df=5):
    x = np.asarray(x, float); q = np.quantile(x, np.linspace(0, 1, df + 1)[1:-1])
    return np.column_stack([np.ones_like(x), x, x ** 2, x ** 3] + [np.clip(x - k, 0, None) ** 3 for k in q])


def cqn(C, negs, poss, feat):
    """Per-library GC/length spline residuals; response = mean(pos) - mean(neg)."""
    lib = C.sum(); f = feat.reindex(C.index)
    keep = ((C[negs].sum(axis=1) >= 60) & (C[poss].sum(axis=1) >= 50)
            & f.seq.str.match(r"^chr[0-9X]+$").fillna(False) & f.gc_gene.notna() & (f.gene_name != "NEAT1"))
    K = f[keep]; y = np.log2(C[keep].div(lib, axis=1) * 1e6 + 0.5)
    X = np.column_stack([bs(K.gc_gene.values), bs(np.log10(K.len_gene.values))[:, 1:]])
    R = pd.DataFrame(index=y.index, columns=y.columns, dtype=float)
    for s in y.columns:
        b, *_ = np.linalg.lstsq(X, y[s].values, rcond=None); R[s] = y[s].values - X @ b
    return R, K


feat = load_feat = pd.read_csv("work/gene_features_gencode_v39.tsv", sep="\t", index_col=0)
feat.index = feat.index.str.split(".").str[0]; feat = feat[~feat.index.duplicated()]
alu = pd.read_pickle("work/alu_counts.pkl"); alu.index = alu.index.str.split(".").str[0]
alu = alu[~alu.index.duplicated()]

C6 = load("work/counts_ngs6_genespan.tsv")
NEGSETS = ("NEG", "N3O", "N3E", "N3T", "NM", "N80")
n6 = [c for c in C6.columns if c.split("_REP")[0] in NEGSETS]
p6 = {}
for c in C6.columns:
    if c not in n6:
        p6.setdefault(c.rsplit("_REP", 1)[0], []).append(c)
print(f"NGS6: {len(n6)} negative libraries, {sum(len(v) for v in p6.values())} positive libraries in {len(p6)} sets")

cpm6 = C6.div(C6.sum(), axis=1) * 1e6
negpool = cpm6[n6].mean(axis=1)
MATCH = {"3_ODD": "N3O", "3_EVEN": "N3E", "3_TOT": "N3T", "MED": "NM", "TOTAL": "N80"}
L6 = pd.DataFrame({k: np.log2((cpm6[v].mean(axis=1) + 0.5) / (negpool + 0.5)) for k, v in p6.items()})
L6m = pd.DataFrame({k: np.log2((cpm6[p6[k]].mean(axis=1) + 0.5)
                               / (cpm6[[c for c in C6.columns if c.startswith(MATCH[k] + "_")]].mean(axis=1) + 0.5))
                    for k in p6})

# --- 1/3. the probe-signature genes, and TUB
best = pd.read_pickle("work/probe_pattern_exon.pkl")
sig = [i.split(".")[0] for i in best[(best.b > 1.0) & (best.r2 > 0.8) & (best["max"] > 1.0)].index]
sub = L6.reindex([g for g in sig if g in L6.index])
print("\nthe NGS4 probe-signature genes in NGS6 (n=%d): median" % len(sub),
      {k: round(float(sub[k].median()), 2) for k in L6})
name2id = {v: k for k, v in feat.gene_name.items()}
for nm in ["TUB", "P3H4", "MTRFR", "SIRT7", "DHRS4-AS1", "MAP3K13", "TMEM9B-AS1", "H3-3A", "PSRC1", "MCU"]:
    g = name2id.get(nm)
    if g in L6.index:
        print(f"  {nm:12s}", {k: round(float(L6.loc[g, k]), 2) for k in L6})

# --- 2/4/5. corrected response, 10-vs-12 permutation
R6, K6 = cqn(C6, n6, [c for v in p6.values() for c in v], feat)
resp = R6[[c for v in p6.values() for c in v]].mean(axis=1) - R6[n6].mean(axis=1)
mid = (K6.start + K6.end) / 2
near = ((K6.seq == "chr11") & ((mid - NPOS).abs() < 2e6)).values
lnc = (K6.gene_type == "lncRNA").values
nalu = alu.n_sense.reindex(K6.index).fillna(0).values
true = dict(near=resp[near].mean() - resp[~near].mean(), lnc=resp[lnc].mean() - resp[~lnc].mean(),
            gc=stats.spearmanr(K6.gc_gene, resp)[0], alu=resp[nalu >= 2].mean() - resp[nalu == 0].mean())
cols = list(R6.columns); null = {k: [] for k in true}
for _ in range(2000):
    pp = list(rng.choice(cols, 10, replace=False)); nn = [c for c in cols if c not in pp]
    r = R6[pp].mean(axis=1) - R6[nn].mean(axis=1)
    null["near"].append(r[near].mean() - r[~near].mean()); null["lnc"].append(r[lnc].mean() - r[~lnc].mean())
    null["gc"].append(stats.spearmanr(K6.gc_gene, r)[0]); null["alu"].append(r[nalu >= 2].mean() - r[nalu == 0].mean())
print(f"\nNGS6 corrected response ({len(K6)} genes), 10 pos vs 12 neg, 2,000 random splits:")
for k, v in true.items():
    a = np.array(null[k]); print(f"  {k:5s} true {v:+.3f}  null {a.mean():+.3f} +- {a.std():.3f}  p {(a >= v).mean():.4f}")
d = (mid - NPOS).abs() / 1e6
print("cis decay:", {f"{r}Mb": round(float(resp[((K6.seq == 'chr11') & (d < r)).values].mean()), 3)
                     for r in [0.5, 1, 2, 3, 5, 10]}, "genome %+.3f" % resp.mean())
S = []
for ch, g in pd.DataFrame({"mid": mid, "v": resp, "chr": K6.seq}).groupby("chr"):
    pos = g.mid.values; val = g.v.values
    for cc in np.arange(pos.min(), pos.max(), 1e6):
        m = np.abs(pos - cc) < 2e6
        if 40 <= m.sum() <= 300: S.append(val[m].mean() - resp.mean())
S = np.array(S); obs = resp[near].mean() - resp.mean()
print(f"2 Mb windows: NEAT1 window {obs:+.3f}, rank {int((S >= obs).sum())}/{len(S)}")
nk = pd.read_pickle("work/nascent_cqn.pkl"); top = nk[nk.near].sort_values("yI", ascending=False).head(15)
sel = resp.reindex([i.split(".")[0] for i in top.index]).dropna()
print(f"the 15 strongest NGS4 cis-intron genes in NGS6: median {sel.median():+.3f} (genome {resp.median():+.3f})")

# --- 6. cross-batch replication
C4 = load("work/counts_genespan.tsv")
cpm4 = C4.div(C4.sum(), axis=1) * 1e6; neg4 = cpm4[["neg_Rep1", "neg_Rep2"]].mean(axis=1)
L4 = pd.DataFrame({s: np.log2((cpm4[[f"{s}_Rep1", f"{s}_Rep2"]].mean(axis=1) + 0.5) / (neg4 + 0.5))
                   for s in ["20old", "40old", "40odd", "40even", "80"]})
idx = L4.index.intersection(L6.index)
ok = (C6.loc[idx, n6].sum(axis=1) >= 60) & (C4.loc[idx, ["neg_Rep1", "neg_Rep2"]].sum(axis=1) >= 20)
A, B, Bm = L4.loc[idx][ok], L6.loc[idx][ok], L6m.loc[idx][ok]
R4, K4 = cqn(C4, ["neg_Rep1", "neg_Rep2"], [c for c in C4.columns if not c.startswith("neg_")], feat)
r4 = R4[[c for c in C4.columns if not c.startswith("neg_")]].mean(axis=1) - R4[["neg_Rep1", "neg_Rep2"]].mean(axis=1)
both = r4.index.intersection(resp.index)
print("\nNGS4 vs NGS6 corrected response: r = %+.3f genome-wide, %+.3f on chr11"
      % (np.corrcoef(r4[both], resp[both])[0, 1],
         np.corrcoef(r4[both][feat.seq.reindex(both) == "chr11"], resp[both][feat.seq.reindex(both) == "chr11"])[0, 1]))
ov = pd.read_pickle("work/probe_hits_genes.pkl"); ov["gid"] = ov.gene_id_g.str.split(".").str[0]
site = ov[ov.sense & (ov.mm <= 3)].groupby("gid").mm.min()
print("enriched in ALL five sets of BOTH batches (pooled and matched negatives):")
for thr in [1.0, 0.75, 0.5, 0.25]:
    m = (A.min(axis=1) > thr) & (B.min(axis=1) > thr) & (Bm.min(axis=1) > thr)
    names = list(feat.gene_name.reindex(A.index[m]).dropna())
    print(f"  log2 > {thr}: {int(m.sum())}" + (f"  {names}" if m.sum() <= 6 else ""))
m = (A.min(axis=1) > 0.5) & (B.min(axis=1) > 0.5) & (Bm.min(axis=1) > 0.5)
print(pd.DataFrame({"name": feat.gene_name.reindex(A.index[m]), "type": feat.gene_type.reindex(A.index[m]),
                    "NGS4_min": A[m].min(axis=1).round(2), "NGS6_min": B[m].min(axis=1).round(2),
                    "probe_site_mm": site.reindex(A.index[m])}).to_string())
L6.to_pickle("work/ngs6_log2_vs_neg.pkl"); L6m.to_pickle("work/ngs6_log2_vs_matched.pkl")
pd.DataFrame({"resp": resp, "near": near, "gene_name": K6.gene_name}).to_pickle("work/ngs6_cqn_response.pkl")
