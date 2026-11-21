"""The MALAT1 anchor as a methodological control for the NEAT1 anchor.

The folder contains three MALAT1 rtRT-seq pulldowns (NGS2 POS/NEG, NGS3 40probes
and 40block) and one SON antibody pulldown. All of the MALAT1 and SON libraries
are HFF, not HeLa: chrY is 617-946 CPM and COL1A2 1677-2801 CPM in them, against
7-10 and 0 in the HeLa NEAT1 libraries. So the MALAT1 arm is NOT a matched
control for the HeLa NEAT1 experiments; it is a control on the method.

Two documented label problems are corrected here:
  - NGS3 MALAT1 neg_Rep2 <-> 40probes_Rep1 are swapped in the count sheets
    (stated in the lab's MALAT1_coverage_NGS3/resolve_labels.R, and visible in
    the MALAT1 CPM: the "negative" carries 349k CPM).
  - NGS2 and NGS3 both have a library called MALAT1_neg/NEG, differing only in
    case, so the BAMs are kept in per-batch directories.

Counts come from featureCounts -s 1 on the same SAFs used for NEAT1:
  featureCounts -T 6 -s 1 -F SAF -a work/{gene_span,exon_union,intron}.saf \
      -o work/malat1_<kind>.txt data/rtRT_bam/ngs2/*.bam data/rtRT_bam/ngs3/*.bam
"""
import pandas as pd, numpy as np, itertools
from scipy import stats

NEAT1, MALAT1 = "ENSG00000245532", "ENSG00000251562"
NPOS = (65422774 + 65445540) / 2
SWAP = {"MALAT1_neg_Rep2": "MALAT1_40probes_Rep1", "MALAT1_40probes_Rep1": "MALAT1_neg_Rep2"}

feat = pd.read_csv("work/gene_features_gencode_v39.tsv", sep="\t", index_col=0)
feat.index = feat.index.str.split(".").str[0]; feat = feat[~feat.index.duplicated()]


def load(p, swap=False):
    d = pd.read_csv(p, sep="\t", index_col=0)
    d.index = d.index.astype(str).str.split(".").str[0]
    d = d.groupby(level=0).sum()
    return d.rename(columns=SWAP) if swap else d


def sizefactors(C):
    """DESeq median-of-ratios, robust to the anchor taking a third of the library."""
    lg = np.log(C.replace(0, np.nan))
    return np.exp((lg.sub(lg.mean(axis=1), axis=0)).median(axis=0))


def lfc(C, pos, neg):
    C = C.loc[(C[pos + neg] > 0).all(axis=1)]
    n = C[pos + neg].div(sizefactors(C[pos + neg]), axis=1)
    return np.log2(n[pos].mean(axis=1) / n[neg].mean(axis=1))


def spliceosomal(idx):
    nm = feat.gene_name.reindex(idx)
    return np.array([isinstance(x, str) and not x.endswith("P")
                     and x.startswith(("RNU1-", "RNVU1-", "RNU2-", "RNU4-", "RNU5", "RNU11", "RNU12", "RNU6ATAC"))
                     for x in nm])


def snrna_test(tag, L):
    f = feat.reindex(L.index); sp = spliceosomal(L.index); pc = (f.gene_type == "protein_coding").values
    p = stats.mannwhitneyu(L[sp], L[pc]).pvalue if sp.sum() >= 4 else np.nan
    print(f"  {tag:46s} anchor {L.get(MALAT1, L.get(NEAT1)):+.2f}  "
          f"snRNA n={sp.sum():2d} med {L[sp].median():+.3f}  mRNA {L[pc].median():+.3f}  p={p:.1e}")


print("=== cell line of every library (chrY CPM, COL1A2 CPM)")
gm = load("work/counts_malat1_gene_span.tsv", swap=True)
C4 = load("work/counts_genespan.tsv")
Y = [g for g in feat.index[feat.seq == "chrY"]]
col1a2 = feat.index[feat.gene_name == "COL1A2"][0]
for tag, C in [("MALAT1/SON arm", gm), ("NEAT1 arm (NGS4)", C4[["neg_Rep1", "80_Rep1"]])]:
    cpm = C.div(C.sum(), axis=1) * 1e6
    print(f"  {tag:18s} chrY {cpm.reindex([g for g in Y if g in cpm.index]).sum().mean():7.1f}   "
          f"COL1A2 {cpm.loc[col1a2].mean() if col1a2 in cpm.index else 0:7.0f}")

print("\n=== spliceosomal snRNA recovery, each pulldown against its own no-probe control")
M = {}
for tag, pos, neg in [
        ("NGS2 MALAT1 POS (HFF)", ["MALAT1_POS_Rep1", "MALAT1_POS_Rep2"], ["MALAT1_NEG_Rep1", "MALAT1_NEG_Rep2"]),
        ("NGS3 MALAT1 40probes (HFF)", ["MALAT1_40probes_Rep1", "MALAT1_40probes_Rep2"], ["MALAT1_neg_Rep1", "MALAT1_neg_Rep2"]),
        ("NGS3 MALAT1 40block (HFF)", ["MALAT1_40block_Rep1", "MALAT1_40block_Rep2"], ["MALAT1_neg_Rep1", "MALAT1_neg_Rep2"])]:
    M[tag] = lfc(gm, pos, neg); snrna_test(tag, M[tag])
snrna_test("NGS4 NEAT1, all 10 positives (HeLa)", lfc(C4, [c for c in C4.columns if not c.startswith("neg_")], ["neg_Rep1", "neg_Rep2"]))
C6 = load("work/counts_ngs6_genespan.tsv")
n6 = [c for c in C6.columns if c.split("_REP")[0] in ("NEG", "N3O", "N3E", "N3T", "NM", "N80")]
snrna_test("NGS6 NEAT1, all positives (HeLa)", lfc(C6, [c for c in C6.columns if c not in n6], n6))

print("\n=== genes above 2-fold in all three MALAT1 pulldowns")
D = pd.DataFrame(M).dropna(); hit = D[D.min(axis=1) > 1]
print(f"  {len(hit)} genes; by class: "
      + ", ".join(f"{k} {v}" for k, v in feat.gene_type.reindex(hit.index).value_counts().items()))
print("  " + ", ".join(feat.gene_name.reindex(hit.index).sort_values().dropna().astype(str)))
print(f"  pairwise r between the three: {D.corr().values[np.triu_indices(3,1)].round(3)}")

print("\n=== the locus-proximity halo, both anchors, GC/length corrected, intron vs exon")


def bs(x, df=5):
    x = np.asarray(x, float); q = np.quantile(x, np.linspace(0, 1, df + 1)[1:-1])
    return np.column_stack([np.ones_like(x), x, x ** 2, x ** 3] + [np.clip(x - k, 0, None) ** 3 for k in q])


def halo(tag, gsf, exf, itf, pos, neg, swap=False):
    gs, ex, it = load(gsf, swap), load(exf, swap), load(itf, swap)
    lib = gs.sum(); f = feat.reindex(ex.index)
    keep = (f.on_primary_chrom & (f.seq != "chrM")).fillna(False)
    keep &= (ex[neg].sum(axis=1) >= 20) & (it[neg].sum(axis=1) >= 20) & (ex[pos].sum(axis=1) >= 20) & (it[pos].sum(axis=1) >= 20)
    keep &= f.gc_exon.notna() & f.gc_intron.notna() & (f.len_intron > 200)
    K = f[keep]; out = {}
    for nm, cnt, gc, ll in [("E", ex, K.gc_exon.values, K.log10_len_exon.values),
                            ("I", it, K.gc_intron.values, K.log10_len_intron.values)]:
        y = np.log2(cnt.loc[K.index].div(lib, axis=1) * 1e6 + 0.5)
        X = np.column_stack([bs(gc), bs(ll)[:, 1:]])
        R = pd.DataFrame(index=y.index, columns=y.columns, dtype=float)
        for s in y.columns:
            b, *_ = np.linalg.lstsq(X, y[s].values, rcond=None); R[s] = y[s].values - X @ b
        out[nm] = R[pos].mean(axis=1) - R[neg].mean(axis=1)
    mid = (K.start + K.end) / 2
    parts = []
    for r in [0.5, 1, 2, 5]:
        m = ((K.seq == "chr11") & ((mid - NPOS).abs() < r * 1e6)).values
        parts.append(f"{r}Mb I{out['I'][m].mean()-out['I'].mean():+.2f}/E{out['E'][m].mean()-out['E'].mean():+.2f}")
    print(f"  {tag:34s} " + "  ".join(parts))


halo("NEAT1 anchor NGS4 (HeLa)", "work/counts_genespan.tsv", "work/counts_exon.tsv", "work/counts_intron.tsv",
     [c for c in C4.columns if not c.startswith("neg_")], ["neg_Rep1", "neg_Rep2"])
for tag, pos, neg in [("MALAT1 NGS2 (HFF)", ["MALAT1_POS_Rep1", "MALAT1_POS_Rep2"], ["MALAT1_NEG_Rep1", "MALAT1_NEG_Rep2"]),
                      ("MALAT1 NGS3 40block (HFF)", ["MALAT1_40block_Rep1", "MALAT1_40block_Rep2"], ["MALAT1_neg_Rep1", "MALAT1_neg_Rep2"])]:
    halo(tag, "work/counts_malat1_gene_span.tsv", "work/counts_malat1_exon.tsv", "work/counts_malat1_intron.tsv", pos, neg, swap=True)

print("\n=== cross-anchor: does either anchor capture the other RNA?")
print("  NEAT1 in the MALAT1 pulldowns: " + ", ".join(f"{k.split()[1]} {v.get(NEAT1):+.2f}" for k, v in M.items()))
L4 = lfc(C4, [c for c in C4.columns if not c.startswith("neg_")], ["neg_Rep1", "neg_Rep2"])
print(f"  MALAT1 in the NEAT1 pulldown (NGS4): {L4.get(MALAT1):+.2f}")
