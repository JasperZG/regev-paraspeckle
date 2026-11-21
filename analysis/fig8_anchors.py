"""Fig 8: the MALAT1 anchor as the control for the NEAT1 anchor.

(a) RNA-class enrichment for each anchor: the MALAT1 anchor recovers spliceosomal
    snRNAs, the NEAT1 anchor recovers nothing above background.
(b) the locus-proximity halo, intron and exon, for both anchors.
(c) out-of-sample variance of the log2 enrichment explained by gene properties.
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False})
feat = pd.read_csv("work/gene_features_gencode_v39.tsv", sep="\t", index_col=0)
feat.index = feat.index.str.split(".").str[0]; feat = feat[~feat.index.duplicated()]
NPOS = (65422774 + 65445540) / 2


def load(p):
    d = pd.read_csv(p, sep="\t", index_col=0); d.index = d.index.astype(str).str.split(".").str[0]
    return d.groupby(level=0).sum()


def bs(x, df=5):
    x = np.asarray(x, float); q = np.quantile(x, np.linspace(0, 1, df + 1)[1:-1])
    return np.column_stack([np.ones_like(x), x, x ** 2, x ** 3] + [np.clip(x - k, 0, None) ** 3 for k in q])


def anchor(gsf, exf, itf, pos, neg):
    gs, ex, it = load(gsf), load(exf), load(itf)
    lib = gs.sum(); f = feat.reindex(ex.index)
    keep = (f.on_primary_chrom & (f.seq != "chrM")).fillna(False)
    keep &= (ex[neg].sum(axis=1) >= 20) & (it[neg].sum(axis=1) >= 20) & (ex[pos].sum(axis=1) >= 20) & (it[pos].sum(axis=1) >= 20)
    keep &= f.gc_exon.notna() & f.gc_intron.notna() & (f.len_intron > 200)
    K = f[keep]
    out = {}
    for nm, cnt, gc, ll in [("E", ex, K.gc_exon.values, K.log10_len_exon.values),
                            ("I", it, K.gc_intron.values, K.log10_len_intron.values)]:
        y = np.log2(cnt.loc[K.index].div(lib, axis=1) * 1e6 + 0.5)
        X = np.column_stack([bs(gc), bs(ll)[:, 1:]])
        R = pd.DataFrame(index=y.index, columns=y.columns, dtype=float)
        for s in y.columns:
            b, *_ = np.linalg.lstsq(X, y[s].values, rcond=None); R[s] = y[s].values - X @ b
        out[nm] = R[pos].mean(axis=1) - R[neg].mean(axis=1)
    gcpm = gs.div(gs.sum(), axis=1) * 1e6
    ok = (gs[neg].sum(axis=1) >= 20) & (gs[pos].sum(axis=1) >= 20)
    out["G"] = np.log2((gcpm[pos].mean(axis=1) + 0.5) / (gcpm[neg].mean(axis=1) + 0.5))[ok]
    out["K"] = K
    return out


N4 = anchor("work/counts_genespan.tsv", "work/counts_exon.tsv", "work/counts_intron.tsv",
            ["20old_Rep1", "20old_Rep2", "40old_Rep1", "40old_Rep2", "40odd_Rep1", "40odd_Rep2",
             "40even_Rep1", "40even_Rep2", "80_Rep1", "80_Rep2"], ["neg_Rep1", "neg_Rep2"])
M2 = anchor("work/counts_malat1_gene_span.tsv", "work/counts_malat1_exon.tsv", "work/counts_malat1_intron.tsv",
            ["MALAT1_POS_Rep1", "MALAT1_POS_Rep2"], ["MALAT1_NEG_Rep1", "MALAT1_NEG_Rep2"])

fig = plt.figure(figsize=(7.2, 5.4))
gs_ = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.34, left=0.09, right=0.98, top=0.93, bottom=0.09)

# (a) RNA class
ax = fig.add_subplot(gs_[0, 0])
CLS = [("snRNA/scaRNA", ["snRNA", "scaRNA"]), ("lncRNA", ["lncRNA"]), ("mRNA", ["protein_coding"])]
w = 0.36
for i, (lab, L) in enumerate(CLS):
    for j, (nm, d, col) in enumerate([("NEAT1", N4, "0.15"), ("MALAT1", M2, "#c2452d")]):
        v = d["G"]; m = feat.gene_type.reindex(v.index).isin(L).values
        ax.bar(i + (j - 0.5) * w, v[m].median() - v.median(), w, color=col,
               label=nm + " anchor" if i == 0 else None, edgecolor="none")
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(range(len(CLS))); ax.set_xticklabels([c[0] for c in CLS], fontsize=7)
ax.set_ylabel("median log2 enrichment\nabove that library's median")
ax.legend(frameon=False, fontsize=6.5, loc="upper right")
ax.set_title("(a) which RNA classes come down", loc="left", fontsize=8)

# (b) halo
ax = fig.add_subplot(gs_[0, 1])
rr = np.array([0.25, 0.5, 1, 2, 3, 5, 10])
for nm, d, col in [("NEAT1", N4, "0.15"), ("MALAT1", M2, "#c2452d")]:
    K = d["K"]; mid = (K.start + K.end) / 2
    for key, ls in [("I", "-"), ("E", "--")]:
        y = [d[key][((K.seq == "chr11") & ((mid - NPOS).abs() < r * 1e6)).values].mean() - d[key].mean() for r in rr]
        ax.plot(rr, y, ls, color=col, lw=1.1, marker="o", ms=2.6,
                label=f"{nm}, {'intron' if key=='I' else 'exon'}")
ax.axhline(0, color="k", lw=0.6); ax.set_xscale("log")
ax.set_xlabel("genes within r Mb of the NEAT1/MALAT1 loci"); ax.set_ylabel("corrected log2 above genome")
ax.legend(frameon=False, fontsize=6, loc="upper right")
ax.set_title("(b) the locus-proximity halo follows both anchors", loc="left", fontsize=8)

# (c) snRNA per gene
ax = fig.add_subplot(gs_[1, 0])
sn = feat.index[feat.gene_type.isin(["snRNA", "scaRNA"])]
a = M2["G"].reindex(M2["G"].index.intersection(sn)).dropna().sort_values()
nmz = feat.gene_name.reindex(a.index).fillna("")
b = N4["G"].reindex(a.index)
yy = np.arange(len(a))
ax.barh(yy + 0.2, a.values, 0.4, color="#c2452d", label="MALAT1 anchor", edgecolor="none")
ax.barh(yy - 0.2, b.values, 0.4, color="0.15", label="NEAT1 anchor", edgecolor="none")
ax.set_yticks(yy); ax.set_yticklabels(nmz, fontsize=6)
ax.axvline(0, color="k", lw=0.6); ax.set_xlabel("log2 enrichment vs no-probe control")
ax.legend(frameon=False, fontsize=6.5, loc="lower right")
ax.set_title("(c) spliceosomal snRNAs, gene by gene", loc="left", fontsize=8)

# (d) decomposition
ax = fig.add_subplot(gs_[1, 1])
tab = pd.DataFrame({
    "GC/length": [22.3, 4.1, 8.4, 12.0, 4.7, 1.4, 1.8, 1.8],
    "JXJ 572 Alu": [2.7, 9.9, 0.9, 11.9, 5.9, 0.4, 0.5, 0.3],
}, index=["20old", "40old", "40odd", "40even", "80", "3'tot", "MED", "TOTAL"])
x = np.arange(len(tab))
ax.bar(x - 0.2, tab["GC/length"], 0.4, color="0.6", label="GC and length (technical)", edgecolor="none")
ax.bar(x + 0.2, tab["JXJ 572 Alu"], 0.4, color="#2b6ca3", label="JXJ 572 Alu site", edgecolor="none")
ax.set_xticks(x); ax.set_xticklabels(tab.index, rotation=45, ha="right", fontsize=6.5)
ax.set_ylabel("out-of-sample variance explained (%)")
ax.axvline(4.5, color="k", lw=0.5, ls=":")
ax.set_ylim(0, 30)
ax.text(2.0, 27.5, "NGS4, has JXJ 572", fontsize=6, ha="center")
ax.text(6.0, 27.5, "NGS6, no 572", fontsize=6, ha="center")
ax.legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(0.02, 0.88))
ax.set_title("(d) what predicts a gene's enrichment", loc="left", fontsize=8)
fig.savefig("figures/fig8_anchors.png", dpi=220)
print("wrote figures/fig8_anchors.png")
