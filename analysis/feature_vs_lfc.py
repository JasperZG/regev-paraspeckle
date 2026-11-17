"""
GC content and log-length vs paraspeckle enrichment (DESeq2 log2FoldChange).

Inputs : DATA/DE_condition_*.csv          (DESeq2 results, 67,928 GENCODE v39 genes)
         gene_features_gencode_v39.tsv    (from gene_features.py)
Outputs: OUTDIR/feature_correlations.tsv, OUTDIR/model_r2.tsv, OUTDIR/fig_*.png

Response = normal-prior shrunken LFC (LFC * tau^2 / (tau^2 + lfcSE^2), tau^2 by
method of moments on baseMean > 10 genes) on genes with baseMean > 10 — the set the
concordance analysis recommended.
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

DATA = sys.argv[1] if len(sys.argv) > 1 else "work/data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "work/out"
FEAT = "work/gene_features_gencode_v39.tsv"
os.makedirs(OUT, exist_ok=True)
BM = 10

# ---------------- load ----------------
feat = pd.read_csv(FEAT, sep="\t", index_col=0)
files = sorted(glob.glob(os.path.join(DATA, "DE_condition_*.csv")))
de = {}
for f in files:
    d = pd.read_csv(f).set_index("Geneid")
    assert d.index.equals(feat.index), f"{f}: Geneid order differs from GENCODE v39"
    name = os.path.basename(f)
    ok = (d.baseMean > BM) & d.log2FoldChange.notna()
    v = d.loc[ok, "log2FoldChange"].var() - (d.loc[ok, "lfcSE"] ** 2).mean()
    tau2 = max(v, 1e-3)
    d["lfc_shrunk"] = d.log2FoldChange * tau2 / (tau2 + d.lfcSE ** 2)
    d["keep"] = ok
    de[name] = d
    print(f"{name}: kept {ok.sum():,} genes (baseMean>{BM}); tau={np.sqrt(tau2):.3f}")

main = "DE_condition_80.csv" if "DE_condition_80.csv" in de else list(de)[0]
reps = [k for k in ["DE_condition_40even.csv", "DE_condition_40odd.csv"] if k in de]

FEATS = {
    "gc_exon": "GC, exons (union)",
    "gc_intron": "GC, introns",
    "gc_gene": "GC, gene body",
    "gc_cds": "GC, CDS",
    "gc_utr5": "GC, 5'UTR",
    "gc_utr3": "GC, 3'UTR",
    "gc_prom": "GC, promoter (TSS ±1 kb)",
    "gc_win100k": "GC, 100 kb neighbourhood",
    "log10_len_gene": "log10 gene length",
    "log10_len_exon": "log10 exonic length",
    "log10_len_intron": "log10 intronic length",
    "log10_len_tx": "log10 canonical transcript length",
    "log10_len_cds": "log10 CDS length",
    "log10_len_utr5": "log10 5'UTR length",
    "log10_len_utr3": "log10 3'UTR length",
    "n_exons": "exons in canonical transcript",
}

def rank_resid(y, X):
    """residual of rank(y) after OLS on rank-transformed columns of X"""
    ry = stats.rankdata(y)
    RX = np.column_stack([np.ones(len(y))] + [stats.rankdata(c) for c in X.T])
    beta, *_ = np.linalg.lstsq(RX, ry, rcond=None)
    return ry - RX @ beta

# ---------------- univariate + partial correlations ----------------
rows = []
for fname, d in de.items():
    base = d[d.keep].join(feat, how="inner")
    base = base[base.on_primary_chrom]
    for subset, sub in [("all", base), ("protein_coding", base[base.gene_type == "protein_coding"]),
                        ("lncRNA", base[base.gene_type == "lncRNA"])]:
        for k in FEATS:
            s = sub[[k, "lfc_shrunk", "baseMean"]].dropna()
            if k.startswith("log10_len_utr") or k == "log10_len_cds":
                s = s[s[k] > 0]
            if len(s) < 50:
                continue
            rho, p = stats.spearmanr(s[k], s.lfc_shrunk)
            r = np.corrcoef(s[k], s.lfc_shrunk)[0, 1]
            lb = np.log10(s.baseMean.values)[:, None]
            prho = np.corrcoef(rank_resid(s[k].values, lb), rank_resid(s.lfc_shrunk.values, lb))[0, 1]
            rho_bm = stats.spearmanr(s[k], s.baseMean)[0]
            rows.append(dict(file=fname, subset=subset, feature=k, n=len(s), pearson_r=r,
                             spearman_rho=rho, p=p, partial_rho_given_expr=prho,
                             rho_feature_vs_baseMean=rho_bm))
cor = pd.DataFrame(rows)
cor.to_csv(os.path.join(OUT, "feature_correlations.tsv"), sep="\t", index=False, float_format="%.4g")

# ---------------- multivariate, leave-one-chromosome-out ----------------
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import HistGradientBoostingRegressor

GC = ["gc_exon", "gc_intron", "gc_gene", "gc_cds", "gc_utr5", "gc_utr3", "gc_prom", "gc_win100k"]
LEN = ["log10_len_gene", "log10_len_exon", "log10_len_intron", "log10_len_tx",
       "log10_len_cds", "log10_len_utr5", "log10_len_utr3", "n_exons", "n_transcripts"]

def design(df, cols, with_type):
    X = df[cols].copy()
    for c in cols:  # missing (e.g. UTR GC for non-coding): impute median + indicator
        if X[c].isna().any():
            X[c + "_na"] = X[c].isna().astype(float)
            X[c] = X[c].fillna(X[c].median())
    if with_type:
        top = ["protein_coding", "lncRNA", "processed_pseudogene", "unprocessed_pseudogene"]
        for t in top:
            X["type_" + t] = (df.gene_type == t).astype(float)
    return X.values

def cv_r2(df, cols, model, with_type=False):
    y = df.lfc_shrunk.values
    X = design(df, cols, with_type) if cols or with_type else np.zeros((len(df), 1))
    pred = np.zeros(len(y))
    for chrom in df.seq.unique():
        te = (df.seq == chrom).values
        m = model()
        m.fit(X[~te], y[~te])
        pred[te] = m.predict(X[te])
    return 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()

mrows = []
d = de[main]
df = d[d.keep].join(feat, how="inner")
df = df[df.on_primary_chrom & (df.seq != "chrM")]
blocks = {"gene type only": ([], True), "GC only": (GC, False), "length only": (LEN, False),
          "GC + length": (GC + LEN, False), "GC + length + gene type": (GC + LEN, True)}
for label, (cols, wt) in blocks.items():
    for mname, model in [("linear", LinearRegression),
                         ("boosted trees", lambda: HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05,
                                                                                  max_leaf_nodes=31, random_state=0))]:
        mrows.append(dict(file=main, block=label, model=mname, n=len(df),
                          cv_r2_leave_chrom_out=cv_r2(df, cols, model, wt)))
# reliability ceiling from the split halves
if len(reps) == 2:
    a, b = de[reps[0]], de[reps[1]]
    both = a.keep & b.keep
    r_half = np.corrcoef(a.loc[both, "lfc_shrunk"], b.loc[both, "lfc_shrunk"])[0, 1]
    rel = 2 * r_half / (1 + r_half)  # Spearman-Brown for the full 80-probe set
    mrows.append(dict(file=f"{reps[0]} vs {reps[1]}", block="split-half r", model="-", n=int(both.sum()),
                      cv_r2_leave_chrom_out=r_half))
    mrows.append(dict(file=main, block="reliability ceiling (Spearman-Brown)", model="-", n=int(both.sum()),
                      cv_r2_leave_chrom_out=rel))
model_tab = pd.DataFrame(mrows)
model_tab.to_csv(os.path.join(OUT, "model_r2.tsv"), sep="\t", index=False, float_format="%.4f")

# ---------------- figures (Oded house style) ----------------
def hist2d_feature(x, y, xlabel, ylabel, title, path, bins=100):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.hist2d(x, y, bins=bins, cmap="Blues", norm=LogNorm())
    r = np.corrcoef(x, y)[0, 1]
    ax.text(0.02, 0.98, f"n={len(x)}\nr={r:.3f}", transform=ax.transAxes, va="top", ha="left")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_title(title, loc="center")
    ax.set_box_aspect(1)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

try:
    sys.path.insert(0, DATA); sys.path.insert(0, os.path.dirname(DATA))
    import style  # noqa: F401  (Oded's helper, if present, just to confirm availability)
except Exception:
    pass

base = de[main][de[main].keep].join(feat, how="inner")
base = base[base.on_primary_chrom]
for k in ["gc_exon", "gc_intron", "gc_utr3", "gc_win100k", "log10_len_gene", "log10_len_exon",
          "log10_len_utr3", "log10_len_intron"]:
    s = base[[k, "lfc_shrunk"]].dropna()
    if k in ("log10_len_utr3",):
        s = s[s[k] > 0]
    hist2d_feature(s[k].values, s.lfc_shrunk.values, FEATS[k], "log2FoldChange (shrunk)",
                   f"{k}\n{main}", os.path.join(OUT, f"fig_{k}.png"))

pd.set_option("display.width", 200)
print(cor[(cor.file == main) & (cor.subset == "all")].sort_values("spearman_rho", key=abs, ascending=False)
      .round(3).to_string(index=False))
print(model_tab.round(3).to_string(index=False))
