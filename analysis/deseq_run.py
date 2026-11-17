import pandas as pd, numpy as np, warnings, sys
warnings.filterwarnings("ignore")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
cnt = pd.read_csv("data/counts.tsv", sep="\t", index_col=0)
keep = cnt.sum(axis=1) >= 10
X = cnt.loc[keep].T
cond = pd.Series([c.split("_")[1] for c in cnt.columns], index=cnt.columns).replace({"neg":"0neg"})
inf = DefaultInference(n_cpus=2)
out = {}
# Model A: per-condition
meta = pd.DataFrame({"condition": cond})
dds = DeseqDataSet(counts=X, metadata=meta, design="~condition", refit_cooks=True, inference=inf, quiet=True)
dds.deseq2()
for c in ["20old","40old","40odd","40even","80"]:
    st = DeseqStats(dds, contrast=["condition", c, "0neg"], inference=inf, quiet=True)
    st.summary()
    r = st.results_df.copy()
    try:
        st.lfc_shrink(coeff=f"condition[T.{c}]") if False else None
    except Exception as e:
        pass
    out[c] = r
    print(c, r.shape, flush=True)
# Model B: pooled pulldown vs neg
meta2 = pd.DataFrame({"group": np.where(cond=="0neg","neg","probe")}, index=cnt.columns)
dds2 = DeseqDataSet(counts=X, metadata=meta2, design="~group", refit_cooks=True, inference=inf, quiet=True)
dds2.deseq2()
st2 = DeseqStats(dds2, contrast=["group","probe","neg"], inference=inf, quiet=True)
st2.summary()
r2 = st2.results_df.copy()
print(dds2.varm["LFC"].columns.tolist())
st2.lfc_shrink(coeff=dds2.varm["LFC"].columns[1])
r2s = st2.results_df.copy()
r2["lfc_apeglm"] = r2s["log2FoldChange"]; r2["lfcSE_apeglm"] = r2s["lfcSE"]
r2.to_csv("pooled_probe_vs_neg.tsv", sep="\t")
# per-condition shrunk via reference-level refits (neg is reference so coeffs are condition[T.x]? check)
print(dds.varm["LFC"].columns.tolist())
for c in ["20old","40old","40odd","40even","80"]:
    st = DeseqStats(dds, contrast=["condition", c, "0neg"], inference=inf, quiet=True)
    st.summary()
    coeffs = [k for k in dds.varm["LFC"].columns if k == f"condition[T.{c}]"]
    if coeffs:
        st.lfc_shrink(coeff=coeffs[0])
        out[c]["lfc_apeglm"] = st.results_df["log2FoldChange"]
        out[c]["lfcSE_apeglm"] = st.results_df["lfcSE"]
    out[c].to_csv(f"deseq_{c}_vs_neg.tsv", sep="\t")
# normalized counts
nc = pd.DataFrame(dds.layers["normed_counts"], index=X.index, columns=X.columns).T
nc.to_csv("normed_counts.tsv", sep="\t")
print("size factors", dict(zip(X.index, np.round(dds.obs["size_factors"].values,3))))
