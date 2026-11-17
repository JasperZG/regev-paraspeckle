"""
Per-gene length and GC-content features from GENCODE v39 (chr_patch_hapl_scaff)
and the matching GRCh38.p13 genome.

GENCODE v39 "ALL" has exactly 67,928 genes, the same Geneid list and order as the
DE_condition_*.csv files, so the output joins 1:1 on Geneid.

Regions (all 0-based half-open, strand-aware where it matters):
  gene        whole gene span (start..end)
  exon        union of exons across every transcript of the gene
              (this is what featureCounts counts reads on)
  intron      gene span minus union exons
  tx          exons of the Ensembl_canonical transcript
  utr5 / cds / utr3   canonical transcript split at CDS (CDS includes stop codon)
  prom        canonical TSS +/- 1 kb
  win100k     100 kb window centred on the gene midpoint (isochore context)

GC = (G+C) / (A+C+G+T); N bases excluded. Soft-masked (lowercase) bases counted.
GC is set to NaN where a region has < 20 called bases.
"""
import gzip, re, sys, time
import numpy as np
import pandas as pd
from pyfaidx import Fasta

GTF = "ref/gencode.v39.chr_patch_hapl_scaff.annotation.gtf.gz"
FA = "ref/GRCh38.p13.genome.fa"
OUT = "work/gene_features_gencode_v39.tsv"

t0 = time.time()
cols = ["seq", "src", "feat", "start", "end", "score", "strand", "frame", "attr"]
keep = {"gene", "transcript", "exon", "CDS", "UTR", "stop_codon"}
rows = []
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] in keep:
            rows.append(f)
gtf = pd.DataFrame(rows, columns=cols)
del rows
gtf["start"] = gtf["start"].astype(np.int64)
gtf["end"] = gtf["end"].astype(np.int64)
gtf["gene_id"] = gtf["attr"].str.extract(r'gene_id "([^"]+)"', expand=False)
gtf["transcript_id"] = gtf["attr"].str.extract(r'transcript_id "([^"]+)"', expand=False)
print(f"parsed GTF: {len(gtf):,} rows  ({time.time()-t0:.0f}s)", flush=True)

# ---------- gene table ----------
genes = gtf[gtf.feat == "gene"].copy()
genes["gene_name"] = genes["attr"].str.extract(r'gene_name "([^"]+)"', expand=False)
genes["gene_type"] = genes["attr"].str.extract(r'gene_type "([^"]+)"', expand=False)
genes["level"] = genes["attr"].str.extract(r"level (\d)", expand=False).astype(int)
genes = genes[["gene_id", "gene_name", "gene_type", "level", "seq", "start", "end", "strand"]]
genes = genes.reset_index(drop=True)
genes["order"] = np.arange(len(genes))
assert genes.gene_id.is_unique and len(genes) == 67928

tx = gtf[gtf.feat == "transcript"]
canon_ids = set(tx.loc[tx.attr.str.contains('tag "Ensembl_canonical"'), "transcript_id"])
n_tx = tx.groupby("gene_id").size().rename("n_transcripts")
canon_tx = tx[tx.transcript_id.isin(canon_ids)][["gene_id", "transcript_id", "start", "end", "strand"]]
assert canon_tx.gene_id.is_unique and len(canon_tx) == len(genes)
canon_tx = canon_tx.set_index("gene_id")

# ---------- interval helpers ----------
def merge(df):
    """merge overlapping/adjacent 0-based half-open intervals within each gene"""
    df = df.sort_values(["gene_id", "s"]).reset_index(drop=True)
    prev_max = df.groupby("gene_id")["e"].cummax().groupby(df["gene_id"]).shift(1)
    new = prev_max.isna() | (df["s"] > prev_max)
    df["blk"] = new.cumsum()
    out = df.groupby("blk").agg(gene_id=("gene_id", "first"), seq=("seq", "first"),
                                s=("s", "min"), e=("e", "max")).reset_index(drop=True)
    return out

ivs = []  # long table: gene_id, seq, s, e, region

g = genes.assign(s=genes.start - 1, e=genes.end)
ivs.append(g[["gene_id", "seq", "s", "e"]].assign(region="gene"))

ex = gtf[gtf.feat == "exon"][["gene_id", "transcript_id", "seq", "start", "end", "strand"]].copy()
ex["s"] = ex.start - 1
ex["e"] = ex.end
ivs.append(merge(ex[["gene_id", "seq", "s", "e"]]).assign(region="exon"))

cex = ex[ex.transcript_id.isin(canon_ids)].copy()
ivs.append(cex[["gene_id", "seq", "s", "e"]].assign(region="tx"))
n_exons_canon = cex.groupby("gene_id").size().rename("n_exons")

# CDS (+ stop codon) of the canonical transcript
cds = gtf[gtf.feat.isin(["CDS", "stop_codon"]) & gtf.transcript_id.isin(canon_ids)][
    ["gene_id", "seq", "start", "end"]].copy()
cds["s"] = cds.start - 1
cds["e"] = cds.end
cds_m = merge(cds[["gene_id", "seq", "s", "e"]])
ivs.append(cds_m.assign(region="cds"))
bounds = cds_m.groupby("gene_id").agg(lo=("s", "min"), hi=("e", "max"))

# UTRs: exonic parts of the canonical transcript outside [lo, hi)
c2 = cex.merge(bounds, left_on="gene_id", right_index=True, how="inner")
left = c2.assign(e=np.minimum(c2.e, c2.lo))
left = left[left.e > left.s]
right = c2.assign(s=np.maximum(c2.s, c2.hi))
right = right[right.e > right.s]
plus_l = left.strand == "+"
plus_r = right.strand == "+"
ivs.append(left[plus_l][["gene_id", "seq", "s", "e"]].assign(region="utr5"))
ivs.append(left[~plus_l][["gene_id", "seq", "s", "e"]].assign(region="utr3"))
ivs.append(right[plus_r][["gene_id", "seq", "s", "e"]].assign(region="utr3"))
ivs.append(right[~plus_r][["gene_id", "seq", "s", "e"]].assign(region="utr5"))

# promoter: canonical TSS +/- 1 kb
ct = canon_tx.reset_index().merge(genes[["gene_id", "seq"]], on="gene_id")
tss0 = np.where(ct.strand == "+", ct.start - 1, ct.end - 1)
ivs.append(pd.DataFrame({"gene_id": ct.gene_id, "seq": ct.seq,
                         "s": tss0 - 1000, "e": tss0 + 1000, "region": "prom"}))
# 100 kb window centred on gene midpoint
mid = ((genes.start - 1 + genes.end) // 2).values
ivs.append(pd.DataFrame({"gene_id": genes.gene_id, "seq": genes.seq,
                         "s": mid - 50_000, "e": mid + 50_000, "region": "win100k"}))

iv = pd.concat(ivs, ignore_index=True)
print(f"built {len(iv):,} intervals ({time.time()-t0:.0f}s)", flush=True)
genes.to_pickle("work/_genes.pkl")
iv.to_pickle("work/_iv.pkl")
pd.concat([n_tx, n_exons_canon], axis=1).to_pickle("work/_aux.pkl")
if len(sys.argv) > 1 and sys.argv[1] == "gtf":
    sys.exit(0)

# ---------- genome pass (run as a fresh process: `python gene_features.py genome`) ----------


agg = iv.groupby(["gene_id", "region"])[["len", "gc", "acgt"]].sum().unstack("region")

feat = genes.set_index("gene_id").copy()
feat = feat.join(n_tx).join(n_exons_canon)
for r in ["gene", "exon", "tx", "cds", "utr5", "utr3"]:
    feat[f"len_{r}"] = agg[("len", r)].reindex(feat.index).fillna(0).astype(np.int64)
# intron = gene span minus union exons
feat["len_intron"] = feat.len_gene - feat.len_exon
gc_int = agg[("gc", "gene")] - agg[("gc", "exon")]
n_int = agg[("acgt", "gene")] - agg[("acgt", "exon")]

def frac(gc, n):
    gc = gc.reindex(feat.index); n = n.reindex(feat.index)
    return (gc / n).where(n >= 20)

for r in ["gene", "exon", "tx", "cds", "utr5", "utr3", "prom", "win100k"]:
    feat[f"gc_{r}"] = frac(agg[("gc", r)], agg[("acgt", r)])
feat["gc_intron"] = frac(gc_int, n_int)

for r in ["gene", "exon", "intron", "tx", "cds", "utr5", "utr3"]:
    feat[f"log10_len_{r}"] = np.log10(feat[f"len_{r}"] + 1)

feat["is_canonical_coding"] = feat.len_cds > 0
feat["on_primary_chrom"] = feat.seq.str.match(r"^chr[0-9XYM]+$")
feat = feat.sort_values("order").drop(columns="order")
feat.index.name = "Geneid"
feat.to_csv(OUT, sep="\t", float_format="%.6g")
print(f"wrote {OUT}: {feat.shape}  ({time.time()-t0:.0f}s)")
