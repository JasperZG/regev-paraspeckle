"""Stage 2 of gene_features.py: GC / length per interval from GRCh38.p13 via memmap."""
import sys, time
import numpy as np
import pandas as pd

FA = "ref/GRCh38.p13.genome.fa"
OUT = "work/gene_features_gencode_v39.tsv"
t0 = time.time()

genes = pd.read_pickle("work/_genes.pkl")
iv = pd.read_pickle("work/_iv.pkl")
aux = pd.read_pickle("work/_aux.pkl")

fai = pd.read_csv(FA + ".fai", sep="\t", header=None,
                  names=["name", "length", "offset", "linebases", "linewidth"]).set_index("name")
mm = np.memmap(FA, dtype=np.uint8, mode="r")

LUT_GC = np.zeros(256, np.uint8); LUT_N = np.zeros(256, np.uint8)
for b in b"GCgc": LUT_GC[b] = 1
for b in b"ACGTacgt": LUT_N[b] = 1

def load(name):
    L, off, lb, lw = fai.loc[name, ["length", "offset", "linebases", "linewidth"]]
    nlines = -(-L // lb)
    span = (nlines - 1) * lw + (L - (nlines - 1) * lb)
    chunk = np.asarray(mm[off: off + span])
    keep = chunk != ord("\n")
    keep &= chunk != ord("\r")
    seq = chunk[keep]
    assert seq.size == L, (name, seq.size, L)
    return seq

iv["len"] = 0; iv["gc"] = 0; iv["acgt"] = 0
missing = []
for seqname, idx in iv.groupby("seq").groups.items():
    if seqname not in fai.index:
        missing.append(seqname); continue
    raw = load(seqname)
    L = raw.size
    cgc = np.zeros(L + 1, np.int32); np.cumsum(LUT_GC[raw], dtype=np.int32, out=cgc[1:])
    cn = np.zeros(L + 1, np.int32); np.cumsum(LUT_N[raw], dtype=np.int32, out=cn[1:])
    del raw
    s = np.clip(iv.loc[idx, "s"].values, 0, L)
    e = np.clip(iv.loc[idx, "e"].values, 0, L)
    iv.loc[idx, "len"] = e - s
    iv.loc[idx, "gc"] = cgc[e] - cgc[s]
    iv.loc[idx, "acgt"] = cn[e] - cn[s]
    del cgc, cn
print(f"genome pass done ({time.time()-t0:.0f}s); seqs missing from fasta: {missing}", flush=True)

agg = iv.groupby(["gene_id", "region"])[["len", "gc", "acgt"]].sum().unstack("region")

feat = genes.set_index("gene_id").copy()
feat = feat.join(aux)
for r in ["gene", "exon", "tx", "cds", "utr5", "utr3"]:
    feat[f"len_{r}"] = agg[("len", r)].reindex(feat.index).fillna(0).astype(np.int64)
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
