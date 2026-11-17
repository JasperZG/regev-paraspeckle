"""Short-motif scan on the batch-robust response.

batch 1 = {20old, 40odd} vs neg   (neg sits in this batch: within-batch contrast)
batch 2 = {40old, 40even, 80} vs neg
A motif has to predict BOTH to count.

Tests the published candidate elements too: C-rich / HNRNPK (Jacq et al. 15-mer has
a 5' C stretch), UG repeats (TDP-43), AG-rich (paraspeckle shell, West et al.),
polypyrimidine (SFPQ).
"""
import pandas as pd, numpy as np, itertools, time
from scipy import stats

M = pd.read_pickle("master.pkl")
res = pd.read_pickle("cqn_response.pkl")
y1 = res[["y_20old", "y_40odd"]].mean(axis=1)
y2 = res[["y_40old", "y_40even", "y_80"]].mean(axis=1)
print("batch1 vs batch2 response r = %.3f" % np.corrcoef(y1, y2)[0, 1])

FA = "ref/GRCh38.p13.genome.fa"
fai = pd.read_csv(FA + ".fai", sep="\t", header=None,
                  names=["name", "length", "offset", "linebases", "linewidth"]).set_index("name")
mm = np.memmap(FA, dtype=np.uint8, mode="r")
_c = {}
def chrom(n):
    if n not in _c:
        L, off, lb, lw = fai.loc[n, ["length", "offset", "linebases", "linewidth"]]
        nl = -(-L // lb); span = (nl - 1) * lw + (L - (nl - 1) * lb)
        ch = np.asarray(mm[off:off + span]); ch = ch[(ch != 10) & (ch != 13)]
        if len(_c) > 2: _c.clear()
        _c[n] = ch
    return _c[n]
COMP = np.zeros(256, np.uint8)
for a, b in zip(b"ACGTNacgtn", b"TGCANtgcan"): COMP[a] = b
UP = np.arange(256, dtype=np.uint8)
for a, b in zip(b"acgtn", b"ACGTN"): UP[a] = b
def seq_of(ch, blocks, strand):
    s = UP[np.concatenate([ch[a:b] for a, b in blocks])]
    return (COMP[s][::-1] if strand == "-" else s).tobytes().decode()

genes = pd.read_pickle("work/_genes.pkl").set_index("gene_id")
iv = pd.read_pickle("work/_iv.pkl")
use = list(res.index)
seqs = {}
for region in ["utr3", "tx"]:
    sub = iv[(iv.region == region) & iv.gene_id.isin(use)].copy()
    sub["chrom"] = sub.gene_id.map(genes.seq)
    sub = sub.sort_values(["chrom", "gene_id", "s"])
    d = {}
    t0 = time.time()
    for gid, g in sub.groupby("gene_id"):
        gi = genes.loc[gid]
        d[gid] = seq_of(chrom(gi.seq), list(zip(g.s.values, g.e.values)), gi.strand)
    seqs[region] = d
    print(region, len(d), "seqs", round(time.time() - t0), "s", flush=True)
pd.to_pickle(seqs, "seqs.pkl")

K = 6
for region in ["utr3", "tx"]:
    d = seqs[region]
    ids = [g for g in use if g in d and len(d[g]) >= 100]
    yy1 = y1.reindex(ids).values; yy2 = y2.reindex(ids).values
    gc = M.gc_gene.reindex(ids).values
    ln = np.log10([len(d[g]) for g in ids])
    rows = []
    for mer in ("".join(p) for p in itertools.product("ACGT", repeat=K)):
        f = np.array([d[g].count(mer) for g in ids], float)
        dens = f / np.maximum([len(d[g]) for g in ids], 1) * 1000
        if (f > 0).mean() < 0.05: continue
        rows.append((mer, stats.spearmanr(dens, yy1)[0], stats.spearmanr(dens, yy2)[0],
                     stats.spearmanr(dens, gc)[0]))
    T = pd.DataFrame(rows, columns=["mer", "rho1", "rho2", "rho_gc"])
    T["both"] = np.where(np.sign(T.rho1) == np.sign(T.rho2), np.minimum(T.rho1.abs(), T.rho2.abs()), 0) * np.sign(T.rho1)
    T = T.sort_values("both")
    print(f"\n== {region}: 6-mer density vs response (n={len(ids)}) ==")
    print("  most POSITIVE (consistent in both batches):")
    print(T.tail(12).iloc[::-1].round(3).to_string(index=False))
    print("  most NEGATIVE:")
    print(T.head(8).round(3).to_string(index=False))
    T.to_csv(f"res_motifs_{region}.tsv", sep="\t", index=False)
