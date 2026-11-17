"""Probe-set membership logic.

80 = 40even U 40odd (disjoint halves). So if a transcript is captured by one
particular capture oligo rather than by paraspeckle proximity:
    up in 40even and 80 but not 40odd  -> an even-numbered oligo
    up in 40odd  and 80 but not 40even -> an odd-numbered oligo
and the gene should carry a stretch of NEAT1-sense sequence for that oligo to
anneal to. Paraspeckle proximity cannot produce a half-set-specific pattern.
"""
import pandas as pd, numpy as np, time
from scipy import stats

M = pd.read_pickle("master.pkl")
res = pd.read_pickle("cqn_response.pkl")
Y = res[["y_20old", "y_40old", "y_40odd", "y_40even", "y_80"]]
Z = (Y - Y.mean()) / Y.std()
nm = M.gene_name.reindex(Z.index)

even_only = Z[(Z.y_40even > 2) & (Z.y_80 > 1) & (Z.y_40odd < 0.75)]
odd_only = Z[(Z.y_40odd > 2) & (Z.y_80 > 1) & (Z.y_40even < 0.75)]
both = Z[(Z.y_40even > 1.5) & (Z.y_40odd > 1.5)]
flat = Z[(Z.abs() < 0.4).all(axis=1)]
print("even-oligo pattern: %d genes" % len(even_only))
print("odd-oligo pattern : %d genes" % len(odd_only))
print("both halves       : %d genes" % len(both))
print("flat (control)    : %d genes" % len(flat))
print("\neven-only:", ", ".join(sorted(str(x) for x in nm[even_only.index])[:40]))
print("\nodd-only :", ", ".join(sorted(str(x) for x in nm[odd_only.index])[:40]))
print("\nboth     :", ", ".join(sorted(str(x) for x in nm[both.index])[:40]))

groups = {"even_only": even_only.index, "odd_only": odd_only.index, "both_halves": both.index,
          "flat_control": flat.index}
pd.to_pickle(groups, "probe_groups.pkl")

# ---- sensitive NEAT1 similarity on these groups ----
FA = "ref/GRCh38.p13.genome.fa"
fai = pd.read_csv(FA + ".fai", sep="\t", header=None,
                  names=["name", "length", "offset", "linebases", "linewidth"]).set_index("name")
mm = np.memmap(FA, dtype=np.uint8, mode="r")
_c = {}
def chrom(name):
    if name not in _c:
        L, off, lb, lw = fai.loc[name, ["length", "offset", "linebases", "linewidth"]]
        nl = -(-L // lb); span = (nl - 1) * lw + (L - (nl - 1) * lb)
        ch = np.asarray(mm[off: off + span]); ch = ch[(ch != 10) & (ch != 13)]
        if len(_c) > 2: _c.clear()
        _c[name] = ch
    return _c[name]
COMP = np.zeros(256, np.uint8)
for a, b in zip(b"ACGTNacgtn", b"TGCANtgcan"): COMP[a] = b
UP = np.arange(256, dtype=np.uint8)
for a, b in zip(b"acgtn", b"ACGTN"): UP[a] = b
def seq_of(ch, blocks, strand):
    s = UP[np.concatenate([ch[a:b] for a, b in blocks])]
    return (COMP[s][::-1] if strand == "-" else s).tobytes().decode()

genes = pd.read_pickle("work/_genes.pkl").set_index("gene_id")
nrow = genes[genes.gene_name == "NEAT1"].iloc[0]
neat = seq_of(chrom("chr11"), [(nrow.start - 1, nrow.end)], "+")
K = 16
nset = {neat[i:i+K] for i in range(len(neat) - K + 1) if "N" not in neat[i:i+K]}
rcn = neat.translate(str.maketrans("ACGT", "TGCA"))[::-1]
aset = {rcn[i:i+K] for i in range(len(rcn) - K + 1) if "N" not in rcn[i:i+K]}

iv = pd.read_pickle("work/_iv.pkl")
want = set().union(*[set(v) for v in groups.values()])
tx = iv[(iv.region == "tx") & iv.gene_id.isin(want)].copy()
tx["chrom"] = tx.gene_id.map(genes.seq)
tx = tx.sort_values(["chrom", "gene_id", "s"])
rows = {}
t0 = time.time()
for gid, g in tx.groupby("gene_id"):
    gi = genes.loc[gid]
    s = seq_of(chrom(gi.seq), list(zip(g.s.values, g.e.values)), gi.strand)
    km = [s[i:i+K] for i in range(len(s) - K + 1)]
    ns = sum(1 for x in km if x in nset); na = sum(1 for x in km if x in aset)
    rows[gid] = (ns, na, len(s))
D = pd.DataFrame(rows, index=["sense16", "anti16", "txlen"]).T
D.to_pickle("offtarget_kmer.pkl")
print("\n== NEAT1 16-mer sharing by probe-set group (%.0fs) ==" % (time.time() - t0))
for lab, idx in groups.items():
    d = D.reindex(idx).dropna()
    print(f"  {lab:14s} n={len(d):4d}  any sense hit: {(d.sense16>0).mean():.2f}  "
          f"median sense={d.sense16.median():.0f} mean={d.sense16.mean():.1f}  "
          f"any anti: {(d.anti16>0).mean():.2f}  median txlen={d.txlen.median():.0f}")
