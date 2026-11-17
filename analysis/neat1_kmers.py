"""Shared k-mer content with NEAT1.

Two distinct things this can mean:
  sense  : the transcript contains NEAT1-like sequence in the SAME orientation as
           NEAT1 -> the antisense capture probes can hybridise to it directly.
           This is a pulldown artifact, not biology.
  anti   : the transcript is complementary to NEAT1 -> it could base-pair with the
           NEAT1 lncRNA itself. This is the direct RNA-RNA interaction hypothesis.
"""
import pandas as pd, numpy as np, time
from scipy import stats

FA = "ref/GRCh38.p13.genome.fa"
fai = pd.read_csv(FA + ".fai", sep="\t", header=None,
                  names=["name", "length", "offset", "linebases", "linewidth"]).set_index("name")
mm = np.memmap(FA, dtype=np.uint8, mode="r")
CACHE = {}

def chrom(name):
    if name in CACHE:
        return CACHE[name]
    L, off, lb, lw = fai.loc[name, ["length", "offset", "linebases", "linewidth"]]
    nlines = -(-L // lb)
    span = (nlines - 1) * lw + (L - (nlines - 1) * lb)
    ch = np.asarray(mm[off: off + span])
    ch = ch[(ch != 10) & (ch != 13)]
    if len(CACHE) > 3: CACHE.clear()
    CACHE[name] = ch
    return ch

COMP = np.zeros(256, np.uint8)
for a, b in zip(b"ACGTNacgtn", b"TGCANtgcan"):
    COMP[a] = b
UP = np.arange(256, dtype=np.uint8)
for a, b in zip(b"acgtn", b"ACGTN"):
    UP[a] = b

def seq_of(ch, blocks, strand):
    s = np.concatenate([ch[a:b] for a, b in blocks]) if blocks else np.array([], np.uint8)
    s = UP[s]
    if strand == "-":
        s = COMP[s][::-1]
    return s.tobytes().decode()

def kmers(s, k):
    return {s[i:i + k] for i in range(len(s) - k + 1) if "N" not in s[i:i + k]}

K = 25
neat = pd.read_pickle("work/_genes.pkl")
nrow = neat[neat.gene_name == "NEAT1"].iloc[0]
ch11 = chrom("chr11")
neat_seq = seq_of(ch11, [(nrow.start - 1, nrow.end)], "+")
print("NEAT1 length", len(neat_seq))
n_sense = kmers(neat_seq, K)
rc = neat_seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
n_anti = kmers(rc, K)
print("NEAT1 %d-mers: %d" % (K, len(n_sense)))

res = pd.read_pickle("cqn_response.pkl")
iv = pd.read_pickle("work/_iv.pkl")
genes = pd.read_pickle("work/_genes.pkl").set_index("gene_id")
use = list(res.index)
tx = iv[(iv.region == "tx") & iv.gene_id.isin(use)].copy()
tx["chrom"] = tx.gene_id.map(genes.seq)
tx = tx.sort_values(["chrom", "gene_id", "s"])
t0 = time.time()
out = {}
for i, (gid, g) in enumerate(tx.groupby("gene_id")):
    gi = genes.loc[gid]
    ch = chrom(gi.seq)
    s = seq_of(ch, list(zip(g.s.values, g.e.values)), gi.strand)
    km = kmers(s, K)
    out[gid] = (len(km & n_sense), len(km & n_anti), len(km))
    if i % 3000 == 0:
        print(i, round(time.time() - t0), flush=True)
D = pd.DataFrame(out, index=["k_sense", "k_anti", "k_tot"]).T
D["f_sense"] = D.k_sense / D.k_tot.clip(lower=1)
D["f_anti"] = D.k_anti / D.k_tot.clip(lower=1)
D.to_pickle("neat1_kmer.pkl")
print(D.describe().T.round(4).to_string())
print("genes with any shared sense 25-mer:", (D.k_sense > 0).sum(), " anti:", (D.k_anti > 0).sum())
y = res.y_cqn.reindex(D.index)
for c in ["k_sense", "k_anti", "f_sense", "f_anti"]:
    print(f"  rho(y, {c}) = {stats.spearmanr(D[c], y)[0]:+.3f}")
hit = D[D.k_sense > 0].index
print("median y with shared sense k-mer: %+.3f (n=%d)  without: %+.3f"
      % (y[hit].median(), len(hit), y[~y.index.isin(hit)].median()))
