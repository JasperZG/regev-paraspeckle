"""Genome-wide scan for probe target sites.

Each of the 80 capture oligos is complementary to a 20-nt NEAT1 target site. Any
transcript carrying a near-copy of that target site (sense strand) can anneal to the
oligo directly. Scan GRCh38 primary chromosomes for every target site in both
orientations allowing up to MAXMM mismatches (12-mer seed + Hamming check), assign
hits to overlapping GENCODE v39 genes, and record per gene the best sense hit (gene
strand carries the target site -> capturable) and best antisense hit (control).

Output: work/probe_hits.tsv (one row per hit) and work/probe_best_per_gene.pkl.
"""
import numpy as np, pandas as pd, time, bioframe as bf
from pyfaidx import Fasta

FA = "ref/GRCh38.p13.genome.fa"
MAXMM = 3
K = 12
t0 = time.time()
P = pd.read_csv("work/probes80.tsv", sep="\t")
comp = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(comp)[::-1]
HANDLE_RC = rc("AGCATAAGTTCAACGGGACCGTGGCCGATGTTTCG")[:20]  # RNA that would anneal to the 5' part of the handle
pats = []
for _, r in P.iterrows():
    pats.append((r.probe, "+", r.target))      # genomic + strand carries target  -> sense for + genes
    pats.append((r.probe, "-", rc(r.target)))  # genomic + strand carries rc(target) -> sense for - genes
pats.append(("HANDLE", "+", HANDLE_RC)); pats.append(("HANDLE", "-", rc(HANDLE_RC)))
LUT = np.full(256, 4, np.uint8)
for i, b in enumerate(b"ACGT"): LUT[b] = i
def enc(s): return LUT[np.frombuffer(s.upper().encode(), np.uint8)]
pat_arr = [(p, o, enc(s)) for p, o, s in pats]
# seeds: all K-mers of each pattern (positions 0..20-K)
seed_map = {}
pow4 = (4 ** np.arange(K - 1, -1, -1)).astype(np.int64)
for idx, (p, o, a) in enumerate(pat_arr):
    for st in range(0, len(a) - K + 1):
        code = int((a[st:st + K].astype(np.int64) * pow4).sum())
        seed_map.setdefault(code, []).append((idx, st))
seed_codes = np.array(sorted(seed_map), dtype=np.int64)
print(f"{len(pat_arr)} patterns, {len(seed_codes)} seed codes", flush=True)

fa = Fasta(FA, as_raw=True, sequence_always_upper=True)
chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
hits = []
for ch in chroms:
    seq = str(fa[ch][:])
    a = enc(seq); L = len(a)
    valid = a < 4
    a4 = np.where(valid, a, 0).astype(np.int64)
    # rolling K-mer code via cumulative polynomial: code[i] = sum a[i+j]*4^(K-1-j)
    code = np.zeros(L - K + 1, np.int64)
    for j in range(K):
        code += a4[j:L - K + 1 + j] * int(pow4[j])
    # positions with any N inside window
    nbad = np.convolve((~valid).astype(np.int32), np.ones(K, np.int32), "valid") > 0
    cand = np.nonzero(np.isin(code, seed_codes) & ~nbad)[0]
    for pos in cand:
        for idx, st in seed_map[int(code[pos])]:
            p, o, pa = pat_arr[idx]
            g0 = pos - st
            if g0 < 0 or g0 + len(pa) > L: continue
            w = a[g0:g0 + len(pa)]
            mm = int((w != pa).sum())
            if mm <= MAXMM:
                hits.append((ch, g0, g0 + len(pa), p, o, mm))
    del seq, a, a4, code, valid, nbad
    print(f"{ch}: {len(hits)} cumulative hits ({time.time()-t0:.0f}s)", flush=True)
H = pd.DataFrame(hits, columns=["chrom", "start", "end", "probe", "orient", "mm"]).drop_duplicates()
H.to_csv("work/probe_hits.tsv", sep="\t", index=False)
print("hits:", len(H), H.groupby("mm").size().to_dict())

# assign to genes (any overlap with the gene span) and call sense/antisense relative to gene strand
G = pd.read_csv("work/genes_v39.tsv", sep="\t")
G = G[G.chr.isin(chroms)].rename(columns={"chr": "chrom"})
ov = bf.overlap(H, G[["chrom", "start", "end", "gene_id", "gene_name", "gene_type", "strand"]], how="inner", suffixes=("", "_g"))
ov["sense"] = ov.orient == ov.strand_g
# exon overlap?
E = pd.read_csv("work/exon_union.saf", sep="\t").rename(columns={"Chr": "chrom", "Start": "start", "End": "end", "GeneID": "gene_id"})
E["start"] -= 1
ex_ov = bf.overlap(H, E[["chrom", "start", "end", "gene_id"]], how="inner", suffixes=("", "_e"))
ex_keys = set(zip(ex_ov.chrom, ex_ov.start, ex_ov.probe, ex_ov.orient, ex_ov.gene_id_e))
ov["in_exon"] = [k in ex_keys for k in zip(ov.chrom, ov.start, ov.probe, ov.orient, ov.gene_id_g)]
ov.to_pickle("work/probe_hits_genes.pkl")
best_s = ov[ov.sense].sort_values("mm").groupby("gene_id_g").first()
best_a = ov[~ov.sense].sort_values("mm").groupby("gene_id_g").first()
out = pd.DataFrame({"sense_mm": best_s.mm, "sense_probe": best_s.probe, "sense_in_exon": best_s.in_exon,
                    "anti_mm": best_a.mm, "anti_probe": best_a.probe})
out.to_pickle("work/probe_best_per_gene.pkl")
print("genes with a sense hit <=3 mm:", len(best_s), "; <=1 mm:", (best_s.mm <= 1).sum(), "; antisense <=1 mm:", (best_a.mm <= 1).sum())
print("per-probe number of sense gene hits (<=2 mm):")
print(ov[ov.sense & (ov.mm <= 2)].groupby("probe").gene_id_g.nunique().sort_values(ascending=False).head(25).to_string())
print(f"done ({time.time()-t0:.0f}s)")
