"""Screen the October 2026 probe set (36 oligos, JXJ 1056-1091) the way JXJ 572 was caught:
map each 20-nt target onto NEAT1, flag RepeatMasker overlap, and count genome-wide
sense-strand sites within 2 mismatches (12-mer seed + Hamming), assigned to genes."""
import re, zipfile, numpy as np, pandas as pd, time, bioframe as bf
from pyfaidx import Fasta
z = zipfile.ZipFile("data/ngs4/NEAT1 probe positions_2026_10_batch.xlsx"); ss = z.read("xl/sharedStrings.xml").decode()
strs = [re.sub(r"<[^>]+>", "", s) for s in re.findall(r"<si>(.*?)</si>", ss, flags=re.S)]
d = z.read("xl/worksheets/sheet1.xml").decode(); probes = []
for r in re.findall(r"<row[^>]*>(.*?)</row>", d, flags=re.S):
    row = {}
    for col, rn, attrs, inner in re.findall(r'<c r="([A-Z]+)(\d+)"([^>]*?)(?:/>|>(.*?)</c>)', r, flags=re.S):
        v = re.search(r"<v>(.*?)</v>", inner or ""); val = v.group(1) if v else ""
        if 't="s"' in attrs and val: val = strs[int(val)]
        row[col] = val
    if row.get("A", "").startswith("JXJ"): probes.append((row["A"], row["B"].upper()))
comp = str.maketrans("ACGT", "TGCA"); rc = lambda s: s.translate(comp)[::-1]
HANDLE = "AGCATAAGTTCAACGGGACCGTGGCCGATGTTTCG"
neat1 = open("work/neat1_hg38.fa").read().split("\n", 1)[1].replace("\n", "").upper()
rm = pd.read_csv("work/rmsk_neat1.tsv", sep="\t", header=None, names=["chr", "s", "e", "strand", "name", "cls", "fam", "sw"]); rm["ts"] = rm.s - 65422773; rm["te"] = rm.e - 65422773
rows = []
for p, seq in probes:
    assert seq.endswith(HANDLE), p
    tgt = rc(seq[:-len(HANDLE)]); assert len(tgt) == 20
    pos = neat1.find(tgt); hit = rm[(rm.ts < pos + 20) & (rm.te > pos)] if pos >= 0 else rm.iloc[0:0]
    rows.append(dict(probe=p, target=tgt, pos=pos, gc=(tgt.count("G") + tgt.count("C")) / 20, repeat=";".join(f"{n}({c})" for n, c in zip(hit.name, hit.cls))))
P = pd.DataFrame(rows); print("probes:", len(P), "; all map once to NEAT1 sense:", (P.pos >= 0).all(), "; positions", P.pos.min(), "-", P.pos.max())
alu = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCCGGCTAAAACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAA"
P["alu_mm"] = [min(sum(a != b for a, b in zip(alu[i:i + 20], t)) for i in range(len(alu) - 19)) for t in P.target]
# genome scan
t0 = time.time(); MAXMM = 2; K = 12
pats = [(p, "+", t) for p, t in zip(P.probe, P.target)] + [(p, "-", rc(t)) for p, t in zip(P.probe, P.target)]
LUT = np.full(256, 4, np.uint8)
for i, b in enumerate(b"ACGT"): LUT[b] = i
enc = lambda s: LUT[np.frombuffer(s.encode(), np.uint8)]
pat_arr = [(p, o, enc(s)) for p, o, s in pats]; pow4 = (4 ** np.arange(K - 1, -1, -1)).astype(np.int64); seed_map = {}
for idx, (p, o, a) in enumerate(pat_arr):
    for st in range(0, 20 - K + 1): seed_map.setdefault(int((a[st:st + K].astype(np.int64) * pow4).sum()), []).append((idx, st))
seed_codes = np.array(sorted(seed_map), dtype=np.int64)
fa = Fasta("ref/GRCh38.p13.genome.fa", as_raw=True, sequence_always_upper=True); hits = []
for ch in [f"chr{i}" for i in range(1, 23)] + ["chrX"]:
    a = enc(str(fa[ch][:])); L = len(a); valid = a < 4; a4 = np.where(valid, a, 0).astype(np.int64)
    code = np.zeros(L - K + 1, np.int64)
    for j in range(K): code += a4[j:L - K + 1 + j] * int(pow4[j])
    nbad = np.convolve((~valid).astype(np.int32), np.ones(K, np.int32), "valid") > 0
    for pos in np.nonzero(np.isin(code, seed_codes) & ~nbad)[0]:
        for idx, st in seed_map[int(code[pos])]:
            p, o, pa = pat_arr[idx]; g0 = pos - st
            if g0 < 0 or g0 + 20 > L: continue
            mm = int((a[g0:g0 + 20] != pa).sum())
            if mm <= MAXMM: hits.append((ch, g0, g0 + 20, p, o, mm))
H = pd.DataFrame(hits, columns=["chrom", "start", "end", "probe", "orient", "mm"]).drop_duplicates(); print("scan done", f"{time.time()-t0:.0f}s", len(H), "hits")
G = pd.read_csv("work/genes_v39.tsv", sep="\t"); G = G[G.chr.str.match(r"^chr[0-9X]+$")].rename(columns={"chr": "chrom"})
ov = bf.overlap(H, G[["chrom", "start", "end", "gene_id", "gene_name", "strand"]], how="inner", suffixes=("", "_g")); ov["sense"] = ov.orient == ov.strand_g
E = pd.read_csv("work/exon_union.saf", sep="\t").rename(columns={"Chr": "chrom", "Start": "start", "End": "end", "GeneID": "gene_id"}); E["start"] -= 1
exov = bf.overlap(H, E[["chrom", "start", "end", "gene_id"]], how="inner", suffixes=("", "_e")); exov = exov.merge(G[["gene_id", "strand"]], left_on="gene_id_e", right_on="gene_id"); exov["sense"] = exov.orient == exov.strand
P = P.set_index("probe")
P["genome_sites_le2mm"] = H.groupby("probe").size().reindex(P.index).fillna(0).astype(int)
P["sense_genes_le2mm"] = ov[ov.sense].groupby("probe").gene_id_g.nunique().reindex(P.index).fillna(0).astype(int)
P["sense_exonic_genes_le2mm"] = exov[exov.sense].groupby("probe").gene_id_e.nunique().reindex(P.index).fillna(0).astype(int)
P["sense_exonic_genes_le1mm"] = exov[exov.sense & (exov.mm <= 1)].groupby("probe").gene_id_e.nunique().reindex(P.index).fillna(0).astype(int)
pd.set_option("display.width", 200); print(P.sort_values("sense_exonic_genes_le2mm", ascending=False).to_string())
P.to_csv("work/probes_oct36_screen.tsv", sep="\t")
old = pd.read_pickle("work/probe_hits_genes.pkl"); oldc = old[old.sense & (old.mm <= 2)].groupby("probe").gene_id_g.nunique()
print("\nfor comparison, the old 80 (sense genes <=2 mm): median %d, 90th pct %d, JXJ 572 %d" % (oldc.median(), oldc.quantile(.9), oldc.get("JXJ 572", 0)))
