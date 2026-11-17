"""Direct test of probe cross-hybridisation.

For each gene in the probe-pattern / uniform / flat groups (sequences fetched from
Ensembl REST, work/pattern_gene_seqs.pkl), scan every one of the 80 NEAT1 target
sites (the 20-nt NEAT1 sense sequence each probe is complementary to) against the
gene's transcripts (cDNA, sense strand) and pre-mRNA (genomic, sense strand), and
record the best ungapped match (fewest mismatches over the 20 nt, and the longest
exact shared substring). A gene captured by direct annealing to probe P should carry
a near-perfect copy of P's target site, and P should belong to the probe sets in
which the gene is enriched.

Control: the same scan on the reverse complement of the target sites (what the RNA
would need to look like to base-pair with NEAT1 itself rather than with the oligo),
and on the flat genes.
"""
import re, zipfile, pandas as pd, numpy as np

seq = open("work/neat1_hg38.fa").read().split("\n", 1)[1].replace("\n", "")
comp = str.maketrans("ACGTacgt", "TGCAtgca")
rc = lambda s: s.translate(comp)[::-1]

def parse(f):
    z = zipfile.ZipFile(f); ss = z.read("xl/sharedStrings.xml").decode()
    strs = [re.sub(r"<[^>]+>", "", s) for s in re.findall(r"<si>(.*?)</si>", ss, flags=re.S)]
    d = z.read("xl/worksheets/sheet1.xml").decode()
    out = []
    for r in re.findall(r"<row[^>]*>(.*?)</row>", d, flags=re.S):
        row = {}
        for col, rn, attrs, inner in re.findall(r'<c r="([A-Z]+)(\d+)"([^>]*?)(?:/>|>(.*?)</c>)', r, flags=re.S):
            v = re.search(r"<v>(.*?)</v>", inner or ""); val = v.group(1) if v else ""
            if 't="s"' in attrs and val: val = strs[int(val)]
            row[col] = val
        out.append(row)
    return out

rows = parse("data/ngs4/NEAT1 probe positions_2026_08_batch.xlsx")
region = None; probes = []
for r in rows:
    if r.get("A"): region = r["A"]
    if r.get("B", "").startswith("JXJ"):
        p = r["B"]; s = r["C"][:20].upper(); tgt = rc(s)
        pos = seq.find(tgt)
        num = int(p.split()[1])
        reg = region
        if p == "JXJ 572": reg = "Region8"
        if p == "JXJ 810": reg = "Region7"
        rn = int(reg.replace("Region", ""))
        series = "4xx" if num < 500 else ("5xx" if num < 700 else "new")
        parity = "odd" if rn % 2 else "even"
        sets = {"20old": series == "4xx", "40old": series in ("4xx", "5xx"),
                "40odd": parity == "odd", "40even": parity == "even", "80": True}
        probes.append(dict(probe=p, region=reg, parity=parity, series=series, pat=f"JXJ{series}_{parity}" if series != "new" else f"new_{parity}",
                           target=tgt, pos=pos, **sets))
P = pd.DataFrame(probes)
P.to_csv("work/probes80.tsv", sep="\t", index=False)
assert (P.pos >= 0).all()
HANDLE = "AGCATAAGTTCAACGGGACCGTGGCCGATGTTTCG"  # common 3' tail of every probe (DNA, as written)

LUT = np.full(256, 4, np.uint8)
for i, b in enumerate(b"ACGT"): LUT[b] = i
def enc(s): return LUT[np.frombuffer(s.upper().encode(), np.uint8)]

def best_match(gene_arr, tgt):
    """min mismatches of tgt (20-mer) over all ungapped positions in gene_arr."""
    k = len(tgt); t = enc(tgt)
    if len(gene_arr) < k: return 99
    W = np.lib.stride_tricks.sliding_window_view(gene_arr, k)
    mm = (W != t).sum(axis=1)
    return int(mm.min())

def longest_common(gene_seq, tgt):
    """longest exact substring of tgt present in gene_seq."""
    L = len(tgt)
    for k in range(L, 9, -1):
        for i in range(L - k + 1):
            if tgt[i:i + k] in gene_seq: return k
    return 0

S = pd.read_pickle("work/pattern_gene_seqs.pkl")
out = []
for gid, rec in S.items():
    seqs = {"cdna": [s for s in rec.get("cdna", []) if s], "genomic": [s for s in rec.get("genomic", []) if s and len(s) < 1_500_000]}
    if not seqs["cdna"] and not seqs["genomic"]: continue
    arrs = {k: [enc(s) for s in v] for k, v in seqs.items()}
    ups = {k: [s.upper() for s in v] for k, v in seqs.items()}
    for _, pr in P.iterrows():
        for orient, tgt in [("sense", pr.target), ("anti", rc(pr.target))]:
            row = dict(gene=gid, name=rec["name"], group=rec["group"], gene_pat=rec["pat"], probe=pr.probe, probe_pat=pr.pat, orient=orient)
            for k in ["cdna", "genomic"]:
                if arrs[k]:
                    row[f"mm_{k}"] = min(best_match(a, tgt) for a in arrs[k])
                    row[f"lcs_{k}"] = max(longest_common(u, tgt) for u in ups[k])
                else:
                    row[f"mm_{k}"] = np.nan; row[f"lcs_{k}"] = np.nan
            out.append(row)
    # handle
    for orient, tgt in [("handle_rc", rc(HANDLE)), ("handle", HANDLE)]:
        row = dict(gene=gid, name=rec["name"], group=rec["group"], gene_pat=rec["pat"], probe="HANDLE", probe_pat="all", orient=orient)
        for k in ["cdna", "genomic"]:
            if arrs[k]:
                row[f"mm_{k}"] = min(best_match(a, tgt) for a in arrs[k]); row[f"lcs_{k}"] = max(longest_common(u, tgt) for u in ups[k])
            else:
                row[f"mm_{k}"] = np.nan; row[f"lcs_{k}"] = np.nan
        out.append(row)
R = pd.DataFrame(out)
R.to_pickle("work/probe_scan.pkl")

pd.set_option("display.width", 250)
print("genes scanned:", R.gene.nunique(), R.groupby("group").gene.nunique().to_dict())
sense = R[(R.orient == "sense") & (R.probe != "HANDLE")]
# best probe per gene (fewest mismatches on cDNA, then genomic)
best = sense.sort_values(["mm_cdna", "mm_genomic"]).groupby("gene").first()
print("\nbest probe match per gene (sense, cDNA mismatches over 20 nt): distribution by group")
print(pd.crosstab(best.group, best.mm_cdna.clip(upper=8)))
print("\nmatched probe set pattern == gene's enrichment pattern?")
best["pat_match"] = best.probe_pat == best.gene_pat
print(pd.crosstab(best.group, best.pat_match))
print("\n--- pattern genes, best probe ---")
print(best[best.group == "pattern"].sort_values("mm_cdna")[["name", "gene_pat", "probe", "probe_pat", "mm_cdna", "lcs_cdna", "mm_genomic", "lcs_genomic"]].to_string())
anti = R[(R.orient == "anti") & (R.probe != "HANDLE")].sort_values(["mm_cdna", "mm_genomic"]).groupby("gene").first()
print("\nantisense control, cDNA mismatches by group:")
print(pd.crosstab(anti.group, anti.mm_cdna.clip(upper=8)))
h = R[R.probe == "HANDLE"]
print("\nhandle: min mismatches over 35 nt by group/orient")
print(h.groupby(["group", "orient"])[["mm_cdna", "mm_genomic"]].min())
