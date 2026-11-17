"""Download ENCODE eCLIP IDR peaks for the paraspeckle proteins and controls,
and count peaks per gene."""
import json, subprocess, os, pandas as pd, numpy as np, bioframe as bf

EXPS = {
    "NONO_K562": "ENCSR861PAR", "SFPQ_HepG2": "ENCSR965DLL", "FUS_K562": "ENCSR069EVH",
    "FUS_HepG2": "ENCSR464OSH", "TARDBP_K562": "ENCSR720BJU", "TARDBP_HepG2": "ENCSR187VEQ",
    "HNRNPK_K562": "ENCSR953ZOA", "HNRNPK_HepG2": "ENCSR828ZID", "MATR3_K562": "ENCSR440SUX",
    "MATR3_HepG2": "ENCSR290VLT", "HNRNPU_K562": "ENCSR520BZQ", "HNRNPU_HepG2": "ENCSR240MVJ",
    "EWSR1_K562": "ENCSR887LPK", "TAF15_K562": "ENCSR568DZW", "AKAP8L_K562": "ENCSR206RXT",
    "ZC3H11A_K562": "ENCSR712IAG",
    # controls: abundant splicing / general RBPs, not paraspeckle proteins
    "PTBP1_K562": "ENCSR981WKN", "U2AF2_K562": "ENCSR893RAV", "SRSF1_K562": "ENCSR321PWZ",
    "SUGP2_HepG2": "ENCSR506UPY",
}
os.makedirs("/home/claude/ps/ext/eclip", exist_ok=True)
got = {}
for lab, acc in EXPS.items():
    out = f"/home/claude/ps/ext/eclip/{lab}.bed.gz"
    if os.path.exists(out) and os.path.getsize(out) > 1000:
        got[lab] = out; continue
    q = (f"https://www.encodeproject.org/search/?type=File&dataset=/experiments/{acc}/"
         "&file_type=bed+narrowPeak&assembly=GRCh38"
         "&status=released&format=json&limit=all&field=accession&field=href&field=assembly")
    try:
        j = json.loads(subprocess.check_output(["curl", "-sS", "-m", "90", "-H", "Accept: application/json", q]))
        fs = [f for f in j["@graph"] if f.get("assembly") == "GRCh38"]
        if not fs:
            print("  no IDR peaks for", lab); continue
        href = sorted(fs, key=lambda f: f["accession"])[0]["href"]
        subprocess.check_call(["curl", "-sS", "-L", "-m", "300", "-o", out, "https://www.encodeproject.org" + href])
        got[lab] = out
        print("  got", lab, os.path.getsize(out))
    except Exception as e:
        print("  FAIL", lab, e)

M = pd.read_pickle("/home/claude/ps/master.pkl")
res = pd.read_pickle("/home/claude/ps/cqn_response.pkl")
iv = pd.read_pickle("/home/claude/work/_iv.pkl")
use = set(res.index)
G = iv[(iv.region == "gene") & iv.gene_id.isin(use)].rename(columns={"seq": "chrom", "s": "start", "e": "end"})[
    ["chrom", "start", "end", "gene_id"]]
E = iv[(iv.region == "exon") & iv.gene_id.isin(use)].rename(columns={"seq": "chrom", "s": "start", "e": "end"})[
    ["chrom", "start", "end", "gene_id"]]
out = pd.DataFrame(index=sorted(use))
for lab, path in got.items():
    try:
        pk = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 6],
                         names=["chrom", "start", "end", "sig"])
    except Exception as e:
        print("  parse fail", lab, e); continue
    pk = pk[pk.chrom.str.match(r"^chr[0-9XY]+$")]
    for reg, df in [("gene", G), ("exon", E)]:
        ov = bf.overlap(df, pk, how="inner", suffixes=("", "_p"))
        n = ov.groupby("gene_id").size()
        out[f"{lab}_{reg}_npk"] = n.reindex(out.index).fillna(0)
    print(f"  {lab}: {len(pk)} peaks -> {(out[f'{lab}_gene_npk']>0).sum()} genes with >=1")
out.to_pickle("/home/claude/ps/eclip_features.pkl")
print(out.shape)
