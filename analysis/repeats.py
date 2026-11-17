import pandas as pd, numpy as np, bioframe as bf, time
t0=time.time()
cols = ["bin","sw","div","dele","ins","chrom","start","end","left","strand","name","cls","fam","rs","re","rl","id"]
rm = pd.read_csv("ext/rmsk.txt.gz", sep="\t", header=None, names=cols, usecols=["chrom","start","end","strand","name","cls","fam","div"])
rm = rm[rm.chrom.str.match(r"^chr[0-9XY]+$")]
print("rmsk", len(rm), time.time()-t0, flush=True)
iv = pd.read_pickle("work/_iv.pkl")
genes = pd.read_pickle("work/_genes.pkl").set_index("gene_id")
M = pd.read_pickle("master.pkl")
use = M.index[M.on_primary_chrom & (M.seq!="chrM") & (M.baseMean.fillna(0) > 1)]
iv = iv[iv.gene_id.isin(use) & iv.region.isin(["exon","utr3","gene","cds","utr5"])].rename(columns={"seq":"chrom","s":"start","e":"end"})
iv["gstrand"] = iv.gene_id.map(genes.strand)
out = pd.DataFrame(index=use)
def cls_of(r):
    if r.fam=="Alu": return "Alu"
    return r.cls
rm["grp"] = np.select([rm.fam=="Alu", rm.cls=="LINE", rm.cls=="LTR", rm.cls=="DNA", rm.cls=="Simple_repeat", rm.cls=="Low_complexity", rm.cls=="SINE"],
                      ["Alu","LINE","LTR","DNA","Simple","LowComplex","otherSINE"], "other")
for reg in ["gene","exon","utr3","cds","utr5"]:
    R = iv[iv.region==reg][["chrom","start","end","gene_id","gstrand"]]
    ov = bf.overlap(R, rm[["chrom","start","end","strand","grp","name"]], how="inner", suffixes=("","_r"))
    ov["bp"] = np.minimum(ov.end, ov.end_r) - np.maximum(ov.start, ov.start_r)
    ov["sense"] = ov.strand_r == ov.gstrand
    bp = ov.groupby(["gene_id","grp_r"]).bp.sum().unstack(fill_value=0)
    L = R.assign(L=R.end-R.start).groupby("gene_id").L.sum()
    for g in ["Alu","LINE","LTR","DNA","Simple","LowComplex"]:
        if g in bp: out[f"{reg}_frac_{g}"] = (bp[g].reindex(use).fillna(0) / L.reindex(use)).fillna(0)
    out[f"{reg}_frac_anyrep"] = (ov.drop_duplicates(["gene_id","start","end","start_r","end_r"]).groupby("gene_id").bp.sum().reindex(use).fillna(0)/L.reindex(use)).fillna(0)
    alu = ov[ov.grp_r=="Alu"].drop_duplicates(["gene_id","start_r","end_r"])
    ns = alu.groupby("gene_id").sense.agg(n_sense="sum", n="size")
    ns["n_anti"] = ns.n - ns.n_sense
    out[f"{reg}_nAlu"] = ns.n.reindex(use).fillna(0)
    out[f"{reg}_nAlu_sense"] = ns.n_sense.reindex(use).fillna(0)
    out[f"{reg}_nAlu_anti"] = ns.n_anti.reindex(use).fillna(0)
    out[f"{reg}_IRAlu_pairs"] = out[f"{reg}_nAlu_sense"]*out[f"{reg}_nAlu_anti"]
    out[f"{reg}_has_IRAlu"] = ((out[f"{reg}_nAlu_sense"]>0)&(out[f"{reg}_nAlu_anti"]>0)).astype(int)
    print(reg, "done", round(time.time()-t0), flush=True)
# close inverted pairs (<=2 kb apart) inside the gene body (pre-mRNA dsRNA potential) and in exons
alu_all = rm[rm.grp=="Alu"][["chrom","start","end","strand"]].sort_values(["chrom","start"]).reset_index(drop=True)
for reg in ["gene","exon"]:
    R = iv[iv.region==reg][["chrom","start","end","gene_id"]]
    ov = bf.overlap(R, alu_all, how="inner", suffixes=("","_r"))
    a = ov[["gene_id","chrom","start_r","end_r","strand_r"]].drop_duplicates(["gene_id","start_r","end_r"]).sort_values(["gene_id","start_r"])
    cnt = {}
    for gid, g in a.groupby("gene_id"):
        s = g.start_r.values; e = g.end_r.values; st = g.strand_r.values
        if len(s) < 2: cnt[gid]=0; continue
        paired = 0
        for i in range(len(s)):
            j = np.searchsorted(s, s[i]-2300); k = np.searchsorted(s, e[i]+2000)
            if np.any(st[j:k] != st[i]): paired += 1
        cnt[gid] = paired
    out[f"{reg}_nAlu_inverted_within2kb"] = pd.Series(cnt).reindex(use).fillna(0)
    print(reg, "pairs done", round(time.time()-t0), flush=True)
out.to_pickle("repeat_features.pkl")
print(out.describe().T[["mean","50%","max"]].round(3).to_string())
