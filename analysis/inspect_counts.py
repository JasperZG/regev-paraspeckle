import pandas as pd, numpy as np
pd.set_option("display.width", 250)
M = pd.read_pickle("master.pkl")
cnt = pd.read_csv("data/counts.tsv", sep="\t", index_col=0)
res = pd.read_pickle("cqn_response.pkl")
ps = pd.read_pickle("probeset_structure.pkl")
short = {c: c.replace("NEAT1_", "").replace("_Rep", "r") for c in cnt.columns}
C = cnt.rename(columns=short)
cpm = (C / C.sum() * 1e6).round(1)

setspec = ["PSRC1", "H3-3A", "P3H4", "RHBDD3", "SIRT7", "SIRT3", "MTRFR", "TUB", "OLMALINC", "RABL2B"]
consist = ["EMSY-DT", "CCDC88B", "RABL2A", "H2AZ1-DT", "HIC2", "MARCHF9", "IL15RA", "CMTM6", "TMEM9B-AS1", "PXN-AS1"]
for lab, names in [("SET-SPECIFIC (large spread)", setspec), ("CONSISTENT across all sets", consist)]:
    ids = [M.index[M.gene_name == n][0] for n in names if (M.gene_name == n).any()]
    print("\n==", lab, "raw counts ==")
    d = C.loc[ids]; d.index = [M.gene_name[i] for i in ids]
    print(d.to_string())
    print("-- CPM --")
    d2 = cpm.loc[ids]; d2.index = [M.gene_name[i] for i in ids]
    print(d2.to_string())

print("\n== library totals (M reads) ==")
print((C.sum() / 1e6).round(2).to_string())
print("\n== how many genes are driven by one sample? ==")
K = M.loc[res.index]
sub = C.loc[res.index]
frac_max = sub.max(axis=1) / sub.sum(axis=1)
print("  genes where one library holds >40% of all counts:", (frac_max > 0.4).sum())
top200 = res.y_cqn.sort_values(ascending=False).head(200).index
print("  among the top 200 by response:", (frac_max[top200] > 0.4).sum())
print("  median frac_max overall %.3f ; top200 %.3f" % (frac_max.median(), frac_max[top200].median()))
