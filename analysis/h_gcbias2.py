"""Permutation test on the sample labels.

With 2 neg and 10 probe libraries there are C(12,2)=66 ways to split the samples
2-vs-10. Compute rho(GC, log-ratio) for every split and ask where the true split
(neg vs probe) sits. If the true split is not extreme, the GC association is not
separable from library-level technical variation.
"""
import pandas as pd, numpy as np, itertools
from scipy import stats

M = pd.read_pickle("master.pkl")
cnt = pd.read_csv("data/counts.tsv", sep="\t", index_col=0)
K = M[M.keep]
c = cnt.loc[K.index]
cpm = c / c.sum() * 1e6
ok = c.min(axis=1) >= 5
cp = cpm[ok]
gc = K.gc_gene[ok].values
ln = K.log10_len_gene[ok].values
samples = list(cnt.columns)
negs = ["NEAT1_neg_Rep1", "NEAT1_neg_Rep2"]
print("genes:", ok.sum())

rows = []
for pair in itertools.combinations(samples, 2):
    rest = [s for s in samples if s not in pair]
    x = np.log2((cp[rest].mean(axis=1) + 1) / (cp[list(pair)].mean(axis=1) + 1)).values
    rows.append(dict(pair=" + ".join(p.replace("NEAT1_", "") for p in pair),
                     is_true=set(pair) == set(negs),
                     rho_gc=stats.spearmanr(gc, x)[0], rho_len=stats.spearmanr(ln, x)[0],
                     same_cond=pair[0].rsplit("_", 1)[0] == pair[1].rsplit("_", 1)[0]))
P = pd.DataFrame(rows)
true = P[P.is_true].iloc[0]
print(f"\nTRUE split (neg vs rest): rho_GC = {true.rho_gc:+.3f}   rho_len = {true.rho_len:+.3f}")
for col in ["rho_gc", "rho_len"]:
    v = P[col].values
    t = true[col]
    print(f"  {col}: rank {int((np.abs(v) >= abs(t)).sum())}/66 by |rho|   "
          f"one-sided p = {( (v <= t).mean() if t < 0 else (v >= t).mean() ):.3f}   "
          f"null spread: SD={v.std():.3f} range [{v.min():+.3f},{v.max():+.3f}]")
print("\n  the 8 most GC-negative splits:")
print(P.reindex(P.rho_gc.sort_values().index).head(8)[["pair", "is_true", "same_cond", "rho_gc", "rho_len"]].round(3).to_string(index=False))

# restrict the null to splits that pair the two replicates of one condition -> like-for-like
rep = P[P.same_cond]
print("\n  replicate-pair splits only (n=%d):" % len(rep))
print(rep[["pair", "is_true", "rho_gc", "rho_len"]].round(3).to_string(index=False))

# the same test for NEAT1 itself (positive control)
neat_id = M.index[M.gene_name == "NEAT1"][0]
nn = cpm.loc[neat_id]
rows = []
for pair in itertools.combinations(samples, 2):
    rest = [s for s in samples if s not in pair]
    rows.append((set(pair) == set(negs), np.log2(nn[rest].mean() / nn[list(pair)].mean())))
Q = pd.DataFrame(rows, columns=["is_true", "lfc"])
print("\nNEAT1 positive control: true split log2 = %+.2f ; rank %d/66 (most extreme = 1)"
      % (Q[Q.is_true].lfc.iloc[0], int((Q.lfc >= Q[Q.is_true].lfc.iloc[0]).sum())))
