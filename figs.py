import pandas as pd, numpy as np, itertools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import stats
from style import hist2d_compare

M = pd.read_pickle("master.pkl")
res = pd.read_pickle("cqn_response.pkl")
K = M.loc[res.index].copy()
K["y"] = res.y_cqn
K["y1"] = res[["y_20old", "y_40odd"]].mean(axis=1)
K["y2"] = res[["y_40old", "y_40even", "y_80"]].mean(axis=1)
K["mid"] = (K.start + K.end) / 2
n1 = M[M.gene_name == "NEAT1"].iloc[0]; NPOS = (n1.start + n1.end) / 2

cnt = pd.read_csv("data/counts.tsv", sep="\t", index_col=0)
c = cnt.loc[K.index]; ok = c.min(axis=1) >= 5
cp = (c / c.sum() * 1e6)[ok]; Ks = K[ok]
samples = list(cnt.columns); negs = ["NEAT1_neg_Rep1", "NEAT1_neg_Rep2"]

# ---------------- Fig 1: GC looks like the answer, but it is a library axis ----------
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.6))
hist2d_compare(K.gc_gene.values, K.lfc.values, "gene GC fraction", "log2FoldChange (apeglm)",
               xlabel="gene GC fraction", ylabel="log2FoldChange, probe vs neg", ax=axes[0])
axes[0].set_title("DE_condition_80.csv etc. vs\nGENCODE v39 GC content")

rows = []
for pair in itertools.combinations(samples, 2):
    rest = [s for s in samples if s not in pair]
    x = np.log2((cp[rest].mean(axis=1) + 1) / (cp[list(pair)].mean(axis=1) + 1))
    rows.append((set(pair) == set(negs), stats.spearmanr(Ks.gc_gene, x)[0],
                 pair[0].rsplit("_", 1)[0] == pair[1].rsplit("_", 1)[0]))
P = pd.DataFrame(rows, columns=["is_true", "rho", "same"])
t = P[P.is_true].rho.iloc[0]
ax = axes[1]
ax.hist(P.rho, bins=22, color="0.75", edgecolor="k", linewidth=.4)
ax.axvline(t, color="k", linewidth=1.6)
ax.text(0.02, 0.98, f"all 66 two-vs-ten splits\ntrue split rho={t:.3f}\nrank 16/66, p=0.106",
        transform=ax.transAxes, va="top", ha="left")
ax.set_xlabel("rho(gene GC, log ratio)"); ax.set_ylabel("number of splits")
ax.set_title("sample-label permutation\n(which 2 libraries are the control)")
ax.set_box_aspect(1)

ax = axes[2]
h = ax.hist2d(K.y1.values, K.y2.values, bins=60, cmap="Blues", norm=LogNorm())
r = np.corrcoef(K.y1, K.y2)[0, 1]
ax.text(0.02, 0.98, f"n={len(K)}\nr={r:.3f}", transform=ax.transAxes, va="top", ha="left")
ax.set_xlabel("response, batch 1 (20old+40odd)"); ax.set_ylabel("response, batch 2 (40old+40even+80)")
ax.set_title("CQN-corrected response\nreproducibility across batches")
ax.set_box_aspect(1)
fig.colorbar(h[3], ax=ax)
fig.tight_layout(); fig.savefig("fig1_gc_is_technical.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ---------------- Fig 2: distance decay from the NEAT1 locus ----------------
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.6))
rng = np.random.default_rng(1)
rad = [0.5, 1, 2, 3, 5, 10, 20]
mu, lo, hi, nn = [], [], [], []
for r_ in rad:
    v = K.y[(K.seq == "chr11") & ((K.mid - NPOS).abs() < r_ * 1e6)].values
    bs = [np.mean(rng.choice(v, len(v))) for _ in range(2000)]
    mu.append(v.mean()); lo.append(np.percentile(bs, 2.5)); hi.append(np.percentile(bs, 97.5)); nn.append(len(v))
ax = axes[0]
ax.errorbar(rad, mu, yerr=[np.array(mu) - lo, np.array(hi) - np.array(mu)], fmt="o-", color="k",
            capsize=3, markersize=4, linewidth=1)
ax.axhline(K.y.mean(), color="0.5", linestyle="--", linewidth=1)
ax.set_xscale("log")
ax.text(0.02, 0.98, "genes within r Mb of NEAT1\nerror bars 95% bootstrap\ndashed = genome mean",
        transform=ax.transAxes, va="top", ha="left")
ax.set_xlabel("radius around NEAT1 locus (Mb)"); ax.set_ylabel("mean CQN response")
ax.set_title("chr11:65.43 Mb (NEAT1)\ndistance decay"); ax.set_box_aspect(1)

stat = []
for ch, g in K.groupby("seq"):
    p = g.mid.values; v = g.y.values
    for cc in np.arange(p.min(), p.max(), 1e6):
        m = np.abs(p - cc) < 2e6
        if 60 <= m.sum() <= 160:
            rest = np.concatenate([K.y.values[K.seq.values != ch], v[~m]])
            stat.append((ch, cc, v[m].mean() - rest.mean()))
S = pd.DataFrame(stat, columns=["chrom", "centre", "eff"])
obs = 0.0999
ax = axes[1]
ax.hist(S.eff, bins=18, color="0.75", edgecolor="k", linewidth=.4)
ax.axvline(obs, color="k", linewidth=1.6)
ax.text(0.02, 0.98, f"all 2 Mb windows with 60-160 genes\nn={len(S)} windows\nNEAT1 window +{obs:.3f}, 1/55",
        transform=ax.transAxes, va="top", ha="left")
ax.set_xlabel("window mean - rest of genome"); ax.set_ylabel("number of windows")
ax.set_title("2 Mb sliding window\nnull across the genome"); ax.set_box_aspect(1)

ax = axes[2]
c11 = K[K.seq == "chr11"]
ax.scatter(c11.mid / 1e6, c11.y, s=4, color="0.4", alpha=.5)
w = 4e6
xs = np.arange(c11.mid.min(), c11.mid.max(), 5e5)
ys = [c11.y[(c11.mid > x - w / 2) & (c11.mid < x + w / 2)].mean() for x in xs]
ax.plot(xs / 1e6, ys, color="k", linewidth=1.4)
ax.axvline(NPOS / 1e6, color="k", linestyle=":", linewidth=1.2)
ax.axhline(0, color="0.6", linewidth=.8)
ax.text(0.02, 0.98, f"chr11, n={len(c11)} genes\n4 Mb running mean\ndotted = NEAT1",
        transform=ax.transAxes, va="top", ha="left")
ax.set_xlabel("chr11 position (Mb)"); ax.set_ylabel("CQN response"); ax.set_ylim(-1.2, 1.2)
ax.set_title("response along chr11\n"); ax.set_box_aspect(1)
fig.tight_layout(); fig.savefig("fig2_cis_decay.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ---------------- Fig 3: what survives, and how little it explains -------------
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.6))
def bs_(x, df=5):
    x = np.asarray(x, float); q = np.quantile(x, np.linspace(0, 1, df + 1)[1:-1])
    return np.column_stack([np.ones_like(x), x, x**2, x**3] + [np.clip(x - k, 0, None)**3 for k in q])
X = np.column_stack([bs_(Ks.gc_gene.values), bs_(Ks.log10_len_gene.values)[:, 1:], bs_(Ks.log10_len_exon.values)[:, 1:]])
yl = np.log2(cp + 1)
R = pd.DataFrame({s: yl[s].values - X @ np.linalg.lstsq(X, yl[s].values, rcond=None)[0] for s in samples}, index=cp.index)
is_lnc = (Ks.gene_type == "lncRNA").values
near = ((Ks.seq == "chr11") & ((Ks.start + Ks.end) / 2 - NPOS).abs().lt(2e6)).values
lab_bm = np.log10(Ks.baseMean.values)
P2 = []
for pair in itertools.combinations(samples, 2):
    rest = [s for s in samples if s not in pair]
    d = R[rest].mean(axis=1).values - R[list(pair)].mean(axis=1).values
    P2.append((set(pair) == set(negs), d[is_lnc].mean() - d[~is_lnc].mean(),
               d[near].mean() - d[~near].mean(), stats.spearmanr(lab_bm, d)[0]))
P2 = pd.DataFrame(P2, columns=["is_true", "lnc", "cis", "bm"])
for ax, col, lab, p in zip(axes[:2], ["lnc", "cis"],
                           ["lncRNA - other", "within 2 Mb of NEAT1 - rest"], ["1/66, p=0.015", "1/66, p=0.015"]):
    t = P2[P2.is_true][col].iloc[0]
    ax.hist(P2[col], bins=22, color="0.75", edgecolor="k", linewidth=.4)
    ax.axvline(t, color="k", linewidth=1.6)
    ax.text(0.02, 0.98, f"true split {t:+.3f}\nrank {p}", transform=ax.transAxes, va="top", ha="left")
    ax.set_xlabel(f"{lab} (CQN residual)"); ax.set_ylabel("number of splits")
    ax.set_title(f"permutation null\n{lab}"); ax.set_box_aspect(1)
ax = axes[2]
names = ["GC", "gene length", "abundance", "Alu / IRAlu", "3'UTR 6-mers", "RBP eCLIP\n(NONO/SFPQ/TDP43/FUS)",
         "nuclear:cyto", "lncRNA class", "cis <2 Mb"]
vals = [0.106, 0.091, 0.121, 0.30, 0.50, 0.40, 0.30, 0.015, 0.015]
colr = ["0.8" if v > 0.05 else "0.3" for v in vals]
ax.barh(range(len(names)), [-np.log10(v) for v in vals], color=colr, edgecolor="k", linewidth=.4)
ax.axvline(-np.log10(0.05), color="k", linestyle="--", linewidth=1)
ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=8)
ax.set_xlabel("-log10 p (sample-label permutation)")
ax.set_title("what survives the permutation\n(dashed = 0.05; floor is 1/66)")
ax.set_box_aspect(1)
fig.tight_layout(); fig.savefig("fig3_what_survives.png", dpi=160, bbox_inches="tight"); plt.close(fig)
print("figures written")
