"""Fig 7: intron vs exon response after region-level GC/length correction, the cis decay
of the intron response, and the 66-split permutation summary. House style (style.py)."""
import pandas as pd, numpy as np, itertools
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
K = pd.read_pickle("work/nascent_cqn.pkl"); P = pd.read_pickle("work/nascent_perm66.pkl")
NPOS = (65422774 + 65445540) / 2
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.6))
# (a) exon vs intron response, 2d hist
ax = axes[0]
h = ax.hist2d(K.yE, K.yI, bins=80, range=[[-2, 2], [-2, 2]], cmap="Blues", norm=matplotlib.colors.LogNorm())
ax.plot([-2, 2], [-2, 2], color="0.5", linewidth=0.8); ax.axhline(0, color="0.7", linewidth=0.6); ax.axvline(0, color="0.7", linewidth=0.6)
r = np.corrcoef(K.yE, K.yI)[0, 1]
ax.text(0.02, 0.98, f"n={len(K)}\nr={r:.2f}\nmedian exon {K.yE.median():+.2f}\nmedian intron {K.yI.median():+.2f}", transform=ax.transAxes, va="top", ha="left", fontsize=8)
ax.set_xlabel("exon response, log2 (probe / neg), GC-corrected"); ax.set_ylabel("intron response, same gene")
ax.set_title("counts_exon.tsv vs\ncounts_intron.tsv"); ax.set_box_aspect(1)
# (b) cis decay: intron vs exon
ax = axes[1]
rs = [0.5, 1, 2, 3, 5, 10, 20]
chr11 = K.seq == "chr11"; d = (K.mid - NPOS).abs()
rng = np.random.default_rng(0)
for col, c, lab in [("yI", "k", "intron"), ("yE", "0.6", "exon")]:
    m = []; lo = []; hi = []
    for r_ in rs:
        v = K.loc[chr11 & (d < r_ * 1e6), col].values
        bs = [rng.choice(v, len(v)).mean() for _ in range(1000)]
        m.append(v.mean()); lo.append(np.percentile(bs, 2.5)); hi.append(np.percentile(bs, 97.5))
    ax.errorbar(rs, m, yerr=[np.array(m) - lo, np.array(hi) - m], fmt="o-", color=c, capsize=3, markersize=4, linewidth=1, label=lab)
    ax.axhline(K[col].mean(), color=c, linestyle="--", linewidth=0.8)
ax.set_xscale("log"); ax.set_xlabel("radius around NEAT1 (Mb)"); ax.set_ylabel("mean corrected response, genes within radius")
ax.legend(frameon=False, fontsize=8); ax.set_title("cis decay is intronic\n(dashed: genome means)"); ax.set_box_aspect(1)
ax.text(0.02, 0.02, "n per radius: " + ", ".join(str(int((chr11 & (d < r_ * 1e6)).sum())) for r_ in rs), transform=ax.transAxes, fontsize=7, va="bottom")
# (c) 66-split summary
ax = axes[2]
labels = [("med_dIE", "intron - exon, genome median"), ("near_I", "intron, within 2 Mb of NEAT1"), ("near_E", "exon, within 2 Mb of NEAT1"),
          ("lncA_E", "lncRNA - protein-coding, exon,\nAlu-free genes"), ("lncA_I", "lncRNA - protein-coding, intron,\nAlu-free genes"), ("alu_E", "sense-Alu >=2 minus none, exon")]
t = P[P["true"]].iloc[0]
for i, (c, lab) in enumerate(labels):
    null = P.loc[~P["true"], c].values
    ax.scatter(null, np.full(len(null), i), s=8, color="0.6", alpha=0.6)
    ax.scatter([t[c]], [i], s=40, color="k", zorder=3)
    rank = int((P[c] >= t[c]).sum())
    ax.text(ax.get_xlim()[1] if False else max(P[c].max(), 0) + 0.02, i, f"{rank}/66", fontsize=8, va="center")
ax.set_yticks(range(len(labels))); ax.set_yticklabels([l for _, l in labels], fontsize=7); ax.invert_yaxis()
ax.axvline(0, color="0.7", linewidth=0.6); ax.set_xlabel("contrast (black = true split, grey = the other 65)")
ax.set_title("66-split permutation after\nGC/length correction"); ax.set_box_aspect(1)
fig.tight_layout(); fig.savefig("figures/fig7_nascent_cqn.png", dpi=160, bbox_inches="tight"); print("fig7 written")
