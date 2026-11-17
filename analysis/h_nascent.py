"""Does paraspeckle enrichment track nascent / chromatin-associated RNA?

Discriminator between two readings of the GC-negative, length-positive signal:
  (a) real sequence/compartment preference of the paraspeckle
  (b) the pulldown simply captures intron-rich nascent transcripts, and GC/length
      are proxies for intron content
"""
import pandas as pd, numpy as np
from scipy import stats

M = pd.read_pickle("master.pkl")
K = M[M.keep].copy()

def rr(x, covs):
    rx = stats.rankdata(x)
    X = np.column_stack([np.ones(len(x))] + [stats.rankdata(c) for c in covs])
    b, *_ = np.linalg.lstsq(X, rx, rcond=None)
    return rx - X @ b

def prho(a, b, covs):
    return np.corrcoef(rr(a, covs), rr(b, covs))[0, 1]

# --- intron content ---
K["intron_frac"] = K.len_intron / K.len_gene
K["n_introns"] = (K.n_exons - 1).clip(lower=0)
K["intronless"] = (K.n_introns == 0).astype(int)
print("== intron content vs LFC ==")
for c in ["intron_frac", "n_introns", "log10_len_intron", "log10_len_exon"]:
    print(f"  {c:20s} rho={stats.spearmanr(K[c], K.lfc)[0]:+.3f}  "
          f"even={stats.spearmanr(K[c], K.lfc_40even)[0]:+.3f}  odd={stats.spearmanr(K[c], K.lfc_40odd)[0]:+.3f}")

g = K.groupby("intronless").agg(n=("lfc", "size"), med=("lfc", "median"), gc=("gc_gene", "median"),
                                len=("log10_len_gene", "median"))
print(g.round(3).to_string())
sub = K[K.n_introns > 0]
print("  single-intron vs multi:",
      sub.groupby(pd.cut(sub.n_introns, [0, 1, 2, 5, 10, 20, 50, 400])).lfc.median().round(3).to_dict())

# Is length effect just intron length, or does exon length matter at fixed intron length?
cov = [K.gc_gene.values, K.log10_baseMean.values]
print("  partial | GC,expr : log10_len_intron %.3f   log10_len_exon %.3f   intron_frac %.3f"
      % (prho(K.log10_len_intron.values, K.lfc.values, cov),
         prho(K.log10_len_exon.values, K.lfc.values, cov),
         prho(K.intron_frac.values, K.lfc.values, cov)))

# within intronless genes only, does GC still predict?
il = K[K.intronless == 1]
print("  intronless only (n=%d): rho(GC, LFC)=%.3f  rho(len_gene, LFC)=%.3f"
      % (len(il), stats.spearmanr(il.gc_gene, il.lfc)[0], stats.spearmanr(il.log10_len_gene, il.lfc)[0]))
ml = K[K.n_introns >= 5]
print("  >=5 introns (n=%d):    rho(GC, LFC)=%.3f  rho(len_gene, LFC)=%.3f"
      % (len(ml), stats.spearmanr(ml.gc_gene, ml.lfc)[0], stats.spearmanr(ml.log10_len_gene, ml.lfc)[0]))

# --- histone genes: intronless, non-polyadenylated, replication-dependent ---
hist = K[K.gene_name.str.match(r"^H[1-4][A-C]?C\d+$|^H2[AB][A-Z]\d*$", na=False)]
print("\n== histone genes ==", len(hist), "median LFC", round(hist.lfc.median(), 3) if len(hist) else None)

# --- snoRNA/scaRNA hosts and other classes ---
for pat, lab in [(r"^SNHG", "SNHG host"), (r"^MIR", "MIR"), (r"^LINC", "LINC"), (r"-AS1$", "antisense"), (r"-DT$", "divergent")]:
    s = K[K.gene_name.str.contains(pat, na=False, regex=True)]
    if len(s) > 15:
        print(f"  {lab:12s} n={len(s):5d} median LFC={s.lfc.median():+.3f} (all {K.lfc.median():+.3f})")
