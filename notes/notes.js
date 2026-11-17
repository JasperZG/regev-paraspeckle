const fs = require("fs");
const d = require("docx");
const { Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
        WidthType, ShadingType, BorderStyle, ImageRun, AlignmentType, LevelFormat,
        ExternalHyperlink } = d;

const FONT = "Arial";
const SZ = 22; // 11pt

function run(t, o = {}) { return new TextRun({ text: t, font: FONT, size: SZ, ...o }); }

// bullet level helper: Google-Docs-like glyph progression
function b(level, children, opts = {}) {
  return new Paragraph({
    numbering: { reference: "gbul", level },
    children: Array.isArray(children) ? children : [run(children)],
    spacing: { after: 0, line: 276 },
    ...opts,
  });
}
function p(children, opts = {}) {
  return new Paragraph({
    children: Array.isArray(children) ? children : [run(children)],
    spacing: { after: 120, line: 276 },
    ...opts,
  });
}
function img(path, w, h) {
  return new Paragraph({
    children: [new ImageRun({ type: "png", data: fs.readFileSync(path), transformation: { width: w, height: h } })],
    spacing: { before: 160, after: 60 },
  });
}
function cap(parts) {
  return new Paragraph({ children: parts, spacing: { after: 220, line: 276 } });
}
const BRD = { style: BorderStyle.SINGLE, size: 4, color: "000000" };
function tbl(rows, widths) {
  return new Table({
    columnWidths: widths,
    width: { size: widths.reduce((a, c) => a + c, 0), type: WidthType.DXA },
    rows: rows.map(cells => new TableRow({
      children: cells.map((c, i) => new TableCell({
        width: { size: widths[i], type: WidthType.DXA },
        margins: { top: 60, bottom: 60, left: 100, right: 100 },
        borders: { top: BRD, bottom: BRD, left: BRD, right: BRD },
        children: [new Paragraph({ children: [run(c.t, { bold: !!c.b })], spacing: { after: 0, line: 276 } })],
      })),
    })),
  });
}

const numbering = {
  config: [{
    reference: "gbul",
    levels: [0, 1, 2, 3, 4, 5].map(i => ({
      level: i,
      format: LevelFormat.BULLET,
      text: ["●", "○", "■", "●", "○", "■"][i],
      alignment: AlignmentType.LEFT,
      style: { paragraph: { indent: { left: 720 * (i + 1), hanging: 360 } },
               run: { font: FONT, size: SZ } },
    })),
  }],
};

const kids = [];
const P = x => kids.push(x);

// ============================ 1. framing ============================
P(b(0, [run("What determines paraspeckle enrichment?", { bold: true })]));
P(b(1, "Goal from Oded: work out what about a gene sets its log2FoldChange in the NEAT1 pulldown. Started from the simplest hypotheses and worked down the list."));
P(b(1, "Data is all 12 libraries now, not just the two DE tables: NEAT1_input.xlsx (gene counts, 67,928 x 12) + ColData_NEAT1.xlsx. 6 conditions x 2 reps: neg, 20old, 40old, 40odd, 40even, 80. neg is the no-probe control."));
P(b(1, "Re-ran DESeq2 myself in pydeseq2 off the raw counts so I could change the design. My MLE log2FoldChange reproduces the given CSVs at r = 1.000 on the baseMean > 10 genes, so the tables and my pipeline agree and anything new is coming from the model, not the input."));
P(b(1, "Main response I use below = all 10 probe libraries vs the 2 neg, apeglm-shrunk, baseMean > 10, primary chromosomes, chrM and NEAT1 itself dropped. n = 12,235."));
P(b(2, "Sanity: NEAT1 itself is +3.61 (top gene in the whole table), so the pulldown definitely worked. MT-ATP6 and MTATP6P1 are the most depleted, which is what you want from a nuclear body pulldown."));
P(b(2, "Annotation: the 67,928 Geneids are exactly GENCODE v39 chr_patch_hapl_scaff. v38 has 67,049 and v40 has 67,939, and v42 onward reorders chr1, so v39 is the only release that lines up 1:1 and in order. Genome is GRCh38.p13 to match. Checked GAPDH gene-span GC = 0.599793 against the UCSC API, identical."));

// ============================ 2. simple hypotheses ============================
P(p([run("H1-H4, the simple stuff (abundance, length, GC, biotype)", { bold: true })]));
P(b(0, "Spearman rho against the shrunk LFC, and separately against 40even and 40odd so I can see if it replicates:"));
P(new Paragraph({ spacing: { before: 60, after: 120 } }));
P(tbl([
  [{ t: "feature", b: true }, { t: "rho vs pooled LFC", b: true }, { t: "40even", b: true }, { t: "40odd", b: true }],
  [{ t: "GC, gene body" }, { t: "-0.381" }, { t: "-0.377" }, { t: "-0.326" }],
  [{ t: "GC, 100 kb around gene" }, { t: "-0.360" }, { t: "-0.357" }, { t: "-0.317" }],
  [{ t: "GC, exons" }, { t: "-0.347" }, { t: "-0.335" }, { t: "-0.305" }],
  [{ t: "log10 intron length" }, { t: "+0.298" }, { t: "+0.331" }, { t: "+0.192" }],
  [{ t: "log10 gene length" }, { t: "+0.295" }, { t: "+0.326" }, { t: "+0.192" }],
  [{ t: "log10 exonic length" }, { t: "+0.112" }, { t: "+0.113" }, { t: "+0.083" }],
  [{ t: "log10 baseMean" }, { t: "+0.028" }, { t: "+0.011" }, { t: "+0.038" }],
  [{ t: "number of exons" }, { t: "+0.046" }, { t: "+0.044" }, { t: "+0.037" }],
], [3200, 2100, 1400, 1400]));
P(new Paragraph({ spacing: { after: 100 } }));
P(b(0, "So GC looks like the whole answer. rho = -0.38, replicates in both halves, and it survives everything I throw at it: partial rho | log baseMean is still -0.385, it holds inside intronless genes only (-0.312, n=440) and inside genes with 5+ introns (-0.403), so it isn't an intron-content proxy."));
P(b(0, "Length is the mirror image and it is intron length specifically. Partial | GC and expression: intron length +0.152, exonic length +0.048. GC and length are themselves correlated at rho = -0.512 (isochores), and gc_exon vs gc_win100k is +0.841, so these are all one axis."));
P(b(0, "lncRNA sits above protein-coding, median LFC +0.044 vs -0.001, MWU p = 2e-41."));
P(b(0, "At this point I had written down \"GC-poor, long, intron-rich genes go to paraspeckles\" and was about to go look for a mechanism. Then I checked whether the GC axis is even real."));

// ============================ 3. the permutation ============================
P(p([run("The GC result does not survive a label permutation", { bold: true })]));
P(b(0, "The worry: GC bias is the single most common library-prep artifact, and it is a per-library property, not a per-condition one. With only 2 neg libraries against 10 probe libraries, any library-level GC drift between those two groups becomes a \"biological\" GC effect."));
P(b(0, "Per-library GC profile (log2 of that library against the mean of all 12, rho vs gene GC). These should be ~0 if there is no GC drift:"));
P(b(1, "neg_Rep1 +0.416, neg_Rep2 +0.253, 20old_Rep1 -0.302, 20old_Rep2 -0.431, 40old_Rep1 +0.353, 40old_Rep2 -0.069, 40odd_Rep1 -0.222, 40odd_Rep2 +0.194, 40even_Rep1 -0.216, 40even_Rep2 -0.300, 80_Rep1 +0.165, 80_Rep2 +0.030."));
P(b(1, "They are all over the place, and both negs happen to land on the GC-high side. That is the entire \"biology\"."));
P(b(0, "Null contrasts, where there should be no signal at all:"));
P(b(1, "neg_Rep1 vs neg_Rep2: rho(GC) = +0.273. Two replicates of the same condition, and the GC slope is already 3/4 of the real effect."));
P(b(1, "80 vs 40even, both pulldowns: rho(GC) = +0.295."));
P(b(0, "Proper test. There are C(12,2) = 66 ways to pick which 2 libraries are \"the control\". Compute rho(GC, log ratio) for every one of them and see where the true split lands."));
P(b(1, "True split (neg vs rest) rho = -0.366. Rank 16 of 66 by |rho|. One-sided p = 0.106. Null SD 0.279, range -0.606 to +0.537."));
P(b(1, "Length is the same: true +0.281, rank 15/66, p = 0.091."));
P(b(1, "Restricting the null to the 6 splits that pair the two replicates of one condition (the like-for-like nulls) does not rescue it: 20old gives +0.429, 40even +0.324, 40old -0.178, 80 -0.114, 40odd +0.011. The true split is the most negative but it is inside the spread."));
P(b(0, [run("So the strongest-looking determinant in the whole dataset is not separable from library noise. ", { bold: true }),
        run("Not that GC is proven irrelevant, just that with n = 2 controls this design cannot tell. Anything correlated with GC or length inherits the problem, and that includes Alu content, which is why I did not chase IRAlu first.")]));
P(img("fig1_gc_is_technical.png", 630, 220));
P(cap([run("Fig 1. GC content looks like the answer and then does not survive label permutation. ", { bold: true }),
       run("(a) Shrunk log2FoldChange (10 probe libraries vs 2 neg) against gene-body GC fraction, r = -0.350, n = 11,931. (b) The same statistic for all 66 ways of choosing which two libraries are the control. The true split (vertical line) is rank 16/66, one-sided p = 0.106. (c) After conditional-quantile normalisation on GC and length, the response measured in batch 1 against batch 2, r = 0.617.")]));

// ============================ 4. CQN ============================
P(p([run("Fix: CQN, then redo everything", { bold: true })]));
P(b(0, "Took out the per-library GC/length curves the way cqn/EDASeq do: for each library fit log2 CPM ~ spline(GC) + spline(log gene length) + spline(log exonic length) and keep the residual. Response is then mean(residual over 10 probe) - mean(residual over 2 neg). Residual GC association per library drops to about -0.02."));
P(b(0, "The response gets more reproducible, not less, which is the point:"));
P(b(1, "Rep A libraries vs Rep B libraries, fully independent sets: r = 0.697."));
P(b(1, "Batch 1 (20old+40odd) vs batch 2 (40old+40even+80): r = 0.617, so reliability by Spearman-Brown is about 0.76. That is the ceiling on R2 for any gene feature."));
P(b(0, "Found a batch split while doing this that matters for the design. PCA of the libraries puts {40old, 40even, 80} on one side and {neg, 20old, 40odd} on the other (PC1 = 36% of variance), and hierarchical clustering agrees."));
P(b(1, "The negs sit with 20old and 40odd. So 20old/40odd vs neg is a within-batch contrast, and 40old/40even/80 vs neg is confounded with the batch axis. I carry both separately below."));
P(b(1, "Not probe identity: 40old is the older design and still groups with 40even. 40even and 40odd are disjoint halves of the same 80 and land on opposite sides. So this is handling, not oligos."));
P(b(1, "Chased the probe-cross-hybridisation idea anyway since it would have been the cleanest artifact: 208 genes show an \"even-oligo\" pattern (up in 40even and 80, flat in 40odd), e.g. TUB at 258 CPM in 40even and 52 in neg, and 33 show the odd version, e.g. PSRC1. But they carry no more NEAT1 sequence than controls (shared 16-mers with NEAT1 sense: 46% of even-only genes vs 42% of flat controls), so they are not annealing to a NEAT1 oligo. It is the batch."));

// ============================ 5. what died ============================
P(p([run("Hypotheses that died", { bold: true })]));
P(b(0, [run("Nuclear retention. ", { bold: true }), run("ENCODE HeLa-S3 nuclear vs cytosolic polyA RNA (ENCSR000CPQ/CPP, 2 reps each, nuc:cyt reproducible at r = 0.811). rho(response, log2 nuc:cyt) = -0.032 raw, +0.051 partialled on GC/length/expression. MALAT1 sits at nuc:cyt = 5.6 and has LFC -0.16. Being nuclear does not get you into a paraspeckle.")]));
P(b(0, [run("Inverted Alus in the 3'UTR, the textbook model. ", { bold: true }), run("Built it properly off RepeatMasker: Alu counts split by orientation relative to the gene, sense x antisense pair counts, inverted pairs within 2 kb, per region. Raw rho for 3'UTR IRAlu is +0.057, and genes with a 3'UTR IRAlu pair do have a higher median (+0.031 vs +0.003, p = 2e-10) - but after CQN it falls to +0.045 and it tracks the GC/length axis. Only 742 of 12,235 genes have one at all, so it cannot be the general rule here even if it is real for individual transcripts like NICN1 or CTN-RNA.")]));
P(b(0, [run("A sequence motif in the 3'UTR. ", { bold: true }), run("Scanned all 4,096 6-mers as density per kb in 3'UTRs (n = 10,221) and in whole canonical transcripts, scoring each against batch 1 and batch 2 separately. Best consistent 6-mer is GTTATC at rho = +0.040. Nothing. In particular no C-rich/HNRNPK signal (the Jacq et al. 15-mer has a 5' C stretch), no UG repeat for TDP-43, no AG-rich for the West et al. shell.")]));
P(b(0, [run("Sequence similarity to NEAT1. ", { bold: true }), run("Exact 25-mers shared with NEAT1 in sense orientation (probes would anneal) vs antisense (would base-pair with the lncRNA). 1,376 genes share a sense 25-mer, rho = +0.077; antisense rho = +0.036. The small sense > antisense gap is consistent with a little probe cross-reaction but it is tiny and it is Alu-driven - NEAT1 itself is full of L1 and Alu.")]));
P(b(0, [run("RBP binding. ", { bold: true }), run("ENCODE eCLIP peaks for NONO (K562), SFPQ, FUS x2, TARDBP x2, HNRNPK x2, MATR3 x2, HNRNPU x2, EWSR1, TAF15, AKAP8L, ZC3H11A, plus PTBP1/U2AF2/SRSF1/SUGP2 as non-paraspeckle controls. Peak density per kb, partialled on length and expression. Every single one is negative, the strongest being TARDBP_K562 exons at -0.067, and the paraspeckle proteins are indistinguishable from the controls (median partial rho -0.028 vs -0.044). Genes bound by abundant RBPs are slightly depleted, which is probably just that well-covered transcripts are the ones that get peaks.")]));
P(b(0, [run("Genomic position / compartment. ", { bold: true }), run("Distance to telomere looked good raw (rho = +0.172) but goes to +0.005 once GC and length are in the model, so it was the isochore again. Decisive version: spatial autocorrelation of the response between gene pairs on the same chromosome, binned by separation from 0.1 to 50 Mb, against a within-chromosome permutation. Real r never exceeds +0.011 at any separation and the permuted baseline is +0.00x. Neighbouring genes do not share a response.")]));
P(b(1, "That one negative kills the whole class at once - LADs, A/B compartment, speckle proximity, replication timing are all megabase-scale, so if any of them drove this the autocorrelation would be obvious."));
P(b(0, [run("Abundance. ", { bold: true }), run("After CQN there is a monotone trend, median response +0.041 in the lowest expression decile down to -0.049 in the highest. Ran it through the same 66-split permutation: rho = -0.081, rank 21/66, p = 0.121. Not real either. Glad I checked, I nearly wrote this one up.")]));

// ============================ 6. what survives ============================
P(p([run("What survives the permutation", { bold: true })]));
P(b(0, "Ran the 66-split test on everything, so the bar is the same for all of them. Floor on p is 1/66 = 0.015."));
P(new Paragraph({ spacing: { before: 60, after: 120 } }));
P(tbl([
  [{ t: "effect", b: true }, { t: "true split", b: true }, { t: "rank", b: true }, { t: "p", b: true }],
  [{ t: "lncRNA minus everything else" }, { t: "+0.122" }, { t: "1/66" }, { t: "0.015" }],
  [{ t: "antisense/divergent minus rest" }, { t: "+0.137" }, { t: "1/66" }, { t: "0.015" }],
  [{ t: "within 2 Mb of NEAT1 minus rest" }, { t: "+0.103" }, { t: "1/66" }, { t: "0.015" }],
  [{ t: "GC" }, { t: "-0.366" }, { t: "16/66" }, { t: "0.106" }],
  [{ t: "gene length" }, { t: "+0.281" }, { t: "15/66" }, { t: "0.091" }],
  [{ t: "abundance" }, { t: "-0.081" }, { t: "21/66" }, { t: "0.121" }],
], [3900, 1500, 1200, 1200]));
P(new Paragraph({ spacing: { after: 100 } }));
P(b(0, "Three things, and they are the three maximally significant results the design can produce."));
P(b(1, "lncRNA class, +0.124 mean difference. Independent of the NEAT1 neighbourhood (+0.126 with chr11 excluded). Among the 212 genes that replicate in both batches, 26% are lncRNA against 9% background and 13% carry an -AS1/-DT/LINC name against 4% background."));
P(b(1, "The antisense/divergent effect is mostly inside the lncRNA effect - within lncRNAs the extra bump is +0.041, p = 0.18 - so I treat it as one thing."));
P(b(1, "Cis-proximity to the NEAT1 locus, +0.100 overall and +0.127 if you use protein-coding genes only, so it is not the lncRNA effect in disguise."));
P(b(0, "The cis one is the result I believe most, because it has a shape and a proper spatial null:"));
P(b(1, "Decay with radius: +0.173 within 0.5 Mb (n=30), +0.099 within 2 Mb (n=93), +0.060 within 3 Mb, +0.030 within 10 Mb, genome mean 0.000. 95% bootstrap CI excludes 0 at 0.5, 2 and 3 Mb."));
P(b(1, "Slid a 2 Mb window across the whole genome and kept the 55 windows holding 60-160 genes. The NEAT1 window is the top one. Only 1.8% of windows beat it, and the next three best are also chr11 at 66.2, 67.2 and 64.2 Mb - the same neighbourhood."));
P(b(1, "Leave-one-gene-out on the 93 genes: effect stays between +0.081 and +0.107, so it is not two or three outliers."));
P(b(1, "Partial replication of Yap et al. 2022 HyPro-seq: PC is at the 96th percentile and CCDC57 at the 88th; SF3B2, which was their negative, is at the 15th. NAA40 and EML3 do not come along."));
P(b(0, [run("The honest problem with the cis result, and why I still believe it. ", { bold: true }),
        run("It is not visible in the delivered DE tables at all - median LFC near NEAT1 is -0.03 for 80, -0.18 for 40even, -0.21 for 20old, i.e. the wrong way. It only appears after CQN. The reason is that chr11q13 is one of the most GC-rich neighbourhoods in the genome (median gene GC 0.557 against 0.436 genome-wide), so the technical GC axis buries it. That makes the result entangled with the correction, and GC over-correction would produce exactly this, so it needs a GC-matched control, not a genome-wide one.")]));
P(b(1, "Across the 48 windows with median GC > 0.50, the NEAT1 window is rank 1/48, and chr11 66.2 Mb is rank 2. The other GC-rich windows (chr16, chr17, chr19q13, chr8q24) are scattered around zero: window GC vs window effect is only rho = +0.151, nowhere near enough to make +0.10."));
P(b(1, "Gene-by-gene matching: pair each of the 93 window genes with a non-window gene of the same GC (within 0.015), 2000 draws. Difference +0.103, 95% CI [+0.043, +0.162], p = 0.002. Matching on GC and length and expression together: +0.106, CI [+0.045, +0.162], p = 0.001."));
P(b(1, "Inside the window GC runs the other way (rho = +0.116 against y), which it would not if this were leftover GC trend."));
P(b(0, "Checked the normalisation isn't doing the work. Swapped the splines for a crude median-centre within GC x length bins and re-ran the 66-split permutation at four settings (8x8, 12x12, 10x10 bins; min count 2, 5, 10). cis ranks 1/66 and lncRNA ranks 1/66 in every one; GC ranks 32-66/66 in every one. With no correction at all, lncRNA is still 1/66 but cis is 50/66, which is the masking above."));
P(img("fig2_cis_decay.png", 630, 220));
P(cap([run("Fig 2. The response decays with genomic distance from the NEAT1 locus. ", { bold: true }),
       run("(a) Mean CQN response for genes within r Mb of NEAT1 (chr11:65.43 Mb), 95% bootstrap intervals, dashed line is the genome mean. (b) The same 2 Mb window statistic computed at every position genome-wide, keeping windows with 60-160 genes; the NEAT1 window (vertical line) is 1 of 55. (c) Response along chr11 with a 4 Mb running mean; dotted line is NEAT1.")]));
P(img("fig3_what_survives.png", 630, 196));
P(cap([run("Fig 3. Only three effects clear the sample-label permutation. ", { bold: true }),
       run("(a) lncRNA minus non-lncRNA mean response, for all 66 two-vs-ten splits; the true split is the most extreme of the 66. (b) The same for genes within 2 Mb of NEAT1. (c) Every hypothesis on one axis; dashed line is p = 0.05 and the attainable floor is 1/66 = 0.015.")]));

// ============================ 7. hypothesis ============================
P(p([run("Where I land", { bold: true })]));
P(b(0, [run("Hypothesis: what gets labelled is set by where the transcript is made and how long it stays on chromatin, not by anything in the mature RNA.", { bold: true })]));
P(b(1, "Paraspeckles nucleate on nascent NEAT1 at the NEAT1 locus itself. If the labelling radius is short - ARTR-seq is about 45 nm and HyPro is comparable - then what you capture is whatever is physically there: NEAT1, nascent transcripts of neighbouring genes on chr11, and the chromatin-tethered transcripts that linger near their own locus. Those are disproportionately lncRNAs and divergent/antisense transcripts, which is exactly the class that comes out."));
P(b(1, "It explains all three survivors with one mechanism, and it explains the negatives too: there is no motif, no IRAlu requirement, no NONO/SFPQ/TDP-43 binding signature, because nothing is recognising the transcript. It is proximity."));
P(b(1, "It also explains the shape - a decay over ~2-3 Mb centred on the locus, with no broader spatial autocorrelation, is what you get from one nucleation site rather than from a compartment."));
P(b(0, [run("Caveat I cannot get around with this design. ", { bold: true }),
        run("The response is highly reproducible (reliability 0.76) but the three surviving features together give a leave-one-chromosome-out CV R2 of only 0.022. So 97% of the reproducible signal is gene-specific and unexplained by anything I measured. Some of that is real paraspeckle biology on individual transcripts. But some of it is almost certainly the capture chemistry: the neg control is no-probe, so it controls for the beads but not for hybridisation capture itself, and any transcript that sticks to an oligo-loaded bead scores as enriched reproducibly forever. Yap et al. used a scrambled probe set for exactly this reason.")]));

// ============================ 8. next ============================
P(p([run("What I want to do next", { bold: true })]));
P(b(0, "Ask Fei lab for a scrambled / non-targeting probe control, or a probe set against an unrelated abundant nuclear lncRNA. Without it we cannot separate \"near a paraspeckle\" from \"sticks to a capture oligo\", and that is currently the largest term."));
P(b(0, "More neg replicates. Two is why the GC question is unanswerable - with 4-6 negs the permutation null tightens enough that a real GC effect of rho = -0.37 would clear it easily."));
P(b(0, "Intron vs exon counts per gene, the Wu et al. Group A/B split. This is the single best test of the hypothesis above: if the signal is nascent and chromatin-proximal, intronic enrichment should beat exonic enrichment, and it should do so hardest near the NEAT1 locus. I have the three BAMs but they are 477-590 MB and over the staging limit, so I need them split or a featureCounts run with -O on an intron SAF."));
P(b(0, "Where the probes sit on NEAT1_2. If the two halves tile the middle vs the two ends, 40even and 40odd are sampling the core and the shell separately (West et al. core-shell, 5' and 3' ends in the shell, middle in the core, about 110 nm apart), which would make the even/odd difference interpretable rather than a nuisance. Also need to know whether any probe falls in the first 3.7 kb shared with NEAT1_1."));
P(b(0, "Hi-C contact frequency with chr11:65.43 Mb as a predictor, in place of linear distance. Linear distance already works out to 2-3 Mb; contacts would test the trans part and would be the direct version of the hypothesis."));
P(b(0, "Reporter: same promoter, same locus, swap the 3'UTR between a top gene (EMSY-DT, CCDC88B, TMEM9B-AS1) and a flat control, run through the same pulldown. If localisation is positional and not sequence-encoded, swapping the UTR should do nothing - which is the prediction that would actually falsify this."));

// ============================ 9. files ============================
P(p([run("Files", { bold: true })]));
P(b(0, "deseq_run.py - pydeseq2, per-condition and pooled probe-vs-neg, with apeglm shrinkage"));
P(b(0, "cqn.py - per-library GC/length spline removal, builds the response"));
P(b(0, "h_gcbias2.py, abund_perm.py - the 66-split label permutation"));
P(b(0, "repeats.py, neat1_kmers.py, motifs.py, eclip_fetch.py - the feature sets that failed"));
P(b(0, "cis.py, cis_control.py, spatial.py - the positional tests"));
P(b(0, "figs.py, style.py - figures above, same house style"));
P(b(0, "top_replicating_genes.tsv - the 212 genes that replicate in both batches"));

// ============================ 10. second pass: the raw data ============================
P(p([run("Second pass (16 Nov): went back to the BAMs and the probe table", { bold: true })]));
P(b(0, "Everything above was done off the two DE tables. Jingxuan's folder also has the 12 deduplicated STAR BAMs, the DESeq2 inputs, a box note describing the five probe sets, and two probe-position sheets. Went through all of it before adding any new hypothesis, because several of the open questions in \"what I want to do next\" are answerable straight from these files."));
P(b(1, "Box note (10 Sep): 20old = 20 probes over the whole NEAT1, 40old = 40 over the whole thing, 40odd = the 40 probes in the odd-numbered regions, 40even = the even regions, 80 = 40odd + 40even. So odd/even is not alternating oligos, it is alternating ~2 kb regions of NEAT1_2."));
P(b(1, "Probe sheet (2026_08 batch): 80 probes in 10 regions x 8, each a 20-nt target-complementary part plus a common 35-nt 3' tail. Mapped every 20-mer back onto NEAT1 (chr11:65,422,774-65,445,540): all 80 hit exactly once, sense, between 3,846 and 21,577 nt. Nothing in the first 3.7 kb, so NEAT1_1 is never targeted directly, and nothing in the last ~1.2 kb."));
P(b(1, "JXJ 401-420 (16 probes) = 20old. 40old = those + JXJ 535-576 (24). JXJ 783-822 (40) are the new design. Odd regions 1,3,5,7,9 hold 40odd, even regions hold 40even. Two probes (572, 810) are annotated as really sitting in the neighbouring region; used the corrected assignment."));
P(b(1, "So every probe has one of six set-membership vectors over (20old, 40old, 40odd, 40even, 80): {1,1,1,0,1}, {1,1,0,1,1}, {0,1,1,0,1}, {0,1,0,1,1}, {0,0,1,0,1}, {0,0,0,1,1}. This turns out to be the key to the whole batch problem, below."));
P(b(1, "Libraries: single-end, 100 nt after trimming, UMI-deduplicated, rRNA pre-filtered, unique mappers only (STAR --outFilterMultimapNmax 1), 3.9-8.0 M mapped reads after dedup. Stranded, read = transcript strand (GAPDH 8,794 fwd : 0 rev)."));
P(b(1, "Cell line: the raw reads contain HPV18 E6/E7 (1,066 sense reads per 4 M in neg_Rep1, zero antisense), no chrY, no XIST, CDKN2A high. HeLa. Same line as Yap et al. 2022, which helps later."));
P(b(1, "The delivered counts are whole-gene-span counts. featureCounts -s 1 on the GENCODE v39 gene spans reproduces NEAT1_input.xlsx exactly (67,928 x 12, zero differences), whereas exon-only counting gives 5.2 M assigned reads for neg_Rep1 instead of 7.0 M. So every LFC in the first pass mixes intronic and exonic reads, 27-33% of assigned reads are intronic, and that is where the gene-length / intron-length association came from. It also means I can now do the intron-vs-exon test directly."));

// ============================ 11. NEAT1 itself ============================
P(p([run("NEAT1 itself: the capture works on RNA fragments at the probe sites", { bold: true })]));
P(b(0, "Coverage along NEAT1 in 100-bp bins, CPM-normalised, all 12 libraries."));
P(b(1, "neg: 105 CPM/bp over the first 3.7 kb (the NEAT1_1 region), 15 over 3.8-21.6 kb, 31 over the 3' end. About 6x more NEAT1_1 than NEAT1_2 in HeLa total RNA, which is the known picture."));
P(b(1, "Any probe set: flat, ~300-370 CPM/bp over the whole 22.7 kb. Body enrichment log2 4.5-6.5, but the first 3.7 kb only log2 1.3, i.e. the 5' region comes along only as the 5' end of NEAT1_2 molecules. NEAT1_1 is not co-captured. If this were a labelling of everything inside a paraspeckle, or a pulldown of intact paraspeckle contents, the NEAT1_1 that is in them would come along too."));
P(b(1, "40odd vs 40even along NEAT1: log2(40odd/40even) oscillates with region parity, +0.45, -0.46, +0.77, -0.64, +1.10, -1.22, +0.29, -1.42, +1.50, -1.41 for regions 1-10. Perfect alternation, 2-3 fold. 20old against 80 does the same: its enrichment peaks over 8-12 kb and 17-20 kb, which is exactly where the 16 JXJ 4xx probes sit."));
P(b(1, "Toy model: fragment the RNA to length L, keep fragments that contain a probe site plus 10% background, predict the odd/even profile. r = 0.77 at L = 200-400 nt, 0.63 at 1 kb, 0.10 at 2 kb. The fragments are a few hundred nt."));
P(b(0, [run("What this says about the method. ", { bold: true }), run("The RNA is fragmented to a few hundred nt and fragments carrying a probe site are what gets pulled down, with sub-kb resolution along the RNA. That is hybridisation capture of fragmented RNA (ChIRP/RAP-like), not a proximity label with a radius. So a non-NEAT1 fragment has exactly three ways in: (a) it anneals to one of the 80 oligos itself, (b) it is crosslinked or base-paired to a NEAT1 fragment, (c) it sticks to the beads. neg controls for (c) only. (a) and (b) are not separable by neg, but (a) has a signature: it is probe-specific.")]));
P(img("fig4_neat1_profile.png", 630, 220));
P(cap([run("Fig 4. The pulldown captures NEAT1_2 fragments at the probe sites. ", { bold: true }),
       run("(a) Coverage along NEAT1 (100-bp bins, CPM per bp) in neg (grey) and 80 (black); dotted line is the NEAT1_1 3' end at 3.7 kb. (b) log2(40odd / 40even) along NEAT1 with the odd-region probes ticked on top and the even-region probes below. (c) log2(probe set / neg), 500-bp running mean, for 20old, 40old and 80; ticks are the 16 JXJ 4xx probes that make up 20old.")]));

// ============================ 12. the batch is an Alu probe ============================
P(p([run("The \"batch\" was never a batch: JXJ 572 is an Alu probe", { bold: true })]));
P(b(0, "Membership logic. If a gene is captured because its fragment anneals to oligo P, it is enriched in exactly the sets that contain P. A gene that is near a paraspeckle is enriched in all five sets, give or take efficiency (NEAT1 itself is +3.3 to +3.7 at the exon level in every set). Fit each gene's five exon-level log2(set/neg) values to the six membership vectors; keep genes with contrast > 1 log2, R2 > 0.8. 86 genes. Shuffling the five set labels within gene gives 0.7 genes per shuffle."));
P(b(1, "79 are {40old, 40even, 80}, i.e. the JXJ 5xx even-region class (14 probes), which is also exactly the batch split from the first pass. 7 are not: PSRC1, MCU, KCNT2, PIEZO2 are {40odd, 80} only (PSRC1 is +3.3 in 40odd, +2.9 in 80, +0.1 to +0.4 in the other three), TUB is {40even, 80} only (+3.2 / +2.8 vs -0.1), MIB1 {40old, 40odd, 80}, MRPL51 {20old, 40old, 40even, 80}. No handling batch can produce \"up in 40odd and 80, flat in 20old, 40old and 40even\". Those are oligos."));
P(b(1, "They are exon effects: for the 86, exon log2 in the sets that contain the matching probes is +1.57 (median), in the other sets -0.07; intron +0.40 vs +0.06. Mature RNA is being captured, which is what annealing looks like."));
P(b(0, [run("Found the sequences. ", { bold: true }), run("Pulled the read pile-ups in the matching set and scanned them against all 80 target sites. TUB: 558 reads in 40even vs 6 in 40odd on a 200-bp bin at chr11:8,104,844 that carries the JXJ 806 target site with 2 mismatches in 20 (JXJ 806 is a new even-region probe, so {40even, 80} is the right vector). PSRC1: JXJ 818 site, 4 mismatches, 99 vs 1. MCU: JXJ 815, 3 mismatches, 130 vs 17. And every one of the batch-2 genes I looked at is the same probe: MAP3K13 chr3:185,292,341 JXJ 572 site with 1 mismatch (114 vs 1 reads), SIRT7 JXJ 572 1 mismatch (158 vs 3), DHRS4-AS1 JXJ 572 2 mismatches (135 vs 2).")]));
P(b(0, [run("JXJ 572 sits in an Alu. ", { bold: true }), run("RepeatMasker on NEAT1: the JXJ 572 target (position 16,760) is inside an AluSz6 in the same orientation as NEAT1. The target AAATTAGTCAGGCATGGTGG is 3 mismatches from the Alu consensus AAATTAGCCGGGCGTGGTGG (antisense: 10). No other probe is within 6 mismatches of the consensus; 17 probes sit in old L1 fragments and JXJ 814 in a MIR, none of which are young enough to cross-hybridise. So one oligo in the 40old/40even/80 sets anneals to every transcript carrying a sense-oriented Alu, and the odds of a transcript carrying one scale with its Alu count.")]));
P(b(1, "Of the 79 batch-2 signature genes, 78 (99%) carry a sense exonic Alu against 43% of all tested genes and 14% of the other 7 signature genes. Median 3 sense exonic Alus vs 0."));
P(b(1, "Genome-wide it is a dose curve. Mean exon log2({40old,40even,80} / {20old,40odd}) by number of sense exonic Alus: 0: -0.044 (n=5,838), 1: -0.006, 2-3: +0.060, 4+: +0.125 (n=654). Antisense-only Alus: flat to negative at every count. Of all C(10,4) = 210 ways of choosing which 4 probe libraries are \"batch 1\", the true choice gives the largest sense-Alu contrast (rank 1/210); the antisense contrast ranks 188/210."));
P(b(1, "And the batch axis vanishes once Alu-carrying genes are removed. Correlation structure of the per-library log ratios: all genes, mean r within batch +0.11, between batches -0.30, and clustering splits {20old, 40odd} from {40old, 40even, 80}. Genes with no sense exonic Alu: -0.04 / -0.17 and the clustering collapses to replicate pairs. Genes with 2+ sense exonic Alus: +0.36 / -0.51. The \"batch\" is sets-with-JXJ-572 vs sets-without."));
P(b(0, [run("Genome-wide scan, to close it. ", { bold: true }), run("Scanned GRCh38 for all 80 target sites in both orientations at up to 3 mismatches (12-mer seed + Hamming). JXJ 572 has sense-strand sites within 2 mismatches in 10,302 genes; the next most promiscuous probe, JXJ 403 (in an L1), has 265, and the other 78 have 4-45. Of the 79 batch-2 signature genes, 59 carry a JXJ 572 site within 3 mismatches (51 within 1) and the remaining 20 all carry a sense exonic Alu that matches at 4+. Across all expressed genes, exon log2({40old, 40even, 80} / {20old, 40odd}) by best exonic JXJ 572 site: 0 mismatches +0.77 (n=7), 1 mismatch +0.22 (n=267), 2 mismatches +0.03 (n=100), 3 mismatches -0.01, no site -0.04 (n=9,813); fraction of genes above 2-fold 29%, 14%, 9%, 3%, 0.5%. One oligo, one dose curve.")]));
P(b(0, [run("Things in the first pass this overturns. ", { bold: true }), run("The 40old/40even/80 vs neg contrast was never confounded with handling, it was confounded with an Alu probe, which is worse for mature-RNA features and irrelevant for intronic ones. TMEM9B-AS1 (20old +0.89, 40old +2.94, 40odd +0.63, 40even +2.80, 80 +2.48, 4 sense exonic Alus) and EMSY-DT (-0.24, +2.81, +0.28, +2.59, +2.36, 2 Alus), two of the three genes I proposed for the reporter, are JXJ 572 captures. The \"208 even-oligo genes carry no more NEAT1 16-mers than controls\" test missed all of this because a 20-mer with 1-3 mismatches need not share a 16-mer.")]));
P(img("fig6_offtarget_sites.png", 630, 220));
P(cap([run("Fig 6. The batch-2 signature is cross-hybridisation of sense Alus to probe JXJ 572. ", { bold: true }),
       run("(a) Mean exon log2({40old, 40even, 80} / {20old, 40odd}) by number of exonic Alu elements in the sense (black) or only the antisense (grey) orientation, with s.e.m.; true-split rank among the 210 ways of choosing 4 \"batch 1\" libraries. (b) Coverage along TUB in 40even, 40odd and neg; shaded bar is the JXJ 806 target site (2 mismatches). (c) The same for MAP3K13 and the JXJ 572 site (1 mismatch).")]));
P(img("fig5_probe_signatures.png", 630, 215));
P(cap([run("Fig 5. Probe-set signature genes. ", { bold: true }),
       run("(a) Exon log2(set / neg) for the 86 genes whose five-set profile fits one membership vector (contrast > 1 log2, R2 > 0.8), sorted by the vector. (b) The six probe classes and the sets each belongs to. (c) Five examples, one per class: PSRC1 and TUB are the new-design classes, DHRS4-AS1 the JXJ 572 class.")]));

// ============================ 13. nascent RNA directly ============================
P(p([run("Nascent RNA, tested directly: intron vs exon", { bold: true })]));
P(b(0, "Counted exon-union and intron reads per gene separately (featureCounts -s 1; intron = gene span minus the union of same-strand exons of any gene). 8,912 genes with enough reads of both kinds. This is the Wu et al. Group A/B test I asked for in the first pass."));
P(b(0, "Raw, it looks great: median exon log2(probe/neg) -0.15, intron +0.09, intron minus exon +0.23 and positive for 76% of genes; holds in every library; near NEAT1 it is +0.45. 66-split rank 1/66 genome-wide. I nearly stopped here."));
P(b(0, "Then the same GC worry as before, but sharper: introns are AT-rich and the negs are the two most GC-high libraries, so intron-over-exon is exactly what the GC artifact produces too. Two GC-free checks:"));
P(b(1, "Reads crossing an exon-intron boundary unspliced vs reads spliced across the same boundary (same loci, same sequence). Unspliced fraction 5.2-5.5% in both negs, 5.8-7.1% in 8 of 10 pulldown libraries, but 5.1-5.2% in both 80 libraries. 66-split rank 5/66, p = 0.08. Near NEAT1 vs elsewhere on this statistic: rank 36/66, nothing (but only 41 near genes have junction reads)."));
P(b(1, "Region-level CQN: for each library fit log2 CPM ~ spline(GC) + spline(log length) separately for the exon part and the intron part of each gene, keep residuals, response = mean residual over 10 probe minus mean over 2 neg. Residual GC association drops to about -0.01 for both region types."));
P(b(0, [run("After the correction the genome-wide intron excess is gone completely. ", { bold: true }), run("Median corrected intron minus exon -0.011, fraction positive 0.486, 66-split rank 56/66. So \"the pulldown is enriched for unspliced RNA genome-wide\" was the GC axis again, in a new costume. The unspliced-junction statistic (which has no GC confound) agrees that whatever is there is small and absent in the 80 set. The per-library intron fraction tracks sequencing depth (rho = +0.85 across the 10 pulldowns), not NEAT1 capture efficiency (+0.30), so what nascent component exists looks like prep-to-prep chromatin carry-over, not a property of genes.")]));
P(b(0, [run("What survives is the cis part, and it is intronic. ", { bold: true }), run("Corrected intron response for genes within 2 Mb of NEAT1 minus the rest: +0.267, rank 1/66 (null -0.19 to +0.27). Exon response for the same genes: +0.079, also rank 1/66 but a third the size. Restricted to Alu-free, non-signature genes: intron +0.241 rank 1/66, exon +0.067 rank 2/66. Protein-coding only: intron +0.233, exon +0.053. 39 of the 59 genes have a positive intron response; leave-one-out range +0.233 to +0.277.")]));
P(b(1, "It is in every probe set, including the two that have no nascent signal genome-wide: 20old +0.20, 40old +0.34, 40odd +0.16, 40even +0.38, 80 +0.26 (near-minus-rest, intron); every single library is +0.14 to +0.42. The exon version is +0.04 to +0.10 in each set."));
P(b(1, "Shape: intron +0.50 within 0.5 Mb (n=17), +0.22 at 0.5-1 Mb, +0.15 at 1-2 Mb, then -0.03, -0.13, +0.07, +0.07 and -0.01 beyond 20 Mb. Exon: +0.18 within 0.5 Mb and nothing after 1 Mb. Sliding a 2 Mb window across all chromosomes (74 windows with 40-160 genes): chr11 65.2, 66.2, 64.2 and 67.2 Mb are the top four, i.e. the NEAT1 window and its three neighbours, then chr17 5.6 Mb at +0.14."));
P(b(1, "Top genes by corrected intron response within 2 Mb: VEGFB +1.55 (exon +0.12), SCYL1 +1.39 (exon -0.15), VPS51 +1.06, DPF2 +0.90, FOSL1 +0.85, EIF1AD +0.84 (exon -0.32), RAD9A, MAP4K2, PLCB3, MUS81, BRMS1. These are ordinary housekeeping genes whose introns come down with NEAT1 and whose mature mRNAs do not."));
P(b(1, "Not NEAT1 read-through. NEAT1 transcription does run past the annotated 3' end: + strand reads in the first 20 kb downstream are enriched log2 +2.3 to +2.7, then nothing (+0.05 at 30-40 kb, -0.2 to +0.1 from 60 kb on, MALAT1 at +50-60 kb is -0.19). A read-through artifact would only touch + strand genes downstream. The cis intron effect is upstream +0.23 (n=24) and downstream +0.29 (n=35), + strand +0.28 (n=38) and - strand +0.23 (n=21), and +0.25 over the 37 genes that are not + strand downstream. VEGFB, the top gene, is 1.2 Mb upstream."));
P(b(0, "Reproducibility of the corrected responses between the two independent library halves: exon r = 0.67, intron r = 0.50, intron-minus-exon r = 0.56."));
P(img("fig7_nascent_cqn.png", 630, 220));
P(cap([run("Fig 7. After GC correction the only nascent signal is cis to the NEAT1 locus. ", { bold: true }),
       run("(a) Corrected exon response against corrected intron response per gene, n = 8,912. (b) Mean corrected response for genes within r Mb of NEAT1, intron (black) and exon (grey), 95% bootstrap intervals; dashed lines are the genome means. (c) 66-split permutation for the contrasts named; black is the true split, grey the other 65, the rank is printed.")]));

// ============================ 14. lncRNA, Yap ============================
P(p([run("The lncRNA effect, re-examined", { bold: true })]));
P(b(0, "It is an exon effect, not an intron effect: lncRNA minus other, corrected exon +0.115 rank 3/66; intron +0.000 rank 35/66. So lncRNAs are not coming down as nascent transcripts."));
P(b(0, "lncRNA exons are Alu-rich (57% carry a sense exonic Alu vs 43% of mRNAs), so part of the first-pass lncRNA effect was JXJ 572. Removing every gene with a sense exonic Alu and every signature gene: lncRNA minus protein-coding, corrected exon, +0.068, rank 9/66, p = 0.14. Raw (uncorrected) the same contrast is +0.154, rank 3/66, and it is the same size in the sets with JXJ 572 (+0.15) and without (+0.16), so it is not the Alu probe, but after both corrections it is no longer a clean result."));
P(b(1, "By name class, Alu-free genes, raw exon: -DT +0.14 (n=7), unnamed lncRNA +0.05 (81), -AS +0.04 (29), LINC -0.17 (18), SNHG -0.03, protein-coding -0.18. The divergent/antisense/unannotated class is still the one above the mRNAs, and the intergenic LINCs are below them. I keep this as suggestive and nothing more."));
P(b(0, [run("Yap, Chung and Makeyev 2022 (HyPro-seq on NEAT1, same HeLa). ", { bold: true }), run("Their list of 267 NEAT1 targets (DESeq2 vs an enzyme-only, no-probe control) reads like a description of the cis effect above: 121 of 267 are on chr11, \"significant clusters of NEAT1 target genes mapped to a ~16 Mb region on the q-arm of chr11 comprising the NEAT1 locus itself\"; targets have higher intron retention than Pol II-matched controls; they interpret the transcripts as \"likely newly synthesized and retained in the vicinity of their transcription sites\". Gene-by-gene concordance here is weak: 65 of their targets pass my filters, corrected intron response median +0.05 vs -0.02 (MWU p = 0.27), exon +0.02 vs 0.00 (p = 0.10), and the 37 off-chr11 targets show nothing. So the two methods agree on the neighbourhood and not on individual genes, which is what you expect if both are picking up a small cis effect with different noise.")]));
P(b(0, "Hi-C. Only a GM12878 proxy was reachable (4DN HeLa needs a login). Coverage-normalised observed/expected contact between each chr11 gene bin and the NEAT1 bin at 250 kb adds nothing beyond linear distance: in a joint model of log distance and log O/E the O/E term has p = 0.4-0.7 for intron, exon and intron-minus-exon. Wrong cell line, so this only says the proxy does not help; it needs HeLa Hi-C or a 4C from the NEAT1 locus."));

// ============================ 15. where I land now ============================
P(p([run("Where I land now", { bold: true })]));
P(b(0, [run("Hypothesis (revised): nothing in this dataset is enriched because it is in a paraspeckle. The pulldown captures NEAT1_2 fragments at the 80 probe sites; everything else that comes down is either (i) mature RNA annealing directly to one oligo, dominated by sense-Alu transcripts on JXJ 572, or (ii) a small, GC-robust, intronic enrichment of genes transcribed within ~1-2 Mb of the NEAT1 locus. There is no sequence-encoded targeting of mature mRNAs, no genome-wide nascent-RNA enrichment, and the GC/length/abundance axes are technical.", { bold: true })]));
P(b(1, "The strongest statement I can make, and I have checked it every way I can: no transcript other than NEAT1 is enriched more than 2-fold at the exon level in all five probe sets. Three genes clear 1.7-fold (an unnamed chr22 lncRNA, CMTM6, HIC2), none near NEAT1, none with a story. If paraspeckle residents like CTN-RNA, NICN1 or the IRAlu mRNAs were being co-captured with the NEAT1 body, they would be in every set; they are in none. The capture is at the fragment level (Fig 4), which is why: a fragment of a resident RNA has no probe site, and the crosslink to NEAT1 that would carry it along is evidently not what this prep preserves."));
P(b(1, "(i) is now proven rather than inferred: the membership vectors (Fig 5), the pile-ups on the exact target sites (Fig 6b,c), the Alu dose curve (Fig 6a), the genome-wide site scan, and the collapse of the \"batch\" when Alu genes are removed all point at the same oligo. It accounts for the batch split, the 86 signature genes, TMEM9B-AS1 and EMSY-DT, and most of the first-pass lncRNA effect."));
P(b(1, "(ii) is the one piece of biology. It is rank 1/66 before and after GC correction, in Alu-free genes, in protein-coding genes only, in each of the five sets and in each of the ten libraries, it decays over ~1-2 Mb with the right shape, it is the top 2 Mb window in the genome together with its three chr11 neighbours, and it is intronic, with mature mRNAs of the same genes barely moving (+0.08). That is nascent transcripts of NEAT1's neighbours coming down with NEAT1. It is also what Yap et al. saw with an orthogonal method in the same cells (chr11q clustering, intron retention, \"retained in the vicinity of their transcription sites\")."));
P(b(1, "Mechanistically (ii) does not need a paraspeckle either. NEAT1 is transcribed at chr11:65.42 Mb and paraspeckles sit on the locus; nascent transcripts of genes in the same ~1-2 Mb neighbourhood are physically there, and a NEAT1 fragment crosslinked to chromatin or to a neighbouring nascent transcript carries that transcript's intron down with it. The signal is a map of what is near the NEAT1 gene, not of what is in the paraspeckle body."));
P(b(0, [run("What I am not claiming. ", { bold: true }), run("The lncRNA exon effect is +0.07 and p = 0.14 once Alu genes are out and GC is corrected; I carry it as a lead, not a result. The genome-wide unspliced-RNA excess is small (p = 0.08), absent in the 80 set and tracks depth, so I read it as chromatin carry-over. And the cis effect is +0.27 log2 over 59 genes: real, reproducible, and small, which is the honest size of \"paraspeckle localisation\" signal in a hybridisation-capture experiment with a bead-only control.")]));

// ============================ 16. next ============================
P(p([run("What I want to do next", { bold: true })]));
P(b(0, "Drop JXJ 572 from the probe set, or at minimum analyse 20old and 40odd (which lack it) as the clean sets. Also drop JXJ 403 and JXJ 406 (the other two with hundreds of genomic sites). Any future NEAT1 probe design should be screened against RepeatMasker; a probe in an Alu is a transcriptome-wide capture reagent."));
P(b(0, "The control that would settle (i) vs (ii) for good is a scrambled / non-targeting probe set with the same handle, run through the same prep. Beads-only does not control for annealing, and I have now shown that annealing is the largest term."));
P(b(0, "Test (ii) as a locus effect rather than a paraspeckle effect: a NEAT1-knockdown or paraspeckle-disrupting condition (NONO/SFPQ knockdown, or the NEAT1_2-deficient line) should leave the cis intronic signal in place if it is nascent-transcript proximity to the NEAT1 gene, and remove it if it needs the body. Alternatively, pull down a different locus-tethered nascent RNA on chr11 and ask whether NEAT1's neighbours come along."));
P(b(0, "4C or HeLa Hi-C from the NEAT1 promoter, to see whether the 1-2 Mb intronic halo is the contact domain. The GM12878 proxy cannot answer this."));
P(b(0, "More neg replicates (4-6), and a nuclear-extract input lane. Two negs is still why GC remains a one-sided argument everywhere above."));
P(b(0, "If anyone wants to look for genuine paraspeckle residents in this design: the test is \"enriched in all five sets, exon-level, no probe site\", and nothing passes it at 2-fold. A denser probe set with no repeat-overlapping oligos and a formaldehyde-only (no fragmentation) variant would be the way to recover co-captured residents, since whole-body capture needs intact crosslinked complexes."));

// ============================ 17. files (second pass) ============================
P(p([run("Files (second pass)", { bold: true })]));
P(b(0, "work/probes80.tsv - the 80 probes with region, parity, series, target 20-mer, position on NEAT1, set membership; work/probes80_repeats.tsv - the 19 probes that overlap a RepeatMasker element"));
P(b(0, "work/counts_exon.tsv, counts_intron.tsv, counts_genespan.tsv, counts_junc_unspliced.tsv, counts_split.tsv - featureCounts -s 1 per library; counts_genespan reproduces NEAT1_input.xlsx exactly"));
P(b(0, "work/neat1_cov_*.tsv - 100-bp coverage along NEAT1 per library; analysis/probe_scan_genome.py - the genome-wide site scan (work/probe_hits.tsv, probe_best_per_gene.pkl)"));
P(b(0, "work/alu_counts.pkl - sense/antisense Alu counts per gene, exonic and intronic; work/probe_pattern_exon.pkl - the membership-vector fits"));
P(b(0, "analysis/nascent_cqn.py - region-level GC/length correction of exon and intron counts, the 66-split permutation on it, cis decay, window null, HyPro concordance (work/nascent_cqn.pkl, nascent_perm66.pkl); analysis/fig7_nascent.py"));
P(b(0, "work/hic_chr11_gene_feature.pkl - GM12878 Hi-C O/E contact with the NEAT1 bin per chr11 gene; work/uniform_exon_enriched.tsv - the genes above 1.7-fold in all five sets"));
P(b(0, "ext/yap2022/neat1_hypro_targets.tsv - the 267 HyPro-seq targets from Yap et al. 2022 Table S2"));

const doc = new Document({
  numbering,
  styles: { default: { document: { run: { font: FONT, size: SZ } } } },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } },
    children: kids,
  }],
});
Packer.toBuffer(doc).then(buf => { fs.writeFileSync("Regev Lab notes - paraspeckle determinants.docx", buf); console.log("written"); });
