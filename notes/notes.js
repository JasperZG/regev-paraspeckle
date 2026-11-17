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

const doc = new Document({
  numbering,
  styles: { default: { document: { run: { font: FONT, size: SZ } } } },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } },
    children: kids,
  }],
});
Packer.toBuffer(doc).then(buf => { fs.writeFileSync("Regev Lab notes - paraspeckle determinants.docx", buf); console.log("written"); });
