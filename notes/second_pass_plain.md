Paraspeckle localization, second pass

Question stays the same: which property of an RNA changes its tendency to associate with paraspeckles, measured as recovery in the NEAT1-targeted pulldown relative to NEAT1_neg. The first pass ended with GC, length and Alu as the leading features and a working idea that intronic RNA is favoured. Before building on that, go back to the raw material and check what the assay is actually measuring.

Going back to the raw data

The shared folder has more than the DE tables. It has the 12 deduplicated BAM files, the raw count table, a note describing the five designs, and the probe workbook with the 80 sequences.

A BAM file records where each sequencing read landed on the genome. Each library here is single-end, 100 letters per read, deduplicated by molecular barcode, and stranded, meaning the read's orientation tells us which DNA strand the RNA came from. GAPDH gives 8,794 reads on its own strand and 0 on the opposite strand, so strand information is reliable.

The cell line was still unlabelled. The raw reads contain the HPV18 E6/E7 viral genes (1,066 reads per 4 million in NEAT1_neg_Rep1, none on the opposite strand), there are no chrY reads, and CDKN2A is high. HPV18 integration is the signature of HeLa. This agrees with the expression-correlation inference from before and settles it.

The supplied count table turns out to be whole-gene-span counts. Recounting the BAMs over each gene's full annotated span, strand-aware, reproduces NEAT1_input.xlsx exactly, all 67,928 genes by 12 libraries. Counting only exons gives 5.2 million reads for NEAT1_neg_Rep1 instead of 7.0 million. So every log2FoldChange in the DE tables is a mix of intronic and exonic reads, with 27 to 33% of reads intronic. This is why intron length and the intronic read fraction looked important before: they describe what the count is made of, not where the RNA goes.

The probes

A probe is a short DNA oligo with a 20-letter segment complementary to one spot on NEAT1 plus a 35-letter shared tail. The workbook has 80 of them in ten regions of eight. Mapping each 20-letter segment back onto the NEAT1 sequence, every probe matches exactly once, in the sense orientation, between positions 3,846 and 21,577. NEAT1 has a short form, NEAT1_1, which is the first 3.7 kb, and a long form, NEAT1_2, which is the full 22.7 kb. No probe touches the first 3.7 kb, so only NEAT1_2 is targeted.

The experiment note defines the designs by region, not by alternating oligos. 40odd is the 40 probes in regions 1, 3, 5, 7, 9; 40even the other 40; 80 is both. 20old is the 16 oldest probes (JXJ 401 to 420); 40old is those plus 24 more (JXJ 535 to 576); JXJ 783 to 822 are the newest 40 and appear only in 40odd, 40even and 80.

This means each probe is present in a specific combination of designs. There are six possible combinations. A gene that is pulled down because its RNA sticks to one probe should be enriched in exactly the designs that contain that probe. A gene that is pulled down because it is near NEAT1 should be enriched in all five. That difference is testable.

What the assay actually captures

Look at read coverage along NEAT1 itself, in 100-letter bins, normalised to reads per million.

In NEAT1_neg, the first 3.7 kb has about 105 reads per million per letter and the rest about 15. That is the known picture: HeLa makes roughly six times more NEAT1_1 than NEAT1_2.

In any targeted design, coverage is flat at 300 to 370 across the whole 22.7 kb. The NEAT1_2 body is enriched 25 to 90-fold. The first 3.7 kb is enriched only 2.5-fold, which is what you get when that region comes along only as the front end of NEAT1_2 molecules. The abundant NEAT1_1 molecules, which have no probe site, are not recovered.

If the pulldown captured whole paraspeckles, or intact NEAT1 with everything attached, NEAT1_1 would come along, because it is in the same bodies. It does not.

Now compare 40odd with 40even along NEAT1. The ratio alternates with the regions: +0.45, -0.46, +0.77, -0.64, +1.10, -1.22, +0.29, -1.42, +1.50, -1.41 log2 for regions 1 to 10. Each design recovers two to three times more of the stretch its own probes sit on. 20old does the same, peaking over the two stretches that hold its 16 probes.

This is only possible if the RNA has been broken into pieces shorter than a region. A simple simulation makes it concrete. Cut the RNA into fragments of length L, keep the fragments that contain a probe site, and predict the odd/even profile. Fragments of 200 to 400 letters reproduce the observed profile with r = 0.77; 1 kb gives 0.63; 2 kb gives 0.10.

So the assay is hybridisation capture of fragmented RNA. A fragment is recovered if it carries a probe site. A fragment from any other RNA can be recovered in three ways: it anneals to a probe itself, it is crosslinked to a NEAT1 fragment, or it sticks to the beads. NEAT1_neg controls only for the beads.

Fig 4. The pulldown captures NEAT1_2 fragments at the probe sites. (a) Coverage along NEAT1 in neg (grey) and 80 (black); dotted line is the end of NEAT1_1. (b) log2(40odd / 40even) along NEAT1 with odd-region probes ticked above and even-region probes below. (c) log2(design / neg) for 20old, 40old and 80; ticks are the 16 probes that make up 20old.

Testing direct probe binding with the design combinations

Take each gene's five exon-level enrichment values (one per design) and ask whether they match one of the six probe combinations: high in exactly the designs containing one probe class, flat in the others. Keep genes where the fit is clean (contrast above 1 log2, R-squared above 0.8).

86 genes pass. Shuffling the design labels within each gene gives 0.7 genes per shuffle, so these are not chance.

79 of the 86 follow the combination {40old, 40even, 80}. That is also exactly the "batch" split seen in the first pass, so on its own it is ambiguous.

The other 7 are not ambiguous. PSRC1, MCU, KCNT2 and PIEZO2 are enriched only in 40odd and 80 (PSRC1 is +3.3 and +2.9 there, +0.1 to +0.4 elsewhere). TUB is enriched only in 40even and 80 (+3.2, +2.8, flat elsewhere). No handling difference can make a gene go up in 40odd and 80 but not in 40old, 40even or 20old. Only a probe present in exactly those designs can.

These are exon effects. For the 86 genes, exonic enrichment in the designs that contain the matching probe is +1.57 log2 in the median and -0.07 in the others; intronic enrichment is +0.40 versus +0.06. Mature, spliced RNA is being captured, which is what annealing to an oligo looks like.

Finding the binding sites

Pull up read coverage along each of these genes in the matching design and the non-matching design, and look for where the reads pile up.

TUB has 558 reads in 40even and 6 in 40odd over one 200-letter bin. That bin contains a sequence matching the JXJ 806 target site with 2 mismatches in 20 letters. JXJ 806 is a new even-region probe, present in exactly 40even and 80. PSRC1 piles up on a JXJ 818 site (4 mismatches), MCU on a JXJ 815 site (3 mismatches).

Every batch-2 gene checked piles up on the same probe. MAP3K13 on a JXJ 572 site with 1 mismatch, 114 reads in 40even versus 1 in 40odd. SIRT7, JXJ 572, 1 mismatch, 158 versus 3. DHRS4-AS1, JXJ 572, 2 mismatches, 135 versus 2.

JXJ 572 sits in an Alu

RepeatMasker marks an AluSz6 element inside NEAT1 at position 16,760, in the same orientation as NEAT1. The JXJ 572 target is inside it. Its 20 letters, AAATTAGTCAGGCATGGTGG, are 3 mismatches from the Alu consensus AAATTAGCCGGGCGTGGTGG; the reverse complement is 10 mismatches away. No other probe is within 6 of the consensus.

Alu is the most common repeat in the human genome and a large fraction of transcripts carry one. A probe matching Alu will anneal to any RNA fragment carrying an Alu in the sense orientation, with a strength that depends on how close that copy is to the probe.

Check this against the 79 batch-2 genes: 78 of them carry a sense-oriented Alu in their exons. The background rate among all tested genes is 43%. The median count is 3 sense exonic Alus versus 0.

Check it genome-wide as a dose curve. Mean exonic log2({40old, 40even, 80} / {20old, 40odd}) by number of sense exonic Alus: 0 Alus -0.044 (5,838 genes), 1 Alu -0.006, 2 to 3 Alus +0.060, 4 or more +0.125. Antisense-only Alus are flat or negative at every count. Of the 210 ways to choose which 4 targeted libraries count as "batch 1", the true choice gives the largest sense-Alu contrast.

Then check whether the batch exists at all without Alu genes. Among all genes, the libraries cluster into {20old, 40odd} and {40old, 40even, 80}. Among genes with no sense exonic Alu, the clustering collapses to replicate pairs and the two groups are indistinguishable. Among genes with 2 or more sense exonic Alus, the split is stronger than before. The batch is sets-with-JXJ-572 versus sets-without.

Finally, scan the whole genome for every probe's target site at up to 3 mismatches. JXJ 572 has sense-strand sites within 2 mismatches in 10,302 genes. The next most promiscuous probe, JXJ 403, has 265; the rest have 4 to 45. Across all expressed genes, enrichment in the JXJ 572 designs by closeness of the best exonic site: perfect match +0.77 log2, 1 mismatch +0.22, 2 mismatches +0.03, 3 mismatches -0.01, no site -0.04. One oligo, one dose curve.

Fig 6. The batch-2 signature is cross-hybridisation of sense Alus to probe JXJ 572. (a) Mean exon log2({40old, 40even, 80} / {20old, 40odd}) by number of exonic Alus in the sense orientation (black) or only antisense (grey). (b) Coverage along TUB in 40even, 40odd and neg; shaded bar is the JXJ 806 site. (c) The same for MAP3K13 and the JXJ 572 site.

What this changes from the first pass

The Alu result was the probe. The first-pass check for probe binding searched for exact 15- and 20-letter matches to the probe segments, and a 20-letter site with 1 to 3 mismatches shares no 16-letter exact match, so it was missed. The strongest first-pass genes, TMEM9B-AS1 and EMSY-DT, are JXJ 572 captures: enriched +2.5 to +2.9 in the three designs with the probe and +0.9 or less in the two without.

The batch split was never handling. It was one oligo.

Is anything enriched in every design?

NEAT1 itself is +3.3 to +3.7 log2 at the exon level in all five designs. Ask which other genes are enriched in all five.

Exon-level enrichment above 1 log2 (2-fold) in all five designs: 0 genes. Above 0.75: 3 genes (an unnamed lncRNA on chr22, CMTM6, HIC2), none near NEAT1, none with any known paraspeckle connection. Above 0.5: 19.

Known paraspeckle-retained RNAs, such as the inverted-Alu messenger RNAs, are not among them. They would have to be crosslinked to a NEAT1 fragment to come along, and whatever this prep preserves, it does not preserve that.

Nascent RNA, tested directly

Count exonic and intronic reads per gene separately, as before, now for all 12 libraries. Raw, the pulldown looks enriched for intronic RNA: intronic minus exonic enrichment is +0.23 log2 in the median, positive for 76% of genes, and +0.45 for genes within 2 Mb of NEAT1.

The GC problem returns in a sharper form. Introns are AT-rich and the two NEAT1_neg libraries are the two most GC-high libraries of the twelve, so intron-over-exon is exactly what that drift would produce. Two GC-independent checks.

First, reads that cross an exon-intron boundary unspliced versus reads spliced across the same boundary. Same loci, same sequence, so GC cancels. Unspliced fraction is 5.2 to 5.5% in both negs, 5.8 to 7.1% in eight targeted libraries, but 5.1 to 5.2% in both 80 libraries. Against the 66 possible choices of control pair this ranks 5 of 66, p = 0.08, and the 80 design shows nothing.

Second, redo the GC correction at the region level: fit each library's exon counts and intron counts separately against GC and length, keep the residuals, and recompute enrichment. After this the genome-wide intronic excess is gone: median intron-minus-exon -0.011, 48.6% positive, rank 56 of 66. The intron-fraction result from the first pass was the GC axis as well.

Per library, how much intronic RNA comes through tracks sequencing depth (rho = +0.85 across the ten targeted libraries), not how well NEAT1 was captured (+0.30). Whatever nascent component exists looks like chromatin carry-over that varies prep to prep.

What survives is local to the NEAT1 gene

After the correction, genes within 2 Mb of the NEAT1 locus on chromosome 11 have intronic enrichment +0.267 log2 above the rest of the genome, rank 1 of 66. Their exonic enrichment is +0.079, a third of the size. Restricted to genes with no sense Alu, +0.241, still rank 1. Protein-coding only, +0.233. Present in every design (20old +0.20, 40old +0.34, 40odd +0.16, 40even +0.38, 80 +0.26) and in every one of the ten libraries (+0.14 to +0.42).

It decays with distance: +0.50 within 0.5 Mb, +0.22 at 0.5 to 1 Mb, +0.15 at 1 to 2 Mb, then nothing. Sliding a 2 Mb window across the genome, the NEAT1 window and its three chromosome 11 neighbours are the top four of 74.

It is symmetric: +0.23 upstream of NEAT1 and +0.29 downstream, +0.28 on the plus strand and +0.23 on the minus strand. NEAT1 transcription does run 20 kb past its annotated end, but a read-through artifact would only affect plus-strand genes downstream, and the top gene, VEGFB, is 1.2 Mb upstream.

The top genes are ordinary: VEGFB intron +1.55 with exon +0.12, SCYL1 +1.39 with exon -0.15, VPS51, DPF2, FOSL1, EIF1AD. Their unspliced transcripts come down with NEAT1 and their mature messenger RNAs do not.

Yap, Chung and Makeyev (2022) used a proximity-labelling method on NEAT1 in HeLa and reported that their targets cluster in a 16 Mb region of chromosome 11 around NEAT1, show more intron retention, and are "likely newly synthesized and retained in the vicinity of their transcription sites". Gene-by-gene agreement with our list is weak (p = 0.27 for introns), but both methods see the same neighbourhood.

Fig 7. After GC correction the only nascent signal is cis to the NEAT1 locus. (a) Corrected exon response against corrected intron response per gene. (b) Mean corrected response for genes within r Mb of NEAT1, intron (black) and exon (grey). (c) The 66-split permutation for each contrast; black is the true split, grey the other 65.

The lncRNA class

Long noncoding RNA exons are Alu-rich (57% carry a sense exonic Alu versus 43% of messenger RNAs), so part of the first-pass lncRNA effect was JXJ 572. Removing every sense-Alu gene and correcting GC, lncRNA minus protein-coding is +0.068 log2, rank 9 of 66, p = 0.14. Kept as a lead, not a result.

Where this leaves the question

Put the pieces back into the log2FoldChange column. Predicting it out of sample from gene properties, per design: GC, length and abundance explain 4 to 22%, and the permutation test says that is library drift. The Alu match to JXJ 572 explains about 3% in the three designs that contain it and 0% in the two that do not. The distance to NEAT1 rounds to 0% genome-wide because it is 59 genes out of 12,000, in introns. Everything else, 75 to 93% of the column, is not predictable from any gene feature we have.

So the answer to "what property of the gene determines this value" for this dataset is: mostly nothing about the gene. The pulldown recovers NEAT1_2 fragments at the probe sites. Beyond that it recovers mature RNAs that anneal to one Alu-overlapping probe, and a small amount of unspliced RNA from genes transcribed within 1 to 2 Mb of the NEAT1 gene. No mature RNA other than NEAT1 is reliably recovered by every design, and no sequence feature of mature RNA predicts recovery.

This is a statement about the assay as much as about biology. Fragment-level capture with a beads-only control cannot see what is co-resident in a paraspeckle unless the crosslink to NEAT1 survives, and here it does not. To ask the original question, the design has to change: remove JXJ 572 and any other repeat-overlapping probe, add a scrambled-probe control with the same tail so that annealing on its own can be subtracted, and run more than two negative libraries. A NEAT1 or NONO/SFPQ knockdown would then tell us whether the chromosome 11 signal needs the paraspeckle body or only the NEAT1 gene.

---

Third pass, 19 November: the rest of the data arrived, including a batch without the bad probe

The full rtRT-seq folder came through: NGS1 through NGS6, the MALAT1 runs, the April NEAT1/NEAT1_SA sets, FISH images, the probe design files, and Jingxuan's own analysis write-ups. Two things in it change the picture.

They found JXJ 572 independently, and before I did. Their specificity check from 29 October names the same probe, notes that its target is a common Alu variant, finds a near-match in 72% of the enriched genes against 17% of length-matched controls, and shows reads stacking on the match site in P3H4, MTRFR, BRI3BP, VPS51 and GGT4P at 7 to 34 times the sets without the probe. They then left it out of the next design. Two analyses, opposite directions, same answer: they went from the gene list to the probe, I went from the probe-set membership patterns to the probe. The one place they stopped short is the conclusion that direct hybridisation explains at most 5 of the 60 genes, which comes from applying a 16-continuous-base-pair criterion meant for a free oligo in solution. The mismatch dose curve and their own pile-ups say the weaker sites work too.

I also had the method wrong. I described this as hybridisation capture of fragmented RNA. It is rtRT-seq: a reverse transcriptase tethered to each probe labels RNA within its reach. Their own measurement puts the reach at about 64-fold enrichment within 100 nt of a probe, 27-fold at 200 to 400 nt, 17-fold at 400 to 800 nt, settling at a 6 to 10-fold floor, with RT starts 7.6 times higher just upstream of a bound probe at 19 of 20 sites. My 200 to 400 nt fragment estimate was measuring that reach, not a fragment size. The conclusion does not change and is now better supported: the labelling is local along the RNA, so this cannot capture a whole paraspeckle.

NGS6 tests everything I wrote down

NGS6 has no JXJ 572, a new 3' probe design that I had already screened and found clean, five positive sets, and six negative samples instead of one. Every second-pass prediction is testable against it, and I ran them before reading their conclusions.

The probe-signature genes go flat. The 86 genes are at median +0.08 to +0.26 log2 depending on the set, which is the global shift. P3H4 is the most depleted at -0.25 to -0.56. MTRFR, SIRT7, DHRS4-AS1, MAP3K13, TMEM9B-AS1 and H3-3A are all flat. The sense-Alu dose curve is gone: no Alu +0.12, four or more +0.18, rho +0.02.

TUB is the clean prospective confirmation. In NGS4 it was enriched in exactly 40even and 80, which I traced to a 2-mismatch JXJ 806 site. JXJ 806 is a Region 6 probe. NGS6's three 3'-end sets have no Region 6 probes, but MED covers Regions 2 to 6 and TOTAL is a whole-NEAT1 design. Prediction: up in MED and TOTAL only. Result: 3' ODD +0.42, 3' EVEN +0.21, 3' total +0.28, MED +2.78, TOTAL +2.35. The probe explanation predicted the sign and the size of an effect in data I had not seen.

GC is settled. The first pass could not separate it from library drift with two negatives and I left it open. With twelve negative libraries and the same spline correction, the GC association of the response is -0.004 against a null spread of 0.007, p = 0.73 over 2,000 random splits. It was technical. The extra negatives were the whole answer, and that is worth remembering for every future design.

The one positive result I had does not replicate

The cis intronic halo was the thing I kept. NGS6 does not support it. The near-NEAT1 contrast is +0.053 and formally beats a random-split null, but the shape is wrong: within 0.5 Mb it is -0.048, within 1 Mb -0.026, rising to +0.07 out at 10 Mb. In NGS4 it was +0.50 within 0.5 Mb and gone by 3 Mb. A mechanism anchored at the NEAT1 locus does not give you a signal that is absent at the locus and present 10 Mb away. Sliding a 2 Mb window across the genome, the NEAT1 window ranks 44th of 240 in NGS6, behind the chromosome 6 histone cluster; in NGS4 it and its three neighbours were the top four of 74. And the genes do not match: the 15 strongest NGS4 cis genes have a median NGS6 response of -0.087 against a genome median of +0.003. Gene-level agreement between the two batches is r = +0.08 genome-wide and +0.12 on chromosome 11.

The fair caveat is that only the gene-span count table came with NGS6, not its BAMs, so I cannot do the intron/exon split there, and that is exactly where the NGS4 signal lived. Gene-span counts dilute it about threefold. On top of that, their own QC calls NGS6 a failed batch for co-capture: the positives pull NEAT1 down 7 to 16-fold but bring nothing else along, and the negatives are depleted for 160 to 400 genes. So this is a weak test rather than a refutation. But it is the only independent test available and it does not reproduce the effect, so the cis claim moves from "the one thing that survives" to "unconfirmed, needs the NGS6 BAMs".

What is actually enriched, across both batches

The strictest test the data supports: up in all five NGS4 sets and all five NGS6 sets, against pooled negatives and against each set's own matched negative. Ten positive libraries, two preps, two probe designs. Above 2-fold: one gene, NEAT1. Above 1.7-fold: NEAT1 and CMTM6. Above 1.4-fold: NEAT1, CMTM6, RUVBL1, HIC2. Above 1.2-fold: 60.

CMTM6 and RUVBL1 are the only leads worth a FISH probe. Both are up in all ten NGS6 positive libraries and all five NGS4 sets, both are exonic effects, and neither has a probe site within three mismatches. CMTM6 runs 116 to 205 CPM in the NGS6 positives against 11 to 46 in every negative; RUVBL1 219 to 386 against 38 to 58. Against them: the effect size does not match between batches, roughly 1.8-fold in NGS4 against 3 to 5-fold in NGS6, and part of the NGS6 gap comes from the negatives being depleted for these same genes, which is the known failure mode of that batch. So they are candidates to image, not results.

Their own number says why the list is this short. Working from the data, an RNA away from a probe gets labelled at most 15 to 25% as well as the probe site, and about 91% of each library is background, so reaching 2-fold needs roughly 30 to 45% of that RNA's molecules sitting in paraspeckles. Paraspeckles are well under 1% of the nucleus, so even an RNA concentrated 3 to 4-fold in them, which is what P3H4 and MTRFR look like by FISH, should be invisible. Their calculation and my nearly empty cross-batch list are the same statement from two directions.

Where this leaves the question

The assay measures proximity to the probe along the RNA, and at the current background it cannot see paraspeckle residents at all. Of the log2FoldChange column, the explainable part is probe cross-hybridisation, dominated by JXJ 572; the rest is library and background variation. No property of a gene's mature RNA predicts recovery in either batch.

Three things are now settled rather than argued. The GC and length axis is technical, proven with twelve negatives. The batch axis was one Alu-overlapping oligo, proven by the dose curve, the pile-ups, and the effect vanishing when the probe is removed. And no mature RNA other than NEAT1 is reproducibly captured, one gene above 2-fold across ten positive libraries and two designs.

Two are open. Whether the cis intronic halo is real, which needs the NGS6 BAMs or a working batch with intron-level counts. And whether CMTM6 and RUVBL1 are genuine, which needs imaging.

One is now quantitative, and it is the important one. The sensitivity gap is 30 to 45% of molecules, and that is a design problem rather than an analysis problem. It explains every negative result in these notes. Cutting background buys more than any amount of reanalysis.
