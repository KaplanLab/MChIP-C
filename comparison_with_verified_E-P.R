library(dplyr)
library(tidyr)
library(GenomicRanges)
library(ggplot2)
#creating a proper list of tested E-P pairs from Gasperini et al., 2019
gasperini_neg <- read.csv("Documents/Jul2022_MChIPC/fig.4/GSE120861_all_deg_results.at_scale.txt.gz", header = T, sep = "\t")
genes <- unique(gasperini_neg[,c(9,16)])
genes <- filter(genes, !(duplicated(genes$gene_short_name)))
gasperini_neg <- filter(gasperini_neg, (quality_rank_grna == "top_two") & !(outlier_gene) & (fold_change.transcript_remaining <1.1) & (pvalue.empirical.adjusted > 0.1 | pvalue.empirical.adjusted == "not_applicable"))
gasperini_neg <- gasperini_neg[,c(12:14,9,7,16)]
gasperini_neg$Significant <- FALSE
gasperini_neg[,2] <- as.integer(gasperini_neg[,2])
gasperini_neg[,3] <- as.integer(gasperini_neg[,3])
gasperini_pos <- read.csv("Documents/Jul2022_MChIPC/fig.4/gasperini_positive.txt", header = T, sep = "\t")
gasperini_pos <- filter(gasperini_pos, high_confidence_subset)
gasperini_pos <- gasperini_pos[,c(9:11,3,7)]
colnames(gasperini_pos) <- colnames(gasperini_neg)[1:5]
gasperini_pos <- left_join(gasperini_pos,genes)
gasperini_pos$Significant <- TRUE
gasperini_pos[,5] <- as.character(gasperini_pos[,5])
gasperini <- bind_rows(gasperini_pos, gasperini_neg)
write.table(gasperini, "Documents/Jul2022_MChIPC/fig.4/Gasperini_EP_pairs.txt", col.names = T, row.names = F, sep = "\t", quote = F)
rm(list=ls())

#loading functionally tested E-P pairs from Gasperini et al., 2019 and filtering out silencers, long-distance and out of proper baits
gasperini <- read.csv("Documents/Jul2022_MChIPC/fig.4/Gasperini_EP_pairs.txt", header = T, sep = "\t")
gasperini$distance <- abs((gasperini$target_site.start+gasperini$target_site.stop)/2 - gasperini$target_gene.start)
gasperini <- filter(gasperini, Significant == TRUE & distance < 1000000 & distance > 5000)
baits <- read.csv("Documents/Jul2022_MChIPC/proper_baits.txt", header = T, sep = "\t")
gasperini_ranges <- GRanges(seqnames = gasperini$target_site.chr, ranges=IRanges(start=gasperini$target_gene.start, end=gasperini$target_gene.start+1, enh_id = seq(1,nrow(gasperini))))
bait_ranges <- GRanges(seqnames = baits$bait_chr, ranges=IRanges(start=baits$bait_start, end=baits$bait_end), enh_id = seq(1,nrow(baits)))
gasperini <- gasperini[as.data.frame(findOverlaps(gasperini_ranges, bait_ranges))$queryHits,]
#overlapping with MChIPC loops
loops <- read.csv("Documents/Jul2022_MChIPC/interaction calling/loops.min.6.log10.2.bedpe", header = F, sep = "\t")
gasperini_ranges <- GRanges(seqnames = gasperini$target_site.chr, ranges=IRanges(start=gasperini$target_gene.start, end=gasperini$target_gene.start+1, enh_id = seq(1,nrow(gasperini))))
gasperini_OE_ranges <- GRanges(seqnames = gasperini$target_site.chr, ranges=IRanges(start=gasperini$target_site.start, end=gasperini$target_site.stop), enh_id = seq(1,nrow(gasperini)))
loop_bait_ranges <- GRanges(seqnames = loops$V1, ranges=IRanges(start=loops$V2, end=loops$V3), enh_id = seq(1,nrow(loops)))
loop_OE_ranges <- GRanges(seqnames = loops$V4, ranges=IRanges(start=loops$V5-500, end=loops$V6+500), enh_id = seq(1,nrow(loops)))
gasperini$id <- seq(1,nrow(gasperini))
gasperini$overlap.loop <- gasperini$id %in% unique(as.data.frame(findOverlaps(gasperini_OE_ranges, loop_OE_ranges))[paste(as.data.frame(findOverlaps(gasperini_OE_ranges, loop_OE_ranges))$queryHits, as.data.frame(findOverlaps(gasperini_OE_ranges, loop_OE_ranges))$subjectHits) %in% paste(as.data.frame(findOverlaps(gasperini_ranges, loop_bait_ranges))$queryHits, as.data.frame(findOverlaps(gasperini_ranges, loop_bait_ranges))$subjectHits),]$queryHits)
#marking CTCF
CTCF.peaks <- read.csv("Documents/Nov2021_MCC/Chips/peaks/ENCFF002CEL.bed.gz", header = F, sep = "\t")
CTCF.peaks <- filter(CTCF.peaks, V5 > 250)
CTCF_overlaps <- as.data.frame(findOverlaps(GRanges(seqnames=gasperini$target_site.chr, ranges=IRanges(gasperini$target_site.start, gasperini$target_site.stop)), GRanges(seqnames = CTCF.peaks$V1, ranges=IRanges(CTCF.peaks$V2-1000,CTCF.peaks$V3+1000))))
gasperini$overlaps.CTCF <- gasperini$id %in% CTCF_overlaps$queryHits
gasperini$id <- NULL

#loading functionally tested E-P pairs from Fulco et al., 2019 (supplementary table 6a) and filtering out silencers, long-distance and out of proper baits
fulco <- read.csv("Documents/Jul2022_MChIPC/fig.4/fulco_functional_test_summary.txt", header = T, sep = "\t")
fulco$distance <- abs((fulco$start+fulco$end)/2 - fulco$Gene.TSS)
fulco <- filter(fulco, class != "promoter" & Fraction.change.in.gene.expr < 0 & Significant == TRUE & distance < 1000000 & distance > 5000)
baits <- read.csv("Documents/Jul2022_MChIPC/proper_baits.txt", header = T, sep = "\t")
fulco_ranges <- GRanges(seqnames = fulco$chr, ranges=IRanges(start=fulco$Gene.TSS, end=fulco$Gene.TSS+1, enh_id = seq(1,nrow(fulco))))
bait_ranges <- GRanges(seqnames = baits$bait_chr, ranges=IRanges(start=baits$bait_start, end=baits$bait_end), enh_id = seq(1,nrow(baits)))
fulco <- fulco[as.data.frame(findOverlaps(fulco_ranges, bait_ranges))$queryHits,]
#overlapping with MChIPC loops
fulco_ranges <- GRanges(seqnames = fulco$chr, ranges=IRanges(start=fulco$Gene.TSS, end=fulco$Gene.TSS+1, enh_id = seq(1,nrow(fulco))))
fulco_OE_ranges <- GRanges(seqnames = fulco$chr, ranges=IRanges(start=fulco$start, end=fulco$end), enh_id = seq(1,nrow(fulco)))
fulco$id <- seq(1,nrow(fulco))
fulco$overlap.loop <- fulco$id %in% unique(as.data.frame(findOverlaps(fulco_OE_ranges, loop_OE_ranges))[paste(as.data.frame(findOverlaps(fulco_OE_ranges, loop_OE_ranges))$queryHits, as.data.frame(findOverlaps(fulco_OE_ranges, loop_OE_ranges))$subjectHits) %in% paste(as.data.frame(findOverlaps(fulco_ranges, loop_bait_ranges))$queryHits, as.data.frame(findOverlaps(fulco_ranges, loop_bait_ranges))$subjectHits),]$queryHits)
#marking CTCF
CTCF.peaks <- read.csv("Documents/Nov2021_MCC/Chips/peaks/ENCFF002CEL.bed.gz", header = F, sep = "\t")
CTCF.peaks <- filter(CTCF.peaks, V5 > 250)
CTCF_overlaps <- as.data.frame(findOverlaps(GRanges(seqnames=fulco$chr, ranges=IRanges(fulco$start, fulco$end)), GRanges(seqnames = CTCF.peaks$V1, ranges=IRanges(CTCF.peaks$V2-1000,CTCF.peaks$V3+1000))))
fulco$overlaps.CTCF <- fulco$id %in% CTCF_overlaps$queryHits
fulco$id <- NULL

# choosing functionally unverified sequences
#from gasperini
gasperini_neg <- read.csv("Documents/Jul2022_MChIPC/fig.4/Gasperini_EP_pairs.txt", header = T, sep = "\t")
gasperini_neg$distance <- abs((gasperini_neg$target_site.start+gasperini_neg$target_site.stop)/2 - gasperini_neg$target_gene.start)
gasperini_neg <- filter(gasperini_neg, Significant == FALSE & distance < 1000000 & distance > 5000)
gasperini_neg_ranges <- GRanges(seqnames = gasperini_neg$target_site.chr, ranges=IRanges(start=gasperini_neg$target_gene.start, end=gasperini_neg$target_gene.start+1, enh_id = seq(1,nrow(gasperini_neg))))
gasperini_neg <- gasperini_neg[as.data.frame(findOverlaps(gasperini_neg_ranges, bait_ranges))$queryHits,]
gasperini_neg_ranges <- GRanges(seqnames = gasperini_neg$target_site.chr, ranges=IRanges(start=gasperini_neg$target_gene.start, end=gasperini_neg$target_gene.start+1, enh_id = seq(1,nrow(gasperini_neg))))
gasperini_neg_OE_ranges <- GRanges(seqnames = gasperini_neg$target_site.chr, ranges=IRanges(start=gasperini_neg$target_site.start, end=gasperini_neg$target_site.stop), enh_id = seq(1,nrow(gasperini_neg)))
gasperini_neg$id <- seq(1,nrow(gasperini_neg))
gasperini_neg$overlap.loop <- gasperini_neg$id %in% unique(as.data.frame(findOverlaps(gasperini_neg_OE_ranges, loop_OE_ranges))[paste(as.data.frame(findOverlaps(gasperini_neg_OE_ranges, loop_OE_ranges))$queryHits, as.data.frame(findOverlaps(gasperini_neg_OE_ranges, loop_OE_ranges))$subjectHits) %in% paste(as.data.frame(findOverlaps(gasperini_neg_ranges, loop_bait_ranges))$queryHits, as.data.frame(findOverlaps(gasperini_neg_ranges, loop_bait_ranges))$subjectHits),]$queryHits)
CTCF_overlaps <- as.data.frame(findOverlaps(GRanges(seqnames=gasperini_neg$target_site.chr, ranges=IRanges(gasperini_neg$target_site.start, gasperini_neg$target_site.stop)), GRanges(seqnames = CTCF.peaks$V1, ranges=IRanges(CTCF.peaks$V2-1000,CTCF.peaks$V3+1000))))
gasperini_neg$overlaps.CTCF <- gasperini_neg$id %in% CTCF_overlaps$queryHits
gasperini_neg$id <- NULL
# from fulco
fulco_neg <- read.csv("Documents/Jul2022_MChIPC/fig.4/fulco_functional_test_summary.txt", header = T, sep = "\t")
fulco_neg$distance <- abs((fulco_neg$start+fulco_neg$end)/2 - fulco_neg$Gene.TSS)
fulco_neg <- filter(fulco_neg, class != "promoter" & Significant == FALSE & distance < 1000000 & distance > 5000)
baits <- read.csv("Documents/Jul2022_MChIPC/proper_baits.txt", header = T, sep = "\t")
fulco_neg_ranges <- GRanges(seqnames = fulco_neg$chr, ranges=IRanges(start=fulco_neg$Gene.TSS, end=fulco_neg$Gene.TSS+1, enh_id = seq(1,nrow(fulco_neg))))
bait_ranges <- GRanges(seqnames = baits$bait_chr, ranges=IRanges(start=baits$bait_start, end=baits$bait_end), enh_id = seq(1,nrow(baits)))
fulco_neg <- fulco_neg[as.data.frame(findOverlaps(fulco_neg_ranges, bait_ranges))$queryHits,]
fulco_neg_OE_ranges <- GRanges(seqnames = fulco_neg$chr, ranges=IRanges(start=fulco_neg$start, end=fulco_neg$end), enh_id = seq(1,nrow(fulco_neg)))
fulco_neg$id <- seq(1,nrow(fulco_neg))
fulco_neg$overlap.loop <- fulco_neg$id %in% unique(as.data.frame(findOverlaps(fulco_neg_OE_ranges, loop_OE_ranges))[paste(as.data.frame(findOverlaps(fulco_neg_OE_ranges, loop_OE_ranges))$queryHits, as.data.frame(findOverlaps(fulco_neg_OE_ranges, loop_OE_ranges))$subjectHits) %in% paste(as.data.frame(findOverlaps(fulco_neg_ranges, loop_bait_ranges))$queryHits, as.data.frame(findOverlaps(fulco_neg_ranges, loop_bait_ranges))$subjectHits),]$queryHits)
CTCF_overlaps <- as.data.frame(findOverlaps(GRanges(seqnames=fulco_neg$chr, ranges=IRanges(fulco_neg$start, fulco_neg$end)), GRanges(seqnames = CTCF.peaks$V1, ranges=IRanges(CTCF.peaks$V2-1000,CTCF.peaks$V3+1000))))
fulco_neg$overlaps.CTCF <- fulco_neg$id %in% CTCF_overlaps$queryHits
fulco_neg$id <- NULL

#combining two datasets
gasperini <- gasperini[,c(1:4,6,9,10)]
fulco <- fulco[,c(1,2,3,5,18,26,27)]
colnames(fulco) <- colnames(gasperini)
gasperini <- rbind(gasperini,fulco)
# the same for negative
gasperini_neg <- gasperini_neg[,c(1:4,6,9,10)]
fulco_neg <- fulco_neg[,c(1,2,3,5,18,26,27)]
colnames(fulco_neg) <- colnames(gasperini_neg)
gasperini_neg <- rbind(gasperini_neg,fulco_neg)

rm(list=setdiff(ls(), c("gasperini", "gasperini_neg")))

# plotting heatmaps and average lineplots
# for loop positive enhancers
gasperini_ranges <- GRanges(seqnames = gasperini$target_site.chr, ranges=IRanges(start=gasperini$target_gene.start, end=gasperini$target_gene.start+1, enh_id = seq(1,nrow(gasperini))))
interactions <- read.csv("Documents/Jul2022_MChIPC/interaction calling/interactions.txt", header = T, sep = "\t")
interactions_ranges <- GRanges(seqnames = interactions$bait_chr, ranges=IRanges(start=interactions$bait_start, end=interactions$bait_end, enh_id = seq(1,nrow(interactions))))
interactions <- interactions[unique(as.data.frame(findOverlaps(interactions_ranges, gasperini_ranges))$queryHits),]
interactions_ranges <- GRanges(seqnames = interactions$bait_chr, ranges=IRanges(start=interactions$bait_start, end=interactions$bait_end, enh_id = seq(1,nrow(interactions))))
gasperini_ovrl <- filter(gasperini, overlap.loop==TRUE & overlaps.CTCF==FALSE)
gasperini_ovrl_heatmap <- tibble(dist_rank=seq(-19,19))
for (i in seq(nrow(gasperini_ovrl))){
  loop <- gasperini_ovrl[i,]
  loop_ranges <- GRanges(seqnames = loop$target_site.chr, ranges=IRanges(start=loop$target_gene.start, end=loop$target_gene.start+1, enh_id = seq(1,nrow(loop))))
  loop_interactions <- interactions[as.data.frame(findOverlaps(interactions_ranges, loop_ranges))$queryHits,]
  loop_interactions[loop_interactions$bait_start > loop_interactions$OE_start,]$dist_rank <- -1 * loop_interactions[loop_interactions$bait_start > loop_interactions$OE_start,]$dist_rank
  if (loop$target_gene.start > loop$target_site.start) {loop_interactions <- filter(loop_interactions, dist_rank>floor(((loop$target_site.start+loop$target_site.stop)/2 - bait_start)/250)-20 & dist_rank<floor(((loop$target_site.start+loop$target_site.stop)/2 - bait_start)/250)+20)} else {loop_interactions <- filter(loop_interactions, dist_rank>ceiling(((loop$target_site.start+loop$target_site.stop)/2 - bait_end)/250)-20 & dist_rank<ceiling(((loop$target_site.start+loop$target_site.stop)/2 - bait_end)/250)+20)}
  if (loop$target_gene.start > loop$target_site.start) {loop_interactions$dist_rank <- loop_interactions$dist_rank - floor(((loop$target_site.start+loop$target_site.stop)/2 - loop_interactions$bait_start)/250)} else {loop_interactions$dist_rank <- loop_interactions$dist_rank - ceiling(((loop$target_site.start+loop$target_site.stop)/2 - loop_interactions$bait_end)/250)}
  loop_interactions <- loop_interactions[,c(21,29)]
  colnames(loop_interactions)[2] <- paste0("loop_",i)
  gasperini_ovrl_heatmap <- left_join(gasperini_ovrl_heatmap, loop_interactions)}
gasperini_ovrl_heatmap[is.na(gasperini_ovrl_heatmap)] <- 0
gasperini_ovrl_heatmap <- as.data.frame(t(gasperini_ovrl_heatmap)[2:ncol(gasperini_ovrl_heatmap),])
gasperini_ovrl_heatmap <- arrange(gasperini_ovrl_heatmap, rowMeans(gasperini_ovrl_heatmap))
gasperini_ovrl_heatmap_mod <- gasperini_ovrl_heatmap
gasperini_ovrl_heatmap_mod[gasperini_ovrl_heatmap_mod>20] <- 20
heatmap(as.matrix(gasperini_ovrl_heatmap_mod) , Colv = NA, Rowv = NA, scale = "none", labRow = FALSE, labCol = FALSE)
ggplot() + geom_line(aes(x=seq(-19,19,1),y=colMeans(gasperini_ovrl_heatmap))) + theme_classic() + xlim(-19,19) + ylim(0,17) + xlab("")+ylab("")
# calculating positions of loops for PTGER3, BTG1 and ANTXR2
nrow(gasperini_ovrl_heatmap_mod) - which(rownames(gasperini_ovrl_heatmap_mod) %in% paste0("loop_",which(gasperini_ovrl$gene_short_name=="PTGER3")))
nrow(gasperini_ovrl_heatmap_mod) - which(rownames(gasperini_ovrl_heatmap_mod) %in% paste0("loop_",which(gasperini_ovrl$gene_short_name=="BTG1")))
nrow(gasperini_ovrl_heatmap_mod) - which(rownames(gasperini_ovrl_heatmap_mod) %in% paste0("loop_",which(gasperini_ovrl$gene_short_name=="ANTXR2")))
# for loop-negative enhancers
gasperini_non_ovrl <- filter(gasperini, overlap.loop==FALSE & overlaps.CTCF==FALSE)
gasperini_non_ovrl_heatmap <- tibble(dist_rank=seq(-19,19))
for (i in seq(nrow(gasperini_non_ovrl))){
  loop <- gasperini_non_ovrl[i,]
  loop_ranges <- GRanges(seqnames = loop$target_site.chr, ranges=IRanges(start=loop$target_gene.start, end=loop$target_gene.start+1, enh_id = seq(1,nrow(loop))))
  loop_interactions <- interactions[as.data.frame(findOverlaps(interactions_ranges, loop_ranges))$queryHits,]
  loop_interactions[loop_interactions$bait_start > loop_interactions$OE_start,]$dist_rank <- -1 * loop_interactions[loop_interactions$bait_start > loop_interactions$OE_start,]$dist_rank
  if (loop$target_gene.start > loop$target_site.start) {loop_interactions <- filter(loop_interactions, dist_rank>floor(((loop$target_site.start+loop$target_site.stop)/2 - bait_start)/250)-20 & dist_rank<floor(((loop$target_site.start+loop$target_site.stop)/2 - bait_start)/250)+20)} else {loop_interactions <- filter(loop_interactions, dist_rank>ceiling(((loop$target_site.start+loop$target_site.stop)/2 - bait_end)/250)-20 & dist_rank<ceiling(((loop$target_site.start+loop$target_site.stop)/2 - bait_end)/250)+20)}
  if (loop$target_gene.start > loop$target_site.start) {loop_interactions$dist_rank <- loop_interactions$dist_rank - floor(((loop$target_site.start+loop$target_site.stop)/2 - loop_interactions$bait_start)/250)} else {loop_interactions$dist_rank <- loop_interactions$dist_rank - ceiling(((loop$target_site.start+loop$target_site.stop)/2 - loop_interactions$bait_end)/250)}
  loop_interactions <- loop_interactions[,c(21,29)]
  colnames(loop_interactions)[2] <- paste0("loop_",i)
  gasperini_non_ovrl_heatmap <- left_join(gasperini_non_ovrl_heatmap, loop_interactions)}
gasperini_non_ovrl_heatmap[is.na(gasperini_non_ovrl_heatmap)] <- 0
gasperini_non_ovrl_heatmap <- as.data.frame(t(gasperini_non_ovrl_heatmap)[2:ncol(gasperini_non_ovrl_heatmap),])
gasperini_non_ovrl_heatmap <- arrange(gasperini_non_ovrl_heatmap, rowMeans(gasperini_non_ovrl_heatmap))
gasperini_non_ovrl_heatmap_mod <- gasperini_non_ovrl_heatmap
gasperini_non_ovrl_heatmap_mod[gasperini_non_ovrl_heatmap_mod>20] <- 20
heatmap(as.matrix(gasperini_non_ovrl_heatmap_mod) , Colv = NA, Rowv = NA, scale = "none", labRow = FALSE, labCol = FALSE)
ggplot() + geom_line(aes(x=seq(-19,19,1),y=colMeans(gasperini_non_ovrl_heatmap))) + theme_classic() + xlim(-19,19) + ylim(0,17) + xlab("")+ylab("")
# calculating positions of loops for CITED2, SEMA7A and CAT
nrow(gasperini_non_ovrl_heatmap_mod) - which(rownames(gasperini_non_ovrl_heatmap_mod) %in% paste0("loop_",which(gasperini_non_ovrl$gene_short_name=="CITED2")))
nrow(gasperini_non_ovrl_heatmap_mod) - which(rownames(gasperini_non_ovrl_heatmap_mod) %in% paste0("loop_",which(gasperini_non_ovrl$gene_short_name=="SEMA7A")))
nrow(gasperini_non_ovrl_heatmap_mod) - which(rownames(gasperini_non_ovrl_heatmap_mod) %in% paste0("loop_",which(gasperini_non_ovrl$gene_short_name=="CAT")))

# plotting distance distribution
gasperini$distance <- abs((gasperini$target_site.start+gasperini$target_site.stop)/2 - gasperini$target_gene.start)
ggplot(gasperini) + geom_boxplot(aes(overlap.loop, log10(distance))) +  geom_jitter(aes(overlap.loop, log10(distance)), width=0.25)

rm(list=setdiff(ls(), c("gasperini", "gasperini_neg")))

# for loop-positive non-enhancers
gasperini_neg_ranges <- GRanges(seqnames = gasperini_neg$target_site.chr, ranges=IRanges(start=gasperini_neg$target_gene.start, end=gasperini_neg$target_gene.start+1, enh_id = seq(1,nrow(gasperini_neg))))
interactions <- read.csv("Documents/Jul2022_MChIPC/pairs/interactions.txt", header = T, sep = "\t")
interactions_ranges <- GRanges(seqnames = interactions$bait_chr, ranges=IRanges(start=interactions$bait_start, end=interactions$bait_end, enh_id = seq(1,nrow(interactions))))
interactions <- interactions[unique(as.data.frame(findOverlaps(interactions_ranges, gasperini_neg_ranges))$queryHits),]
interactions_ranges <- GRanges(seqnames = interactions$bait_chr, ranges=IRanges(start=interactions$bait_start, end=interactions$bait_end, enh_id = seq(1,nrow(interactions))))
gasperini_neg_ovrl <- filter(gasperini_neg, overlap.loop.ns==TRUE  & overlaps.CTCF==FALSE)
gasperini_neg_ovrl_heatmap <- tibble(dist_rank=seq(-19,19))
for (i in seq(nrow(gasperini_neg_ovrl))){
  loop <- gasperini_neg_ovrl[i,]
  loop_ranges <- GRanges(seqnames = loop$target_site.chr, ranges=IRanges(start=loop$target_gene.start, end=loop$target_gene.start+1, enh_id = seq(1,nrow(loop))))
  loop_interactions <- interactions[as.data.frame(findOverlaps(interactions_ranges, loop_ranges))$queryHits,]
  loop_interactions[loop_interactions$bait_start > loop_interactions$OE_start,]$dist_rank <- -1 * loop_interactions[loop_interactions$bait_start > loop_interactions$OE_start,]$dist_rank
  if (loop$target_gene.start > loop$target_site.start) {loop_interactions <- filter(loop_interactions, dist_rank>floor(((loop$target_site.start+loop$target_site.stop)/2 - bait_start)/250)-20 & dist_rank<floor(((loop$target_site.start+loop$target_site.stop)/2 - bait_start)/250)+20)} else {loop_interactions <- filter(loop_interactions, dist_rank>ceiling(((loop$target_site.start+loop$target_site.stop)/2 - bait_end)/250)-20 & dist_rank<ceiling(((loop$target_site.start+loop$target_site.stop)/2 - bait_end)/250)+20)}
  if (loop$target_gene.start > loop$target_site.start) {loop_interactions$dist_rank <- loop_interactions$dist_rank - floor(((loop$target_site.start+loop$target_site.stop)/2 - loop_interactions$bait_start)/250)} else {loop_interactions$dist_rank <- loop_interactions$dist_rank - ceiling(((loop$target_site.start+loop$target_site.stop)/2 - loop_interactions$bait_end)/250)}
  loop_interactions <- loop_interactions[,c(21,29)]
  colnames(loop_interactions)[2] <- paste0("loop_",i)
  gasperini_neg_ovrl_heatmap <- left_join(gasperini_neg_ovrl_heatmap, loop_interactions)}
gasperini_neg_ovrl_heatmap[is.na(gasperini_neg_ovrl_heatmap)] <- 0
gasperini_neg_ovrl_heatmap <- as.data.frame(t(gasperini_neg_ovrl_heatmap)[2:ncol(gasperini_neg_ovrl_heatmap),])
gasperini_neg_ovrl_heatmap <- arrange(gasperini_neg_ovrl_heatmap, rowMeans(gasperini_neg_ovrl_heatmap))
gasperini_neg_ovrl_heatmap_mod <- gasperini_neg_ovrl_heatmap
gasperini_neg_ovrl_heatmap_mod[gasperini_neg_ovrl_heatmap_mod>20] <- 20
heatmap(as.matrix(gasperini_neg_ovrl_heatmap_mod) , Colv = NA, Rowv = NA, scale = "none", labRow = FALSE, labCol = FALSE)
ggplot() + geom_line(aes(x=seq(-19,19,1),y=colMeans(gasperini_neg_ovrl_heatmap))) + theme_classic() + xlim(-19,19) + ylim(0,17) + xlab("")+ylab("")

# for loop-negative non-enhancers (very slow!!!, that's why for now take nrow(gasperini_neg_ovrl))
gasperini_neg_non_ovrl <- filter(gasperini_neg, overlap.loop.ns==FALSE  & overlaps.CTCF==FALSE)
gasperini_neg_non_ovrl <- gasperini_neg_non_ovrl[sample(nrow(gasperini_neg_non_ovrl), nrow(gasperini_neg_ovrl)),]
gasperini_neg_non_ovrl_heatmap <- tibble(dist_rank=seq(-19,19))
for (i in seq(nrow(gasperini_neg_non_ovrl))){
  loop <- gasperini_neg_non_ovrl[i,]
  loop_ranges <- GRanges(seqnames = loop$target_site.chr, ranges=IRanges(start=loop$target_gene.start, end=loop$target_gene.start+1, enh_id = seq(1,nrow(loop))))
  loop_interactions <- interactions[as.data.frame(findOverlaps(interactions_ranges, loop_ranges))$queryHits,]
  loop_interactions[loop_interactions$bait_start > loop_interactions$OE_start,]$dist_rank <- -1 * loop_interactions[loop_interactions$bait_start > loop_interactions$OE_start,]$dist_rank
  if (loop$target_gene.start > loop$target_site.start) {loop_interactions <- filter(loop_interactions, dist_rank>floor(((loop$target_site.start+loop$target_site.stop)/2 - bait_start)/250)-20 & dist_rank<floor(((loop$target_site.start+loop$target_site.stop)/2 - bait_start)/250)+20)} else {loop_interactions <- filter(loop_interactions, dist_rank>ceiling(((loop$target_site.start+loop$target_site.stop)/2 - bait_end)/250)-20 & dist_rank<ceiling(((loop$target_site.start+loop$target_site.stop)/2 - bait_end)/250)+20)}
  if (loop$target_gene.start > loop$target_site.start) {loop_interactions$dist_rank <- loop_interactions$dist_rank - floor(((loop$target_site.start+loop$target_site.stop)/2 - loop_interactions$bait_start)/250)} else {loop_interactions$dist_rank <- loop_interactions$dist_rank - ceiling(((loop$target_site.start+loop$target_site.stop)/2 - loop_interactions$bait_end)/250)}
  loop_interactions <- loop_interactions[,c(21,26)]
  colnames(loop_interactions)[2] <- paste0("loop_",i)
  gasperini_neg_non_ovrl_heatmap <- left_join(gasperini_neg_non_ovrl_heatmap, loop_interactions)}
gasperini_neg_non_ovrl_heatmap[is.na(gasperini_neg_non_ovrl_heatmap)] <- 0
gasperini_neg_non_ovrl_heatmap <- as.data.frame(t(gasperini_neg_non_ovrl_heatmap)[2:ncol(gasperini_neg_non_ovrl_heatmap),])
gasperini_neg_non_ovrl_heatmap <- arrange(gasperini_neg_non_ovrl_heatmap, rowMeans(gasperini_neg_non_ovrl_heatmap))
gasperini_neg_non_ovrl_heatmap_mod <- gasperini_neg_non_ovrl_heatmap
gasperini_neg_non_ovrl_heatmap_mod[gasperini_neg_non_ovrl_heatmap_mod>20] <- 20
heatmap(as.matrix(gasperini_neg_non_ovrl_heatmap_mod) , Colv = NA, Rowv = NA, scale = "none", labRow = FALSE, labCol = FALSE)
ggplot() + geom_line(aes(x=seq(-19,19,1),y=colMeans(gasperini_neg_non_ovrl_heatmap))) + theme_classic() + xlim(-19,19) + ylim(0,17) + xlab("")+ylab("")

# addition of coverage information
interactions <- read.csv("Documents/Jul2022_MChIPC/pairs/interactions.txt", header = T, sep = "\t")
coverage <- distinct(interactions[,c(1:3,20)])
rm(interactions)
coverage_ranges <- GRanges(seqnames = coverage$bait_chr, ranges=IRanges(start=coverage$bait_start, end=coverage$bait_end), enh_id = seq(1,nrow(coverage)))
gasperini_ranges <- GRanges(seqnames = gasperini$target_site.chr, ranges=IRanges(start=gasperini$target_gene.start, end=gasperini$target_gene.start+1), enh_id = seq(1,nrow(gasperini)))
gasperini_neg_ranges <- GRanges(seqnames = gasperini_neg$target_site.chr, ranges=IRanges(start=gasperini_neg$target_gene.start, end=gasperini_neg$target_gene.start+1), enh_id = seq(1,nrow(gasperini_neg)))
gasperini$bait_coverage <- coverage[as.data.frame(findOverlaps(gasperini_ranges, coverage_ranges))[,2],]$bait_coverage
gasperini_neg$bait_coverage <- coverage[as.data.frame(findOverlaps(gasperini_neg_ranges, coverage_ranges))[,2],]$bait_coverage

