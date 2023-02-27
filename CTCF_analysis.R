library(dplyr)
library(tidyr)
library(GenomicRanges)
library(ggplot2)
# loading loops
loops <-read.csv("Documents/Jul2022_MChIPC/interaction calling//loops.min.6.log10.2.bedpe", header = F, sep = "\t")
loops$V7 <- 0
loops$V7[loops$V2>loops$V5] <- -((loops$V2[loops$V2>loops$V5]-loops$V6[loops$V2>loops$V5])/250+1)
loops$V7[loops$V2<loops$V5] <- (loops$V5[loops$V2<loops$V5]-loops$V3[loops$V2<loops$V5])/250+1
# selecting P-PIR loops
baits <- read.csv("Documents/Jul2022_MChIPC/ChIPs/macs/binned_peaks.bed", header = F, sep = "\t")[,c(1:3)]
bait_ranges <- GRanges(seqnames = baits$V1, ranges=IRanges(start=baits$V2, end=baits$V3, enh_id = seq(1,nrow(baits))))
loops_OE_ranges <- GRanges(seqnames = loops$V4, ranges=IRanges(start=loops$V5, end=loops$V6, enh_id = seq(1,nrow(loops))))
loops <- loops[-(as.data.frame(findOverlaps(loops_OE_ranges,bait_ranges))[,1]),]
# separating loops with CTCF in PIR and either with or without CTCF in promoter
CTCF_peaks <- read.csv("Documents/Nov2021_MCC/Chips/peaks/ENCFF002CEL.bed.gz", header = F, sep = "\t") # all peaks
CTCF_peaks_s <- filter(CTCF_peaks, V7 > 25) # strictly speaking peaks
CTCF_ranges_l <- GRanges(seqnames = CTCF_peaks$V1, ranges=IRanges(start=CTCF_peaks$V2, end=CTCF_peaks$V3, enh_id = seq(1,nrow(CTCF_peaks))))
CTCF_ranges_s <- GRanges(seqnames = CTCF_peaks_s$V1, ranges=IRanges(start=CTCF_peaks_s$V2, end=CTCF_peaks_s$V3, enh_id = seq(1,nrow(CTCF_peaks_s))))
PIR_ranges <- GRanges(seqnames = loops$V4, ranges=IRanges(start=loops$V5, end=loops$V6, enh_id = seq(1,nrow(loops))))
loops_PIR_CTCF <- loops[unique(as.data.frame(findOverlaps(PIR_ranges,CTCF_ranges_s))[,1]),]
loop_origin_ranges <- GRanges(seqnames = loops_PIR_CTCF$V1, ranges=IRanges(start=loops_PIR_CTCF$V2, end=loops_PIR_CTCF$V3, enh_id = seq(1,nrow(loops_PIR_CTCF))))
loops_PIR_only_CTCF <- loops_PIR_CTCF[-(unique(as.data.frame(findOverlaps(loop_origin_ranges,CTCF_ranges_l))[,1])),]
loops_PIR_and_bait_CTCF <- loops_PIR_CTCF[(unique(as.data.frame(findOverlaps(loop_origin_ranges,CTCF_ranges_l))[,1])),]
CTCF_peaks_in_baits <- CTCF_peaks[(unique(as.data.frame(findOverlaps(CTCF_ranges_l,loop_origin_ranges))[,1])),]
# compairing with HiC loops
HiC_loops <- read.csv("Documents/Jul2022_MChIPC/ABC-Enhancer-Gene-Prediction/reference/GSE63525_K562_HiCCUPS_looplist_with_motifs.txt.gz", header = T, sep = "\t")
HiC_loops <- filter(HiC_loops, !is.na(motif_x1) & !is.na(motif_x2))
HiC_loops$chr1 <- paste0("chr",HiC_loops$chr1)
HiC_loops$chr2 <- paste0("chr",HiC_loops$chr2)
HiC_loops_x_ranges <- GRanges(seqnames = HiC_loops$chr1, ranges=IRanges(start=HiC_loops$x1, end=HiC_loops$x2, enh_id = seq(1,nrow(HiC_loops))))
HiC_loops_y_ranges <- GRanges(seqnames = HiC_loops$chr2, ranges=IRanges(start=HiC_loops$y1, end=HiC_loops$y2, enh_id = seq(1,nrow(HiC_loops))))
CTCF_bait_ranges <- GRanges(seqnames = loops_PIR_and_bait_CTCF$V1, ranges=IRanges(start=loops_PIR_and_bait_CTCF$V2, end=loops_PIR_and_bait_CTCF$V3, enh_id = seq(1,nrow(loops_PIR_and_bait_CTCF))))
CTCF_PIR_ranges <- GRanges(seqnames = loops_PIR_and_bait_CTCF$V4, ranges=IRanges(start=loops_PIR_and_bait_CTCF$V5, end=loops_PIR_and_bait_CTCF$V6, enh_id = seq(1,nrow(loops_PIR_and_bait_CTCF))))
corresp <- left_join(unique(as.data.frame(findOverlaps(CTCF_bait_ranges,HiC_loops_x_ranges))), unique(as.data.frame(findOverlaps(CTCF_PIR_ranges,HiC_loops_y_ranges))), by="queryHits") %>% filter(subjectHits.x==subjectHits.y) %>%
  bind_rows(left_join(unique(as.data.frame(findOverlaps(CTCF_bait_ranges,HiC_loops_y_ranges))), unique(as.data.frame(findOverlaps(CTCF_PIR_ranges,HiC_loops_x_ranges))), by="queryHits") %>% filter(subjectHits.x==subjectHits.y))
length(unique(corresp$queryHits))
# now finding motifs in bait CTCF peaks and in PIRs
write.table(unique(loops_PIR_CTCF[,4:6]), file="Documents/Jul2022_MChIPC/fig.3/PIRs_CTCF.bed", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(unique(CTCF_peaks_in_baits[,1:3]), file="Documents/Jul2022_MChIPC/fig.3/baits_CTCF.bed", col.names = F, row.names = F, sep = "\t", quote = F)
# here using Homer 'annotatePeaks.pl' to find motif enrichment
# annotatePeaks.pl PIRs_CTCF.bed hg19 -nogene -noann -m ../../../Downloads/homer/motifs/ctcf.motif > motifs_in_PIRs_n
# annotatePeaks.pl baits_CTCF.bed hg19 -nogene -noann -m ../../../Downloads/homer/motifs/ctcf.motif > motifs_in_baits_n
# loop origins CTCF analysis
loop_origins_with_CTCF <- read.csv("Documents/Jul2022_MChIPC/fig.3/motifs_in_baits_n", header = T, sep = "\t")
loop_origins_with_CTCF <- unique(loop_origins_with_CTCF)
loop_origins_with_CTCF_a <- loop_origins_with_CTCF[!(grepl("\\-",loop_origins_with_CTCF$X.sequence.strand.conservation.) & grepl("\\+",loop_origins_with_CTCF$X.sequence.strand.conservation.) | loop_origins_with_CTCF$X.sequence.strand.conservation. %in% ""),]
loop_origins_with_CTCF_a$bait_CTCF_strand <- NA
for (i in seq(nrow(loop_origins_with_CTCF_a))){ if (grepl("\\+",loop_origins_with_CTCF_a[i,4])) {loop_origins_with_CTCF_a[i,5] <- "+"} else {loop_origins_with_CTCF_a[i,5] <- "-"}}
loop_origins_with_CTCF_b <- loop_origins_with_CTCF[grepl("\\-",loop_origins_with_CTCF$X.sequence.strand.conservation.) & grepl("\\+",loop_origins_with_CTCF$X.sequence.strand.conservation.),]
loop_origins_with_CTCF_b$bait_CTCF_strand <- "both"
loop_origins_with_CTCF <- bind_rows(loop_origins_with_CTCF_a[,c(1:3,5)],loop_origins_with_CTCF_b[,c(1:3,5)])
bait_ranges <- GRanges(seqnames = loop_origins_with_CTCF$Chr, ranges=IRanges(start=loop_origins_with_CTCF$Start, end=loop_origins_with_CTCF$End, enh_id = seq(1,nrow(loop_origins_with_CTCF))))
loop_origin_ranges <- GRanges(seqnames = loops_PIR_and_bait_CTCF$V1, ranges=IRanges(start=loops_PIR_and_bait_CTCF$V2, end=loops_PIR_and_bait_CTCF$V3, enh_id = seq(1,nrow(loops_PIR_and_bait_CTCF))))
loops_PIR_and_bait_CTCF$queryHits <- seq(1,nrow(loops_PIR_and_bait_CTCF))
loop_origins_with_CTCF$subjectHits <- seq(1,nrow(loop_origins_with_CTCF))
loops_PIR_and_bait_CTCF <- left_join(loops_PIR_and_bait_CTCF, as.data.frame(findOverlaps(loop_origin_ranges,bait_ranges)))
loops_PIR_and_bait_CTCF <- left_join(loops_PIR_and_bait_CTCF,loop_origins_with_CTCF[,c(4,5)])
# PIR CTCF analysis
loop_PIRs_with_CTCF <- read.csv("Documents/Jul2022_MChIPC/fig.3/motifs_in_PIRs_n", header = T, sep = "\t")
loop_PIRs_with_CTCF <- unique(loop_PIRs_with_CTCF)
loop_PIRs_with_CTCF_a <- loop_PIRs_with_CTCF[!(grepl("\\-",loop_PIRs_with_CTCF$X.sequence.strand.conservation.) & grepl("\\+",loop_PIRs_with_CTCF$X.sequence.strand.conservation.) | loop_PIRs_with_CTCF$X.sequence.strand.conservation. %in% ""),]
loop_PIRs_with_CTCF_a$PIR_CTCF_strand <- NA
for (i in seq(nrow(loop_PIRs_with_CTCF_a))){ if (grepl("\\+",loop_PIRs_with_CTCF_a[i,4])) {loop_PIRs_with_CTCF_a[i,5] <- "+"} else {loop_PIRs_with_CTCF_a[i,5] <- "-"}}
loop_PIRs_with_CTCF_b <- loop_PIRs_with_CTCF[grepl("\\-",loop_PIRs_with_CTCF$X.sequence.strand.conservation.) & grepl("\\+",loop_PIRs_with_CTCF$X.sequence.strand.conservation.),]
loop_PIRs_with_CTCF_b$PIR_CTCF_strand <- "both"
loop_PIRs_with_CTCF <- bind_rows(loop_PIRs_with_CTCF_a[,c(1:3,5)],loop_PIRs_with_CTCF_b[,c(1:3,5)])
loops_PIR_and_bait_CTCF <- left_join(loops_PIR_and_bait_CTCF, loop_PIRs_with_CTCF[,c(1,3,4)], by=c('V4'='Chr', 'V6'='End'))
loops_PIR_only_CTCF <- left_join(loops_PIR_only_CTCF, loop_PIRs_with_CTCF[,c(1,3,4)], by=c('V4'='Chr', 'V6'='End'))
# estimating number of motifs oriented towards and away from promoters in loops with CTCF in PIR only subset
nrow(filter(loops_PIR_only_CTCF, (V7<0&PIR_CTCF_strand=="-" | V7>0&PIR_CTCF_strand=="+" ))) # towards
nrow(filter(loops_PIR_only_CTCF, (V7<0&PIR_CTCF_strand=="+" | V7>0&PIR_CTCF_strand=="-" ))) # away
# counting number of loops with CTCF motif only in PIR (while CTCF is bound to both bait and PIR)
length(unique(filter(loops_PIR_and_bait_CTCF, is.na(bait_CTCF_strand) & !(is.na(PIR_CTCF_strand)))$queryHits)) # with "both"
nrow(filter(loops_PIR_and_bait_CTCF, is.na(bait_CTCF_strand) & !(is.na(PIR_CTCF_strand)) & (V7<0&PIR_CTCF_strand=="-" | V7>0&PIR_CTCF_strand=="+" ) )) # towards
nrow(filter(loops_PIR_and_bait_CTCF, is.na(bait_CTCF_strand) & !(is.na(PIR_CTCF_strand)) & (V7<0&PIR_CTCF_strand=="+" | V7>0&PIR_CTCF_strand=="-" ) )) # away
# counting number of loops with CTCF motif in both bait and PIR
loops_PIR_and_bait_motifs <- filter(loops_PIR_and_bait_CTCF, queryHits %in% (filter(loops_PIR_and_bait_CTCF, !(is.na(bait_CTCF_strand)) & !(is.na(PIR_CTCF_strand))) %>% group_by(queryHits) %>% summarise(N=length(unique(bait_CTCF_strand))) %>% filter(N==1))$queryHits & bait_CTCF_strand != "both" & PIR_CTCF_strand != "both")
loops_PIR_and_bait_motifs <- unique(loops_PIR_and_bait_motifs[,-9])
nrow(filter(loops_PIR_and_bait_motifs, (bait_CTCF_strand=="+"& PIR_CTCF_strand == "-" & V7 < 0) | (bait_CTCF_strand=="-"& PIR_CTCF_strand == "+" & V7 > 0))) # convergent
nrow(filter(loops_PIR_and_bait_motifs, (bait_CTCF_strand=="+"& PIR_CTCF_strand == "-" & V7 > 0) | (bait_CTCF_strand=="-"& PIR_CTCF_strand == "+" & V7 < 0))) # divergent
nrow(filter(loops_PIR_and_bait_motifs, (bait_CTCF_strand=="-"& PIR_CTCF_strand == "-" & V7 < 0) | (bait_CTCF_strand=="+"& PIR_CTCF_strand == "+" & V7 > 0))) # co-linear with PIR motif towards bait
nrow(filter(loops_PIR_and_bait_motifs, (bait_CTCF_strand=="-"& PIR_CTCF_strand == "-" & V7 > 0) | (bait_CTCF_strand=="+"& PIR_CTCF_strand == "+" & V7 < 0))) # co-linear with PIR motif away from bait

# adding expression
expression <- read.csv("Documents/Jul2022_MChIPC/bait_expression.txt", header = T, sep = "\t")
colnames(expression)[c(1,2,3,5)] <- c("V1","V2","V3","transcription_dir")
# plotting loops of CTCF-occupied promoters
# without motif in promoter
loops_PIR_motif <- filter(loops_PIR_and_bait_CTCF, is.na(bait_CTCF_strand) & !(is.na(PIR_CTCF_strand)) & (PIR_CTCF_strand=="-" | PIR_CTCF_strand=="+" ))
loops_PIR_motif <- left_join(loops_PIR_motif, expression[,c(1:3,5)])
loops_PIR_motif <- filter(loops_PIR_motif, !(is.na(transcription_dir)))
loops_PIR_motif_flipped <- loops_PIR_motif
loops_PIR_motif_flipped['V7'][loops_PIR_motif_flipped['transcription_dir']=="-"] <- -1*loops_PIR_motif_flipped['V7'][loops_PIR_motif_flipped['transcription_dir']=="-"]
loops_PIR_motif_flipped['PIR_CTCF_strand'][loops_PIR_motif_flipped['transcription_dir']=="-" & loops_PIR_motif_flipped['PIR_CTCF_strand']=="+"] <- 0
loops_PIR_motif_flipped['PIR_CTCF_strand'][loops_PIR_motif_flipped['transcription_dir']=="-" & loops_PIR_motif_flipped['PIR_CTCF_strand']=="-"] <- "+"
loops_PIR_motif_flipped['PIR_CTCF_strand'][loops_PIR_motif_flipped['transcription_dir']=="-" & loops_PIR_motif_flipped['PIR_CTCF_strand']==0] <- "-"
ggplot(filter(loops_PIR_motif_flipped, PIR_CTCF_strand %in% c("-","+")), aes(V7, fill=PIR_CTCF_strand, alpha=0.5)) + geom_histogram(bins = 100, position="identity") + xlim(-2000,2000)
# with motif
loops_PIR_and_bait_motifs_flipped <- loops_PIR_and_bait_motifs
loops_PIR_and_bait_motifs_flipped$PIR_CTCF_strand_flipped <- loops_PIR_and_bait_motifs_flipped$PIR_CTCF_strand
loops_PIR_and_bait_motifs_flipped$PIR_CTCF_strand_flipped[loops_PIR_and_bait_motifs_flipped$bait_CTCF_strand == "+" & loops_PIR_and_bait_motifs_flipped$PIR_CTCF_strand == "-"] <- "+"
loops_PIR_and_bait_motifs_flipped$PIR_CTCF_strand_flipped[loops_PIR_and_bait_motifs_flipped$bait_CTCF_strand == "+" & loops_PIR_and_bait_motifs_flipped$PIR_CTCF_strand == "+"] <- "-"
loops_PIR_and_bait_motifs_flipped['V7'][loops_PIR_and_bait_motifs_flipped['bait_CTCF_strand']=="+"] <- -1*loops_PIR_and_bait_motifs_flipped['V7'][loops_PIR_and_bait_motifs_flipped['bait_CTCF_strand']=="+"]
ggplot(filter(loops_PIR_and_bait_motifs_flipped, PIR_CTCF_strand_flipped %in% c("-","+")), aes(V7, fill=PIR_CTCF_strand_flipped, alpha=0.5)) + geom_histogram(bins = 100, position="identity") + xlim(-2000,2000)

# plotting loops of CTCF-less promoters
expression_ranges <- GRanges(seqnames = expression$V1, ranges=IRanges(start=expression$V2, end=expression$V3, enh_id = seq(1,nrow(expression))))
expression <- expression[-(unique(as.data.frame(findOverlaps(expression_ranges, CTCF_ranges_l))[,1])),]
expression <- expression %>% mutate(quartile=ntile(TPM,4))
loops_PIR_only_CTCF <- left_join(loops_PIR_only_CTCF, expression[,c(1,2,3,5,9)])
loops_PIR_only_CTCF <- filter(loops_PIR_only_CTCF, !(is.na(transcription_dir)))
# flipping CTCF position such as all transcription in the same direction
loops_PIR_only_CTCF_flipped <- loops_PIR_only_CTCF
loops_PIR_only_CTCF_flipped['V7'][loops_PIR_only_CTCF_flipped['transcription_dir']=="-"] <- -1*loops_PIR_only_CTCF_flipped['V7'][loops_PIR_only_CTCF_flipped['transcription_dir']=="-"]
loops_PIR_only_CTCF_flipped['PIR_CTCF_strand'][loops_PIR_only_CTCF_flipped['transcription_dir']=="-" & loops_PIR_only_CTCF_flipped['PIR_CTCF_strand']=="+"] <- 0
loops_PIR_only_CTCF_flipped['PIR_CTCF_strand'][loops_PIR_only_CTCF_flipped['transcription_dir']=="-" & loops_PIR_only_CTCF_flipped['PIR_CTCF_strand']=="-"] <- "+"
loops_PIR_only_CTCF_flipped['PIR_CTCF_strand'][loops_PIR_only_CTCF_flipped['transcription_dir']=="-" & loops_PIR_only_CTCF_flipped['PIR_CTCF_strand']==0] <- "-"
ggplot(filter(loops_PIR_only_CTCF_flipped, PIR_CTCF_strand %in% c("-","+")), aes(V7, fill=PIR_CTCF_strand, alpha=0.5)) + geom_histogram(bins = 100, position="identity") + xlim(-2000,2000)
# counting relative to the transcription direction
# upstream
nrow(filter(loops_PIR_only_CTCF_flipped, V7<0 & PIR_CTCF_strand=="-")) #towards
nrow(filter(loops_PIR_only_CTCF_flipped, V7<0 & PIR_CTCF_strand=="+")) #away
# downstream
nrow(filter(loops_PIR_only_CTCF_flipped, V7>0 & PIR_CTCF_strand=="+")) #towards
nrow(filter(loops_PIR_only_CTCF_flipped, V7>0 & PIR_CTCF_strand=="-")) #away

# plotting expression quartiles
ggplot(filter(loops_PIR_only_CTCF_flipped, quartile==1 & PIR_CTCF_strand %in% c("-","+")), aes(V7, fill=PIR_CTCF_strand, alpha=0.5)) + geom_histogram(bins = 100, position="identity") + xlim(-2000,2000)+ylim(0,70)
ggplot(filter(loops_PIR_only_CTCF_flipped, quartile==2 & PIR_CTCF_strand %in% c("-","+")), aes(V7, fill=PIR_CTCF_strand, alpha=0.5)) + geom_histogram(bins = 100, position="identity") + xlim(-2000,2000)+ylim(0,70)
ggplot(filter(loops_PIR_only_CTCF_flipped, quartile==3 & PIR_CTCF_strand %in% c("-","+")), aes(V7, fill=PIR_CTCF_strand, alpha=0.5)) + geom_histogram(bins = 100, position="identity") + xlim(-2000,2000)+ylim(0,70)
ggplot(filter(loops_PIR_only_CTCF_flipped, quartile==4 & PIR_CTCF_strand %in% c("-","+")), aes(V7, fill=PIR_CTCF_strand, alpha=0.5)) + geom_histogram(bins = 100, position="identity") + xlim(-2000,2000)+ylim(0,70)
