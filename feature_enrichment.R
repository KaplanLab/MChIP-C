# algorythm for assessing enrichment of feature (ChIP-seq peaks) in OE of loop-bedpe file
# 3 main inputs: (1) loop-bedpe file; (2) bed file with feature peaks; (3) bed file with MChIPC anchors
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
read_bed <- function(x){
  bed <- read.csv(paste0(folder,x[1,2]), header = F, sep = "\t")
  if(sum(bed$V8) < 0) {bed <- filter(bed, V7>25)} else {bed <- filter(bed, V8>20)}
}
# loading files with ChIP peaks
folder <- "Documents/Nov2021_MCC/Chips/peaks/"
list <- read.csv("Documents/Nov2021_MCC/Chips/list_of_beds_and_bigwigs.txt", header = T, sep = "\t")[,1:2]
list <- split(list, list$target)
beds <- lapply(list, read_bed)
loops <-read.csv("Documents/Jul2022_MChIPC/interaction calling/loops.min.6.log10.2.bedpe", header = F, sep = "\t")
loops$V7 <- 0
loops$V7[loops$V2>loops$V5] <- -((loops$V2[loops$V2>loops$V5]-loops$V6[loops$V2>loops$V5])/250+1)
loops$V7[loops$V2<loops$V5] <- (loops$V5[loops$V2<loops$V5]-loops$V3[loops$V2<loops$V5])/250+1
# for the downstream analysis I'm actually using not baits, but all peaks (mask), to filter out potentially artificial loops
baits <- read.csv("Documents/Jul2022_MChIPC/ChIPs/macs/H3K4me3.mask.bed", header = F, sep = "\t")[,c(1:3)]
# discarding promoter-promoter loops
loops_OE_ranges <- GRanges(seqnames = loops$V4, ranges=IRanges(start=loops$V5, end=loops$V6, enh_id = seq(1,nrow(loops))))
bait_ranges <- GRanges(seqnames = baits$V1, ranges=IRanges(start=baits$V2, end=baits$V3, enh_id = seq(1,nrow(baits))))
PIRs <- unique(loops[-(as.data.frame(findOverlaps(loops_OE_ranges,bait_ranges))[,1]),c(4:6)])
write.table(PIRs, file="Documents/Jul2022_MChIPC/fig.3/PIRs.bed", col.names = F, row.names = F, sep = "\t", quote = F)
# here using Homer 'findMotifsGenome.pl' to find motif enrichment
# findMotifsGenome.pl PIRs.bed hg19 homer_output/ -p 3
PIR_ranges <- GRanges(seqnames = PIRs$V4, ranges=IRanges(start=PIRs$V5, end=PIRs$V6, enh_id = seq(1,nrow(PIRs))))
peak_ranges <- lapply(beds, function(bed){GRanges(seqnames = bed$V1, ranges=IRanges(start=bed$V2, end=bed$V3, enh_id = seq(1,nrow(bed))))})
number_of_overlaps <- sapply(peak_ranges, function(range){length(unique(as.data.frame(findOverlaps(PIR_ranges, range))[,1]))})
# sampling random 100 times
set.seed(42)
randomized_numbers <- c()
for (i in seq(100)){
  randomized_OE <- loops[,1:6]
  randomized_OE$V7 <- sample(loops$V7)
  randomized_pos <- filter(randomized_OE, V7 > 0)
  randomized_neg <- filter(randomized_OE, V7 < 0)
  randomized_pos$V6 <- randomized_pos$V3 + randomized_pos$V7 * 250
  randomized_pos$V5 <- randomized_pos$V6 - 250
  randomized_neg$V5 <- randomized_neg$V2 + randomized_neg$V7 *250
  randomized_neg$V6 <- randomized_neg$V5 + 250
  PIRs_randomized <- unique(bind_rows(randomized_pos, randomized_neg)[,c(4:6)])
  PIR_randomized_ranges <- GRanges(seqnames = PIRs_randomized$V4, ranges=IRanges(start=PIRs_randomized$V5, end=PIRs_randomized$V6, enh_id = seq(1,nrow(PIRs_randomized))))
  PIRs_randomized <- sample_n(PIRs_randomized[-unique(as.data.frame(findOverlaps(PIR_randomized_ranges,bait_ranges))[,1]),], nrow(PIRs))
  PIR_randomized_ranges <- GRanges(seqnames = PIRs_randomized$V4, ranges=IRanges(start=PIRs_randomized$V5, end=PIRs_randomized$V6, enh_id = seq(1,nrow(PIRs_randomized))))
  randomized_number_of_overlaps <- sapply(peak_ranges,function(range){length(unique(as.data.frame(findOverlaps(PIR_randomized_ranges, range))[,1]))})
  randomized_numbers <- bind_rows(randomized_numbers, randomized_number_of_overlaps)}

# summarizing the data and plotting the most enriched features
summary <- data_frame(ChIP=colnames(randomized_numbers), overlaps = number_of_overlaps, randomized_mean = colMeans(randomized_numbers), sd = as.numeric(summarise_if(randomized_numbers,is.numeric, sd)))
summary$enrichment <- summary$overlaps / summary$randomized_mean
melted <- melt(summary[,1:4], c("ChIP","sd"))
for (i in seq(nrow(melted))){
  if (melted[i,3]=="overlaps"){melted[i,2] <- 0}}
# reading output of Homer
motif_enrichment <- read.csv("Documents/Jul2022_MChIPC/fig.3/homer_output/knownResults.txt", header = T, sep = "\t")
motifs <- read.csv("Downloads/homer/motifs/table.txt", header = T, sep = "\t")
motifs["Symbol"][motifs["Factor.Name"]=="Gata1"] <- "GATA1"
motif_enrichment <- left_join(motif_enrichment,motifs[c("Name","Symbol")], by=c("Motif.Name"="Name"))
summary <- left_join(summary, motif_enrichment[c("Symbol","Log.P.value")], by=c("ChIP"="Symbol"))
summary["Log.P.value"] <- -(summary["Log.P.value"])
summary$sd <- round(summary$sd,2)
summary$enrichment <- round(log2(summary$enrichment),2)
write.table(summary, file="Documents/Jul2022_MChIPC/fig.3/enrichment_summary.txt", col.names = T, row.names = F, sep = "\t", quote = F)
summary$Log.P.value[summary$Log.P.value > 50] <- 50
summary <- summary %>% arrange(!is.na(Log.P.value), Log.P.value)
summary <- filter(summary, overlaps >100)
mycolors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(50)
ggplot(summary, aes(overlaps*100/nrow(PIRs),enrichment)) + geom_point(aes(colour=Log.P.value), size=3) +  geom_text(aes(label=ifelse(ChIP %in% c("H3K4me1", "H3K27ac", "CTCF","SMC3", "RAD21", "ARID1B","H3K27me3","EP300", "DPF2","ZNF143"),as.character(ChIP),'')),hjust=0,vjust=0) + scale_colour_gradientn(colours=mycolors, limits=c(0, 50)) + xlab("overlaps (%)") + ylab("log2( Observed / Expected )")
