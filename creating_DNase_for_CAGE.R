library("dplyr")
library("tidyr")
library("GenomicRanges")

H3K4me3 <- read.csv("Documents/Jul2022_MChIPC/ChIPs/macs/binned_peaks.bed", header = F, sep = "\t")
DNase <- read.csv("Documents/Nov2021_MCC/Chips/peaks/ENCFF621ZJY.bed.gz", header = F, sep = "\t")
DNase <-filter(DNase, V7>200)

H3K4me3_ranges <- GRanges(seqnames = H3K4me3$V1, ranges=IRanges(start=H3K4me3$V2, end=H3K4me3$V3, enh_id = seq(1,nrow(H3K4me3))))
DNase_ranges <- GRanges(seqnames = DNase$V1, ranges=IRanges(start=DNase$V2, end=DNase$V3, enh_id = seq(1,nrow(DNase))))
DNase <- DNase[-(as.data.frame(findOverlaps(DNase_ranges,H3K4me3_ranges))[,1]),]

write.table(DNase, file="Documents/Jul2022_MChIPC/fig.S1/DNase_outside_promoters.bed", col.names = F, row.names = F, sep = "\t", quote = F)
