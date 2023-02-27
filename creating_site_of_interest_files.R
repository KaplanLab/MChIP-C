library('dplyr')
library('tidyr')
library('GenomicRanges')
# loading CTCF&DNase peaks as well as H3K4me3 mask
CTCF.peaks <- read.csv("Documents/Nov2021_MCC/Chips/peaks/ENCFF002CEL.bed.gz", header = F, sep = "\t")
DNase.peaks <- read.csv("Documents/Nov2021_MCC/Chips/peaks/ENCFF621ZJY.bed.gz", header = F, sep = "\t")
H3K4me3_mask <- read.csv("Documents/Jul2022_MChIPC/ChIPs/macs/H3K4me3.mask.bed", header = F, sep = "\t")
# filtering out weak peaks
CTCF.peaks <- filter(CTCF.peaks, V5 > 250)
DNase.peaks <- filter(DNase.peaks, V7 > 200)
# removing CTCF peaks overlaping H3K4me3
CTCF_ranges <- GRanges(seqnames = CTCF.peaks$V1, ranges=IRanges(start=CTCF.peaks$V2, end=CTCF.peaks$V3, enh_id = seq(1,nrow(CTCF.peaks))))
H3K4me3_ranges <- GRanges(seqnames = H3K4me3_mask$V1, ranges=IRanges(start=H3K4me3_mask$V2-5000, end=H3K4me3_mask$V3+5000, enh_id = seq(1,nrow(H3K4me3_mask))))
CTCF_wo_H3K4me3 <- CTCF.peaks[-unique(as.data.frame(findOverlaps(CTCF_ranges,H3K4me3_ranges))[,1]),]
# removing DNase peaks overlapping CTCF or H3K4me3
DNase_ranges <- GRanges(seqnames = DNase.peaks$V1, ranges=IRanges(start=DNase.peaks$V2, end=DNase.peaks$V3, enh_id = seq(1,nrow(DNase.peaks))))
DNase_wo_H3K4me3 <- DNase.peaks[-unique(c(as.data.frame(findOverlaps(DNase_ranges,H3K4me3_ranges))[,1], as.data.frame(findOverlaps(DNase_ranges,CTCF_ranges))[,1])),]
# writing these all
write.table(CTCF_wo_H3K4me3, file="Documents/Jul2022_MChIPC/fig.2/heatmaps/CTCF_of_interest.bed", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(DNase_wo_H3K4me3, file="Documents/Jul2022_MChIPC/fig.2/heatmaps/DNase_of_interest.bed", col.names = F, row.names = F, sep = "\t", quote = F)

# reading sorted regions (after plotting heatmap)
sorted <- read.csv("Documents/Jul2022_MChIPC/fig.2/heatmaps/sorted.regions.bed", header = T, sep = "\t")
sorted_CTCF <- filter(sorted, grepl("CTCF",deepTools_group))
sorted_DNase <- filter(sorted, grepl("DNase",deepTools_group))
