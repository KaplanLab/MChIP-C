library("dplyr")
library("tidyr")
library("GenomicRanges")
library("ggplot2")
interactions <- read.csv("Documents/Jul2022_MChIPC/interaction calling/interactions.txt", header = T, sep = "\t")
interactions <- interactions[,1:10]
interactions$distance <- log10(abs((interactions$bait_start+interactions$bait_end)/2 - (interactions$OE_start+interactions$OE_end)/2))
interactions <- filter(interactions, distance < 7)
baits <- read.csv("Documents/Jul2022_MChIPC/ChIPs/macs/H3K4me3.mask.bed", header = F, sep = "\t")[,c(1:3)]
interactions_ranges <- GRanges(seqnames = interactions$OE_chr, ranges=IRanges(start=interactions$OE_start, end=interactions$OE_end, enh_id = seq(1,nrow(interactions))))
bait_ranges <- GRanges(seqnames = baits$V1, ranges=IRanges(start=baits$V2, end=baits$V3, enh_id = seq(1,nrow(baits))))
interactions <- interactions[-(unique(as.data.frame(findOverlaps(interactions_ranges,bait_ranges))[,1])),]
interactions <- interactions[,7:11]
interactions <- interactions %>% mutate(bin = cut(distance, breaks=seq(3.5,6.5,0.1), labels=FALSE))
interactions$sum <- interactions$N_rep0 + interactions$N_rep3 + interactions$N_rep4 + interactions$N_rep5
summary <- interactions %>% group_by(bin) %>% summarise(N_rep0 = sum(N_rep0), N_rep3 = sum(N_rep3), N_rep4 = sum(N_rep4), N_rep5 = sum(N_rep5), N_sum=sum(sum))
summary$N_rep0 <- 100*summary$N_rep0/sum(summary$N_rep0)
summary$N_rep3 <- 100*summary$N_rep3/sum(summary$N_rep3)
summary$N_rep4 <- 100*summary$N_rep4/sum(summary$N_rep4)
summary$N_rep5 <- 100*summary$N_rep5/sum(summary$N_rep5)
summary$N_sum <- 100*summary$N_sum/sum(summary$N_sum)
summary$distance <- 10^(3.5+summary$bin*0.1)
ggplot(summary) + geom_smooth(aes(x=distance, y=N_rep0), color=1, span = 0.5) + geom_smooth(aes(x=distance, y=N_rep3), color=2, span = 0.5) +geom_smooth(aes(x=distance, y=N_rep4), color=3, span = 0.5)+geom_smooth(aes(x=distance, y=N_rep5), color=4, span = 0.5) + scale_x_log10(limits=c(10000,2500000)) + ylab("contact frequency (%)") + ylim(0,6.5)

# plotting distance decay for P-P only
interactions <- read.csv("Documents/Jul2022_MChIPC/pairs/interactions.txt", header = T, sep = "\t")
interactions <- interactions[,1:10]
interactions$distance <- log10(abs((interactions$bait_start+interactions$bait_end)/2 - (interactions$OE_start+interactions$OE_end)/2))
interactions <- filter(interactions, distance < 7)
baits <- read.csv("Documents/Jul2022_MChIPC/ChIPs/macs/binned_peaks.bed", header = F, sep = "\t")[,c(1:3)]
interactions_ranges <- GRanges(seqnames = interactions$OE_chr, ranges=IRanges(start=interactions$OE_start, end=interactions$OE_end, enh_id = seq(1,nrow(interactions))))
bait_ranges <- GRanges(seqnames = baits$V1, ranges=IRanges(start=baits$V2, end=baits$V3, enh_id = seq(1,nrow(baits))))
interactions <- interactions[unique(as.data.frame(findOverlaps(interactions_ranges,bait_ranges))[,1]),]
interactions <- interactions[,7:11]
interactions <- interactions %>% mutate(bin = cut(distance, breaks=seq(3.5,6.5,0.1), labels=FALSE))
interactions$sum <- interactions$N_rep0 + interactions$N_rep3 + interactions$N_rep4 + interactions$N_rep5
summary_PP <- interactions %>% group_by(bin) %>% summarise(N_rep0 = sum(N_rep0), N_rep3 = sum(N_rep3), N_rep4 = sum(N_rep4), N_rep5 = sum(N_rep5), N_sum=sum(sum))
summary_PP$N_rep0 <- 100*summary_PP$N_rep0/sum(summary_PP$N_rep0)
summary_PP$N_rep3 <- 100*summary_PP$N_rep3/sum(summary_PP$N_rep3)
summary_PP$N_rep4 <- 100*summary_PP$N_rep4/sum(summary_PP$N_rep4)
summary_PP$N_rep5 <- 100*summary_PP$N_rep5/sum(summary_PP$N_rep5)
summary_PP$N_sum <- 100*summary_PP$N_sum/sum(summary_PP$N_sum)
summary_PP$distance <- 10^(3.5+summary_PP$bin*0.1)
ggplot(summary_PP) + geom_smooth(aes(x=distance, y=N_rep0), color=1, span = 0.5) + geom_smooth(aes(x=distance, y=N_rep3), color=2, span = 0.5) +geom_smooth(aes(x=distance, y=N_rep4), color=3, span = 0.5)+geom_smooth(aes(x=distance, y=N_rep5), color=4, span = 0.5) + scale_x_log10(limits=c(10000,2500000))+ ylim(0,6.5) + ylab("contact frequency (%)")

ggplot() + geom_smooth(data= summary, aes(x=distance, y=N_sum), color=1, span = 0.5) + scale_x_log10(limits=c(10000,2500000)) + ylab("contact frequency (%)") +  geom_smooth(data= summary_PP, aes(x=distance, y=N_sum), color=2, span = 0.5) + scale_x_log10(limits=c(10000,2500000)) + ylab("contact frequency (%)")
