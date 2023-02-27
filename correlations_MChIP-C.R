library("dplyr")
library("tidyr")
library("ggplot2")
library("GGally")
library("GenomicRanges")
GGscatterhex <- function(data, mapping,...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  df <- data.frame(x = x, y = y)
  pp <- ggplot(df, aes(x=x, y=y)) + geom_hex(bins=100) + xlim(0,4) + ylim(0,4) + scale_fill_continuous(limits = c(0, 100), oob = scales::squish)
  return(pp)}
# loading interactions and MultibigwigSummary file
interactions <- read.csv("~/Documents/Jul2022_MChIPC/interaction calling/interactions.txt", header = T, sep = "\t")
interactions <- filter(interactions, dist_rank <= 4000)
multi <- read.csv("~/Documents/Jul2022_MChIPC/fig.S1/replicates_coverage.txt", header = T, sep = "\t")
multi <- filter(multi, (X.MChIPC_rep0.bw.>0&X.MChIPC_rep3.bw.>0&X.MChIPC_rep4.bw.>0&X.MChIPC_rep5.bw.>0))
# loading viewpoints and excluding them from analysis
baits <- read.csv("~/Documents/Jul2022_MChIPC/ChIPs/macs/binned_peaks.bed", header = F, sep = "\t")
bait_ranges <- GRanges(seqnames = baits$V1, ranges=IRanges(start=baits$V2, end=baits$V3, enh_id = seq(1,nrow(baits))))
multi_ranges <- GRanges(seqnames = multi$X..chr., ranges=IRanges(start=multi$X.start., end=multi$X.end., enh_id = seq(1,nrow(multi))))
interactions_ranges <- GRanges(seqnames = interactions$OE_chr, ranges=IRanges(start=interactions$OE_start, end=interactions$OE_end, enh_id = seq(1,nrow(interactions))))
multi <- multi[-(as.data.frame(findOverlaps(multi_ranges,bait_ranges))[,1]),]
interactions <- interactions[-(as.data.frame(findOverlaps(interactions_ranges,bait_ranges))[,1]),]
# plotting graphs
plot_multi <- ggpairs(log10(multi[,4:7]), lower = list(continuous=wrap(GGscatterhex)), upper = list(continuous = wrap("cor", method = "pearson")))
plot_interactions <- ggpairs(log10(interactions[,10:7]), upper = list() , lower = list(continuous=wrap(GGscatterhex)))
# calculating correlations
# merged
print(paste0("r_merged(rep0,rep3)=",cor(multi$X.MChIPC_rep0.bw.,multi$X.MChIPC_rep3.bw.)))
print(paste0("r_merged(rep0,rep4)=",cor(multi$X.MChIPC_rep0.bw.,multi$X.MChIPC_rep4.bw.)))
print(paste0("r_merged(rep0,rep5)=",cor(multi$X.MChIPC_rep0.bw.,multi$X.MChIPC_rep5.bw.)))
print(paste0("r_merged(rep3,rep4)=",cor(multi$X.MChIPC_rep3.bw.,multi$X.MChIPC_rep4.bw.)))
print(paste0("r_merged(rep3,rep5)=",cor(multi$X.MChIPC_rep3.bw.,multi$X.MChIPC_rep5.bw.)))
print(paste0("r_merged(rep4,rep5)=",cor(multi$X.MChIPC_rep4.bw.,multi$X.MChIPC_rep5.bw.)))
# separate bins
print(paste0("r(rep0,rep3)=",cor(interactions$N_rep0, interactions$N_rep3)))
print(paste0("r(rep0,rep4)=",cor(interactions$N_rep0, interactions$N_rep4)))
print(paste0("r(rep0,rep5)=",cor(interactions$N_rep0, interactions$N_rep5)))
print(paste0("r(rep3,rep4)=",cor(interactions$N_rep3, interactions$N_rep4)))
print(paste0("r(rep3,rep5)=",cor(interactions$N_rep3, interactions$N_rep5)))
print(paste0("r(rep4,rep5)=",cor(interactions$N_rep4, interactions$N_rep5)))
