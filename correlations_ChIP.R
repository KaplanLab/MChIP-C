library('dplyr')
library('tidyr')
library("ggplot2")
library("GGally")
GGscatterhex <- function(data, mapping,...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  df <- data.frame(x = x, y = y)
  pp <- ggplot(df, aes(x=x, y=y)) + geom_hex(bins=100) + xlim(-1,3) + ylim(-1,3) + scale_fill_continuous(limits = c(0, 50), oob = scales::squish)
  return(pp)}
# creating filterd DNase peaks
DNase <- read.csv("Documents/Nov2021_MCC/Chips/peaks/ENCFF621ZJY.bed.gz", header = F, sep = "\t")
DNase <- filter(DNase, V7 > 200)
write.table(DNase[,1:3], "Documents/Jul2022_MChIPC/fig.S1/DNase.filtered.bed", quote = F, col.names = F, row.names = F, sep = "\t")
# estimating ChIP coverage in the selected DHSs with deeptools
# multiBigwigSummary BED-file -b ../ChIPs/ChIP_rep0.bw ../ChIPs/ChIP_rep3.bw ../ChIPs/ChIP_rep4.bw ../ChIPs/ChIP_rep5.bw ../ChIPs/ChIP_H3K4me3.bw --outRawCounts H3K4me3_coverage.txt --BED DNase.filtered.bed -o results.npz
coverage <- read.csv("Documents/Jul2022_MChIPC/fig.S1/H3K4me3_coverage.txt", header = T, sep = "\t")
plot <- ggpairs(log10(coverage[,4:8]), lower = list(continuous=wrap(GGscatterhex)), upper = list())
# calculating correlations
# merged
print(paste0("r_merged(rep0,rep3)=",cor(coverage$X.ChIP_rep0.bw.,coverage$X.ChIP_rep3.bw.)))
print(paste0("r_merged(rep0,rep4)=",cor(coverage$X.ChIP_rep0.bw.,coverage$X.ChIP_rep4.bw.)))
print(paste0("r_merged(rep0,rep5)=",cor(coverage$X.ChIP_rep0.bw.,coverage$X.ChIP_rep5.bw.)))
print(paste0("r_merged(rep0,ChIP)=",cor(coverage$X.ChIP_rep0.bw.,coverage$X.ChIP_H3K4me3.bw.)))
print(paste0("r_merged(rep3,rep4)=",cor(coverage$X.ChIP_rep3.bw.,coverage$X.ChIP_rep4.bw.)))
print(paste0("r_merged(rep3,rep5)=",cor(coverage$X.ChIP_rep3.bw.,coverage$X.ChIP_rep5.bw.)))
print(paste0("r_merged(rep3,ChIP)=",cor(coverage$X.ChIP_rep3.bw.,coverage$X.ChIP_H3K4me3.bw.)))
print(paste0("r_merged(rep4,rep5)=",cor(coverage$X.ChIP_rep4.bw.,coverage$X.ChIP_rep5.bw.)))
print(paste0("r_merged(rep4,ChIP)=",cor(coverage$X.ChIP_rep4.bw.,coverage$X.ChIP_H3K4me3.bw.)))
print(paste0("r_merged(rep5,ChIP)=",cor(coverage$X.ChIP_rep5.bw.,coverage$X.ChIP_H3K4me3.bw.)))
