library('dplyr')
library('tidyr')
library('ggplot2')
library('fitdistrplus')
# loading contact sum files and coverage files
interactions_rep0 <- read.csv("Documents/Jul2022_MChIPC/interaction calling/interactions.sum.rep0.txt", header = F, sep = "\t")
interactions_rep3 <- read.csv("Documents/Jul2022_MChIPC/interaction calling/interactions.sum.rep3.txt", header = F, sep = "\t")
interactions_rep4 <- read.csv("Documents/Jul2022_MChIPC/interaction calling/interactions.sum.rep4.txt", header = F, sep = "\t")
interactions_rep5 <- read.csv("Documents/Jul2022_MChIPC/interaction calling/interactions.sum.rep5.txt", header = F, sep = "\t")

coverage_OE_rep0 <- read.csv("Documents/Jul2022_MChIPC/interaction calling/genomic_bins_coverage_rep0.bed", header = F, sep = "\t")
coverage_baits_rep0 <- read.csv("Documents/Jul2022_MChIPC/interaction calling/peaks_coverage_rep0.bed", header = F, sep = "\t")
coverage_baits_rep0 <- coverage_baits_rep0[,c(1,2,3,5)]
coverage_OE_rep3 <- read.csv("Documents/Jul2022_MChIPC/interaction calling/genomic_bins_coverage_rep3.bed", header = F, sep = "\t")
coverage_baits_rep3 <- read.csv("Documents/Jul2022_MChIPC/interaction calling/peaks_coverage_rep3.bed", header = F, sep = "\t")
coverage_baits_rep3 <- coverage_baits_rep3[,c(1,2,3,5)]
coverage_OE_rep4 <- read.csv("Documents/Jul2022_MChIPC/interaction calling/genomic_bins_coverage_rep4.bed", header = F, sep = "\t")
coverage_baits_rep4 <- read.csv("Documents/Jul2022_MChIPC/interaction calling/peaks_coverage_rep4.bed", header = F, sep = "\t")
coverage_baits_rep4 <- coverage_baits_rep4[,c(1,2,3,5)]
coverage_OE_rep5 <- read.csv("Documents/Jul2022_MChIPC/interaction calling/genomic_bins_coverage_rep5.bed", header = F, sep = "\t")
coverage_baits_rep5 <- read.csv("Documents/Jul2022_MChIPC/interaction calling/peaks_coverage_rep5.bed", header = F, sep = "\t")
coverage_baits_rep5 <- coverage_baits_rep5[,c(1,2,3,5)]

colnames(interactions_rep0) <- c("bait_chr","bait_start","bait_end","OE_chr","OE_start","OE_end", "N_rep0")
colnames(interactions_rep3) <- c("bait_chr","bait_start","bait_end","OE_chr","OE_start","OE_end", "N_rep3")
colnames(interactions_rep4) <- c("bait_chr","bait_start","bait_end","OE_chr","OE_start","OE_end", "N_rep4")
colnames(interactions_rep5) <- c("bait_chr","bait_start","bait_end","OE_chr","OE_start","OE_end", "N_rep5")

colnames(coverage_OE_rep0) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep0")
colnames(coverage_baits_rep0) <- c("bait_chr","bait_start","bait_end", "bait_cov_rep0")
colnames(coverage_OE_rep3) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep3")
colnames(coverage_baits_rep3) <- c("bait_chr","bait_start","bait_end", "bait_cov_rep3")
colnames(coverage_OE_rep4) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep4")
colnames(coverage_baits_rep4) <- c("bait_chr","bait_start","bait_end", "bait_cov_rep4")
colnames(coverage_OE_rep5) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep5")
colnames(coverage_baits_rep5) <- c("bait_chr","bait_start","bait_end", "bait_cov_rep5")

# combining all data in one df
interactions <- full_join(interactions_rep0, interactions_rep3)
interactions <- full_join(interactions, interactions_rep4)
interactions <- full_join(interactions, interactions_rep5)
interactions[is.na(interactions)] <- 0
interactions <- interactions %>% arrange(bait_chr,bait_start,OE_start)
interactions$N_sum <- interactions$N_rep0 + interactions$N_rep3 + interactions$N_rep4 + interactions$N_rep5

interactions <- left_join(interactions,coverage_OE_rep0)
interactions <- left_join(interactions,coverage_OE_rep3)
interactions <- left_join(interactions,coverage_OE_rep4)
interactions <- left_join(interactions,coverage_OE_rep5)
interactions <- left_join(interactions,coverage_baits_rep0)
interactions <- left_join(interactions,coverage_baits_rep3)
interactions <- left_join(interactions,coverage_baits_rep4)
interactions <- left_join(interactions,coverage_baits_rep5)

rm(list=ls()[grep("rep", ls())])

# assessing bait coverage and filtering poorly covered baits - in future this should be done after normalization, but now we have strange poorly ChIP covered baits
interactions <- left_join(interactions, group_by(interactions,bait_chr,bait_start,bait_end) %>% summarise(bait_coverage = sum(N_sum)))
interactions <- filter(interactions, bait_coverage>=1000)
# creating dist_rank column and filtering by distance
interactions$dist_rank <- 0
interactions$dist_rank[interactions$bait_start>interactions$OE_start] <- (interactions$bait_start[interactions$bait_start>interactions$OE_start] - interactions$OE_start[interactions$bait_start>interactions$OE_start])/250
interactions$dist_rank[interactions$bait_end<interactions$OE_start] <- (interactions$OE_end[interactions$bait_end<interactions$OE_start] - interactions$bait_end[interactions$bait_end<interactions$OE_start])/250
interactions <- filter(interactions, dist_rank>=20 & dist_rank<=10000)
# filterinf out baits with zero bait_cov (in future I have to figure out where do they come from)
#interactions <- filter(interactions, bait_cov != 0)
# ChIP-signal normalization (want the sum to remain the same)
interactions$N_rep0_norm <- interactions$N_rep0 / (interactions$bait_cov_rep0 + interactions$OE_cov_rep0)
interactions$N_rep0_norm <- interactions$N_rep0_norm*sum(interactions$N_rep0)/sum(interactions$N_rep0_norm)
interactions$N_rep3_norm <- interactions$N_rep3 / (interactions$bait_cov_rep3 + interactions$OE_cov_rep3)
interactions$N_rep3_norm <- interactions$N_rep3_norm*sum(interactions$N_rep3)/sum(interactions$N_rep3_norm)
interactions$N_rep4_norm <- interactions$N_rep4 / (interactions$bait_cov_rep4 + interactions$OE_cov_rep4)
interactions$N_rep4_norm <- interactions$N_rep4_norm*sum(interactions$N_rep4)/sum(interactions$N_rep4_norm)
interactions$N_rep5_norm <- interactions$N_rep5 / (interactions$bait_cov_rep5 + interactions$OE_cov_rep5)
interactions$N_rep5_norm <- interactions$N_rep5_norm*sum(interactions$N_rep5)/sum(interactions$N_rep5_norm)
# summarising replicates
interactions$N_sum_norm <- interactions$N_rep0_norm + interactions$N_rep3_norm + interactions$N_rep4_norm + interactions$N_rep5_norm
#bait coverage after normalization
bait_coverage <- group_by(interactions,bait_chr,bait_start,bait_end) %>% summarise(bait_coverage = sum(N_sum))
bait_coverage <- left_join(bait_coverage, group_by(interactions,bait_chr,bait_start,bait_end) %>% summarise(bait_coverage_norm= sum(N_sum_norm)))
interactions <- left_join(interactions, bait_coverage[,c(1,2,3,5)])
interactions$N_sum_double_norm <- interactions$N_sum_norm*median(interactions$bait_coverage_norm) / interactions$bait_coverage_norm
bait_coverage <- left_join(bait_coverage, group_by(interactions,bait_chr,bait_start,bait_end) %>% summarise(bait_coverage_double_norm= sum(N_sum_double_norm)))
ggplot(bait_coverage) + geom_boxplot(aes(x="coverage",y=bait_coverage)) + geom_boxplot(aes(x="norm_coverage",y=bait_coverage_norm)) + geom_boxplot(aes(x="norm_coverage_double",y=bait_coverage_double_norm)) + ylim(0,10000)
#plot dist function
dist_func <- group_by(interactions, dist_rank) %>% summarise(int = sum(N_sum_double_norm))
ggplot(dist_func) + geom_point(aes(x=dist_rank, y=int)) + xlim(20,1000)
#building models, calculating p-values (this could be definitely changed)
interactions$p_val <- NA
interactions$p_val <- as.numeric(interactions$p_val)
for (i in seq(20,4000)){
  print(paste0("calculating p_val for ",i," dist_rank"))
  interactions_rank <- filter(interactions, dist_rank==i)
  x <- fitdistr(interactions_rank$N_sum_double_norm[interactions_rank$N_sum_double_norm < quantile(interactions_rank$N_sum_double_norm, 0.95)], "weibull")
  interactions_rank$p_val <- 1-pweibull(interactions_rank$N_sum_double_norm, shape = x$estimate["shape"], scale=x$estimate["scale"])
  interactions <- interactions %>% rows_update(interactions_rank, by=c("bait_chr","bait_start","bait_end","OE_chr","OE_start","OE_end"))}
# writing full interaction statistics file
write.table(interactions, file="Documents/Jul2022_MChIPC/interaction calling/interactions.txt", col.names = T, row.names = F, sep = "\t", quote = F)
# true loops p_val<0.01 & 6+ ligation products | interactions with 6+ total ligation products (absolutely arbitrary) further than 1Mbp from viewpoint
loops <- filter(interactions, (p_val<=0.01 & N_sum > 5) | (dist_rank > 4000 & N_sum > 5))
write.table(loops, file="Documents/Jul2022_MChIPC/interaction calling/loops.min.6.log10.2.bedpe", col.names = T, row.names = F, sep = "\t", quote = F)

