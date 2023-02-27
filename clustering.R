library(dplyr)
library(tidyr)
library(dendextend)
library(dendroextras)
library(GenomicRanges)
library(ggplot2)
# first creating a df with DNase sites overlapping PIRs and protein binding sites
read_bed <- function(x){
  bed <- read.csv(paste0(folder,x[1,2]), header = F, sep = "\t")
  if(sum(bed$V8) < 0) {bed <- filter(bed, V7>25)} else {bed <- filter(bed, V8>20)}
}
# loading files with ChIP peaks and PIRs
folder <- "Documents/Nov2021_MCC/Chips/peaks/"
list <- read.csv("Documents/Nov2021_MCC/Chips/list_of_beds_and_bigwigs.txt", header = T, sep = "\t")[,1:2]
list <- split(list, list$target)
beds <- lapply(list, read_bed)
peak_ranges <- lapply(beds, function(bed){GRanges(seqnames = bed$V1, ranges=IRanges(start=bed$V2, end=bed$V3, enh_id = seq(1,nrow(bed))))})
PIRs <- read.csv("Documents/Jul2022_MChIPC/fig.3/PIRs.bed", header = F, sep = "\t")
PIR_ranges <- GRanges(seqnames = PIRs$V1, ranges=IRanges(start=PIRs$V2, end=PIRs$V3, enh_id = seq(1,nrow(PIRs))))
# loading DNase peaks and selecting DNase peaks overlapping PIRs:
DNase <- read.csv("Documents/Nov2021_MCC/Chips/peaks/ENCFF621ZJY.bed.gz", header = F, sep = "\t")
DNase <-filter(DNase, V7>200)
DNase_ranges <- GRanges(seqnames = DNase$V1, ranges=IRanges(start=DNase$V2, end=DNase$V3, enh_id = seq(1,nrow(DNase))))
DNase_in_PIRs <- DNase[unique(as.data.frame(findOverlaps(DNase_ranges,PIR_ranges))[,1]), 1:3]
DNase_in_PIRs_ranges <- GRanges(seqnames = DNase_in_PIRs$V1, ranges=IRanges(start=DNase_in_PIRs$V2, end=DNase_in_PIRs$V3, enh_id = seq(1,nrow(DNase_in_PIRs))))
# finding overlaps with peaks:
overlaps <- sapply(peak_ranges, function(range){unique(as.data.frame(findOverlaps(DNase_in_PIRs_ranges, range))[,1])})
for (i in seq(length(overlaps))){
  DNase_in_PIRs[names(overlaps)[i]] <- 0
  DNase_in_PIRs[names(overlaps)[i]][unlist(overlaps[i]),] <- 1}
colnames(DNase_in_PIRs)[1:3] <- c("DNase_chr","DNase_start","DNase_end")
write.table(DNase_in_PIRs, file="Documents/Jul2022_MChIPC/fig.3/hier_clustering/DNase_overlaps.txt", col.names = T, row.names = F, sep = "\t", quote = F)
rm(list = ls())
# if you are repeating you can start here
DNase <- read.csv("Documents/Jul2022_MChIPC/fig.3/hier_clustering/DNase_overlaps.txt", header = T, sep = "\t")
enrichment <- read.csv("Documents/Jul2022_MChIPC/fig.3/enrichment_summary.txt", header = T, sep = "\t")
enrichment <- filter(enrichment, enrichment > 2)
for (i in seq(ncol(DNase),4)){
  if (sum(DNase[,i]) < 1000) {DNase[,i] <- NULL}
  if (!(colnames(DNase)[i] %in% enrichment$ChIP)) {DNase[,i] <- NULL}}
distance <- dist(DNase[,4:167], method="binary")
cl <- hclust(distance, method="ward.D")
cl$height <- round(cl$height, 6)
dend1 <- as.dendrogram(cl)
dend2 <- color_branches(dend1, groupLabels = T, k=5)
labels(dend2) <- NULL
# plotting a dendrogram:
plot(dend2)
# plotting hetamap
DNase$cluster <- dendroextras::slice(dend1, k=5)[as.character(1:18837)]
DNase <- DNase[labels(dend1),]
order <- as.character(read.csv("Documents/Jul2022_MChIPC/fig.S3/matrix_indeces.txt", header = F)$V1)
heatmap_1 <- heatmap(as.matrix(DNase[1:2500,order]),  Colv = NA,  Rowv = NA, scale = "none", labRow = FALSE, labCol = FALSE, keep.dendro=F)
heatmap_2 <- heatmap(as.matrix(DNase[2501:5000,order]),  Colv = NA,  Rowv = NA, scale = "none", labRow = FALSE, labCol = FALSE, keep.dendro=F)
heatmap_3 <- heatmap(as.matrix(DNase[5001:7500,order]),  Colv = NA, Rowv = NA, scale = "none", labRow = FALSE, labCol = FALSE, keep.dendro=F)
heatmap_4 <- heatmap(as.matrix(DNase[7501:10000,order]),  Colv = NA, Rowv = NA, scale = "none", labRow = FALSE, labCol = FALSE, keep.dendro=F)
heatmap_5 <- heatmap(as.matrix(DNase[10001:12500,order]),  Colv = NA, Rowv = NA, scale = "none", labRow = FALSE, labCol = FALSE, keep.dendro=F)
heatmap_6 <- heatmap(as.matrix(DNase[12501:15000,order]),  Colv = NA, Rowv = NA, scale = "none", labRow = FALSE, labCol = FALSE, keep.dendro=F)
heatmap_7 <- heatmap(as.matrix(DNase[15001:17500,order]),  Colv = NA, Rowv = NA, scale = "none", labRow = FALSE, labCol = FALSE, keep.dendro=F)
heatmap_8 <- heatmap(as.matrix(DNase[17501:18837,order]),  Colv = NA, Rowv = NA, scale = "none", labRow = FALSE, labCol = FALSE, keep.dendro=F)

full_heatmap <- heatmap(as.matrix(DNase[,4:167]), Rowv = NA, scale = "none", labRow = FALSE, labCol = FALSE, keep.dendro=F)

# analyzing feature enrichment in clusters
feature_sum <- DNase[,4:168] %>% group_by(cluster) %>% summarise_each(mean)
feature_sum <- as.data.frame(t(feature_sum[,2:165]))
features <- c()
features <- append(features, rownames(arrange(feature_sum, desc(V1)) %>% filter(V1>0.5))[1:20])
features <- append(features, rownames(arrange(feature_sum, desc(V2)) %>% filter(V2>0.5))[1:20])
features <- append(features, rownames(arrange(feature_sum, desc(V3))  %>% filter(V3>0.5))[1:20])
features <- append(features, rownames(arrange(feature_sum, desc(V4))  %>% filter(V4>0.5))[1:20])
features <- append(features, rownames(arrange(feature_sum, desc(V5))  %>% filter(V5>0.5))[1:20])
features <- unique(features)
features <- features[!(is.na(features))]
# plotting heatmap for selected features
heatmap(as.matrix(DNase[features]), Colv = NA, Rowv = NA, scale = "none", labRow = FALSE, labCol = colnames(DNase[features]))
# analysing enhancers
# counting overlaps with confirmed enhancers
fulco_enh <- read.csv("Documents/Jul2022_MChIPC/fig.4/fulco_full_data.txt", header = T, sep = "\t")
fulco_enh <- filter(fulco_enh, Significant==1 & Fraction.change.in.gene.expr < 0)
gasperini_enh <- read.csv("Documents/Jul2022_MChIPC/fig.4/Gasperini_EP_pairs.txt", header = T, sep = "\t")
gasperini_enh <- filter(gasperini_enh, Significant)
colnames(gasperini_enh)[1:3] <- colnames(fulco_enh)[5:7]
enhancers <- bind_rows(gasperini_enh[,1:3], fulco_enh[5:7])
enhancer_ranges <- GRanges(seqnames = enhancers$chr, ranges=IRanges(start=enhancers$start, end=enhancers$end, enh_id = seq(1,nrow(enhancers))))
DNase_ranges <- GRanges(seqnames = DNase$DNase_chr, ranges=IRanges(start=DNase$DNase_start, end=DNase$DNase_end, enh_id = seq(1,nrow(DNase))))
DNase$number <- 1:nrow(DNase)
DNase$enhancer_overlap <- DNase$number %in% unique(as.data.frame(findOverlaps(DNase_ranges,enhancer_ranges))[,1])
ggplot(filter(DNase, enhancer_overlap)) + geom_jitter(aes(number, enhancer_overlap), size=0.25, height = 0.1)+geom_vline(xintercept=nrow(DNase))
# testing enrichment
binom.test(sum(filter(DNase, cluster==1)$enhancer_overlap), sum(DNase$enhancer_overlap), nrow(filter(DNase, cluster==1))/nrow(DNase), alternative = "g")
binom.test(sum(filter(DNase, cluster==2)$enhancer_overlap), sum(DNase$enhancer_overlap), nrow(filter(DNase, cluster==2))/nrow(DNase), alternative = "g")
binom.test(sum(filter(DNase, cluster==3)$enhancer_overlap), sum(DNase$enhancer_overlap), nrow(filter(DNase, cluster==3))/nrow(DNase), alternative = "g")
binom.test(sum(filter(DNase, cluster==4)$enhancer_overlap), sum(DNase$enhancer_overlap), nrow(filter(DNase, cluster==4))/nrow(DNase), alternative = "g")
binom.test(sum(filter(DNase, cluster==5)$enhancer_overlap), sum(DNase$enhancer_overlap), nrow(filter(DNase, cluster==5))/nrow(DNase), alternative = "g")
# writing bed files for individual clusters
write.table(filter(DNase, cluster==1)[,1:3], file="Documents/Jul2022_MChIPC/fig.3/hier_clustering/cluster1.bed", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(filter(DNase, cluster==2)[,1:3], file="Documents/Jul2022_MChIPC/fig.3/hier_clustering/cluster2.bed", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(filter(DNase, cluster==3)[,1:3], file="Documents/Jul2022_MChIPC/fig.3/hier_clustering/cluster3.bed", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(filter(DNase, cluster==4)[,1:3], file="Documents/Jul2022_MChIPC/fig.3/hier_clustering/cluster4.bed", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(filter(DNase, cluster==5)[,1:3], file="Documents/Jul2022_MChIPC/fig.3/hier_clustering/cluster5.bed", col.names = F, row.names = F, sep = "\t", quote = F)

# number of structural and regulatory loops and their average length
loops <-read.csv("Documents/Jul2022_MChIPC/interaction calling/loops.min.6.log10.2.bedpe", header = F, sep = "\t")
loops$distance <- abs((loops$V2+loops$V3)/2 - (loops$V5+loops$V6)/2)
OE_ranges <-  GRanges(seqnames = loops$V4, ranges=IRanges(start=loops$V5, end=loops$V6, enh_id = seq(1,nrow(loops))))
structural_PIRs <- filter(DNase, cluster==1 | cluster ==2)[,1:3]
structural_ranges <- GRanges(seqnames = structural_PIRs$DNase_chr, ranges=IRanges(start=structural_PIRs$DNase_start, end=structural_PIRs$DNase_end, enh_id = seq(1,nrow(structural_PIRs))))
regulatory_PIRs <- filter(DNase, cluster==3)[,1:3]
regulatory_ranges <- GRanges(seqnames = regulatory_PIRs$DNase_chr, ranges=IRanges(start=regulatory_PIRs$DNase_start, end=regulatory_PIRs$DNase_end, enh_id = seq(1,nrow(regulatory_PIRs))))
structural_loops <- loops[unique(as.data.frame(findOverlaps(OE_ranges, structural_ranges))[,1]),]
regulatory_loops <- loops[unique(as.data.frame(findOverlaps(OE_ranges, regulatory_ranges))[,1]),]
print(median(loops$distance))
print(median(structural_loops$distance))
print(median(regulatory_loops$distance))

# plotting loops length distribution
cl1_PIRs <- filter(DNase, cluster==1)[,1:3]
cl2_PIRs <- filter(DNase, cluster==2)[,1:3]
cl3_PIRs <- filter(DNase, cluster==3)[,1:3]
cl4_PIRs <- filter(DNase, cluster==4)[,1:3]
cl5_PIRs <- filter(DNase, cluster==5)[,1:3]
cl1_ranges <- GRanges(seqnames = cl1_PIRs$DNase_chr, ranges=IRanges(start=cl1_PIRs$DNase_start, end=cl1_PIRs$DNase_end, enh_id = seq(1,nrow(cl1_PIRs))))
cl2_ranges <- GRanges(seqnames = cl2_PIRs$DNase_chr, ranges=IRanges(start=cl2_PIRs$DNase_start, end=cl2_PIRs$DNase_end, enh_id = seq(1,nrow(cl2_PIRs))))
cl3_ranges <- GRanges(seqnames = cl3_PIRs$DNase_chr, ranges=IRanges(start=cl3_PIRs$DNase_start, end=cl3_PIRs$DNase_end, enh_id = seq(1,nrow(cl3_PIRs))))
cl4_ranges <- GRanges(seqnames = cl4_PIRs$DNase_chr, ranges=IRanges(start=cl4_PIRs$DNase_start, end=cl4_PIRs$DNase_end, enh_id = seq(1,nrow(cl4_PIRs))))
cl5_ranges <- GRanges(seqnames = cl5_PIRs$DNase_chr, ranges=IRanges(start=cl5_PIRs$DNase_start, end=cl5_PIRs$DNase_end, enh_id = seq(1,nrow(cl5_PIRs))))
cl1_loops <- loops[unique(as.data.frame(findOverlaps(OE_ranges, cl1_ranges))[,1]),]
cl2_loops <- loops[unique(as.data.frame(findOverlaps(OE_ranges, cl2_ranges))[,1]),]
cl3_loops <- loops[unique(as.data.frame(findOverlaps(OE_ranges, cl3_ranges))[,1]),]
cl4_loops <- loops[unique(as.data.frame(findOverlaps(OE_ranges, cl4_ranges))[,1]),]
cl5_loops <- loops[unique(as.data.frame(findOverlaps(OE_ranges, cl5_ranges))[,1]),]
cl1_loops$cl <- "cl1"
cl2_loops$cl <- "cl2"
cl3_loops$cl <- "cl3"
cl4_loops$cl <- "cl4"
cl5_loops$cl <- "cl5"
loops <- bind_rows(cl1_loops,cl2_loops,cl3_loops,cl4_loops,cl5_loops)
ggplot(loops) + geom_boxplot(aes(cl,log10(distance)))

# clusters and chromatin colors
HMM <- read.csv("Documents/Jul2022_MChIPC/ABC-Enhancer-Gene-Prediction/reference/wgEncodeBroadHmmK562HMM.bed.gz", header = F, sep = "\t")
HMM$subjectHits <- seq(1, nrow(HMM))
for (i in 1:nrow(HMM)){
  HMM$V9[i] <- rgb(as.numeric(strsplit(HMM$V9[i], split=",")[[1]])[1], as.numeric(strsplit(HMM$V9[i], split=",")[[1]])[2], as.numeric(strsplit(HMM$V9[i], split=",")[[1]])[3], maxColorValue = 255)}
HMM_ranges <- GRanges(seqnames = HMM$V1, ranges=IRanges(start=HMM$V2, end=HMM$V3, enh_id = seq(1,nrow(HMM))))
summary <- data.frame()
add_to_summary<- function(x){
  cluster <- read.csv(paste0("Documents/Jul2022_MChIPC/fig.3/hier_clustering/",x,".bed"), header = F, sep = "\t")
  cluster$queryHits <- seq(1,nrow(cluster))
  cluster_ranges <- GRanges(seqnames = cluster$V1, ranges=IRanges(start=cluster$V2, end=cluster$V3, enh_id = seq(1,nrow(cluster))))
  cluster <- left_join(cluster, as.data.frame(findOverlaps(cluster_ranges, HMM_ranges, minoverlap = 76))[!(duplicated(as.data.frame(findOverlaps(cluster_ranges, HMM_ranges, minoverlap = 76))[,1])),])
  cluster <- left_join(cluster, HMM[,c(4,9,10)])
  cluster <- cluster[!(is.na(cluster$V9)),]
  summary <- cluster %>% group_by(V4) %>% summarise(N=length(V4), color=unique(V9), cluster=x)
  return(summary)}
for (i in c("cluster1","cluster2","cluster3","cluster4","cluster5")){
  summary <- rbind(summary, add_to_summary(i))}
summary$V5 <- as.numeric(unlist(lapply(strsplit(summary$V4, split="_"), `[[`, 1)))
summary <- arrange(summary, V5)
ggplot(summary) + geom_col(aes(x=cluster, y=N, fill=factor(V4, levels=unique(V4)))) + scale_fill_manual(values=summary$color[!duplicated(summary$V4)])


