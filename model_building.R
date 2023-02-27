library('dplyr')
library('tidyr')
library('GenomicRanges')
library('ranger')
library('caret')
library("PRROC")
library('reshape2')
library("ggplot2")
baits <- read.csv("Documents/Jul2022_MChIPC/proper_baits.txt", header = T, sep = "\t")
DNase <- read.csv("Documents/Nov2021_MCC/Chips/peaks/ENCFF621ZJY.bed.gz", header = F, sep = "\t")
DNase <- filter(DNase, V7 > 200)
H3K4me3_mask <- read.csv("Documents/Jul2022_MChIPC/ChIPs/macs/H3K4me3.mask.bed", header = F, sep = "\t")
H3K4me3_ranges <- GRanges(seqnames = H3K4me3_mask$V1, ranges=IRanges(start=H3K4me3_mask$V2-5000, end=H3K4me3_mask$V3+5000, enh_id = seq(1,nrow(H3K4me3_mask))))
DNase_ranges <- GRanges(seqnames = DNase$V1, ranges=IRanges(start=DNase$V2, end=DNase$V3, enh_id = seq(1,nrow(DNase))))
DNase <- DNase[-as.data.frame(findOverlaps(DNase_ranges,H3K4me3_ranges))[,1],]
write.table(DNase[,1:3], "Documents/Jul2022_MChIPC/fig.3/random forest/DNase_for_forest.bed", quote = F, col.names = F, row.names = F, sep = "\t")
bait_ranges <- GRanges(seqnames = baits$bait_chr, ranges=IRanges(start=baits$bait_start, end=baits$bait_end, enh_id = seq(1,nrow(baits))))
bait_left_ranges <- GRanges(seqnames = baits$bait_chr, ranges=IRanges(start=baits$bait_start-250000, end=baits$bait_start-5000, enh_id = seq(1,nrow(baits))))
bait_right_ranges <- GRanges(seqnames = baits$bait_chr, ranges=IRanges(start=baits$bait_end+5000, end=baits$bait_end+250000, enh_id = seq(1,nrow(baits))))
DNase_ranges <- GRanges(seqnames = DNase$V1, ranges=IRanges(start=DNase$V2, end=DNase$V3, enh_id = seq(1,nrow(DNase))))
left_overlaps <- countOverlaps(DNase_ranges, bait_left_ranges)
right_overlaps <- countOverlaps(DNase_ranges, bait_right_ranges)
DNase$V4 <- left_overlaps+right_overlaps
DNase <- filter(DNase, V4>0)
DNase_ranges <- GRanges(seqnames = DNase$V1, ranges=IRanges(start=DNase$V2, end=DNase$V3, enh_id = seq(1,nrow(DNase))))
contacts <- read.csv("Documents/Jul2022_MChIPC/pairs/interactions.txt", header = T, sep = "\t")
contacts <- contacts[,c(1:6,11,21,29)]
contacts <- filter(contacts, dist_rank <= 1000)
loops <- data_frame()
for (i in seq(1,nrow(DNase))){
  site <- DNase[i,]
  contact_site <- filter(contacts, bait_chr==site$V1 & OE_start < site$V3 & OE_end > site$V2)
  if (nrow(contact_site)>0){
     contact_site$site_id <- i  
     loops <- rbind(loops,contact_site)}
  else {}
  print(paste0("site ",i," is finished"))
}
loops <- loops %>% group_by(bait_chr,bait_start,bait_end,site_id) %>% summarise(N_sum_double_norm = max(N_sum_double_norm)) %>% arrange(site_id)
DNase$site_id <- seq(1,nrow(DNase))
loops <- left_join(loops, DNase)
loops$distance <- (loops$V2+loops$V3)/2 - (loops$bait_start+loops$bait_end)/2
loops <- loops[,c(1,2,3,6,7,8,5,10)]
# looking for sites overlapping CTCF peaks:
CTCF_peaks <- read.csv("Documents/Nov2021_MCC/Chips/peaks/ENCFF002CEL.bed.gz", header = F, sep = "\t")
CTCF_ranges <- GRanges(seqnames = CTCF_peaks$V1, ranges=IRanges(start=CTCF_peaks$V2, end=CTCF_peaks$V3, enh_id = seq(1,nrow(CTCF_peaks))))
baits_wCTCF <- baits[(unique(as.data.frame(findOverlaps(bait_ranges, CTCF_ranges))[,1])),]
DNase_wCTCF <- DNase[(unique(as.data.frame(findOverlaps(DNase_ranges, CTCF_ranges))[,1])),]
write.table(baits_wCTCF[,1:3], "Documents/Jul2022_MChIPC/fig.3/random forest/baits_wCTCF.txt", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(DNase_wCTCF[,1:3], "Documents/Jul2022_MChIPC/fig.3/random forest/DNase_wCTCF.txt", quote = F, col.names = F, row.names = F, sep = "\t")
# using Homer to find motifs and their orientation in CTCF-overlapping sites (and deleting headers of the output manually)
# annotatePeaks.pl DNase_wCTCF.txt hg19 -nogene -noann -m ../../../../Downloads/homer/motifs/ctcf.motif > motifs_in_DNase
# annotatePeaks.pl baits_wCTCF.txt hg19 -nogene -noann -m ../../../../Downloads/homer/motifs/ctcf.motif > motifs_in_baits
baits_wCTCF <- read.csv("Documents/Jul2022_MChIPC/fig.3/random forest/motifs_in_baits", header = F, sep = "\t")
baits_wCTCF$baits_CTCF_minus <- grepl("\\-",baits_wCTCF$V10)
baits_wCTCF$baits_CTCF_plus <- grepl("\\+",baits_wCTCF$V10)
baits_wCTCF <- baits_wCTCF[,c(2:4,11,12)]
DNase_wCTCF <- read.csv("Documents/Jul2022_MChIPC/fig.3/random forest/motifs_in_DNase", header = F, sep = "\t")
DNase_wCTCF$DNase_CTCF_minus <- grepl("\\-",DNase_wCTCF$V10)
DNase_wCTCF$DNase_CTCF_plus <- grepl("\\+",DNase_wCTCF$V10)
DNase_wCTCF <- DNase_wCTCF[,c(2:4,11,12)]
# joining data
loops[is.na(loops)] <- 0
colnames(baits_wCTCF)[1:3]<- colnames(loops)[1:3]
colnames(DNase_wCTCF)[1:3]<- colnames(loops)[4:6]
loops <- left_join(loops, baits_wCTCF[,c(1,3:5)])
loops <- left_join(loops, DNase_wCTCF[,c(1,3:5)])
loops[is.na(loops)] <- FALSE
# adding CTCF signal on baits and joining OE ChIP-seq data
# coverage for baits is created with deeptools:
# multiBigwigSummary BED-file -b ../Nov2021_MCC/Chips/ENCFF000BWF\ -\ CTCF.bigWig --outRawCounts fig.3/random\ forest/baits_CTCF_coverage.txt --BED ChIPs/macs/binned_peaks.bed -o results.npz
bait_CTCF_coverage <- read.csv("Documents/Jul2022_MChIPC/fig.3/random forest/baits_CTCF_coverage.txt", header = T, sep = "\t")
colnames(bait_CTCF_coverage)[1:3] <- colnames(loops)[1:3]
colnames(bait_CTCF_coverage)[4] <-"bait_CTCF_coverage"
loops <- left_join(loops, bait_CTCF_coverage)
loops$bait_CTCF_coverage[is.na(loops$bait_CTCF_coverage)] <-0
# loading data for ChIP-seq signal in all DNase sites
# created coverage with deeptools:
# multiBigwigSummary BED-file -b *.bigWig -o results.npz --outRawCounts ChIP_signal.txt -p 7 --BED DNase_for_forest.bed 
ChIP_data <- read.csv("Documents/Jul2022_MChIPC/fig.3/random forest/ChIP_signal.txt", header = T, sep = "\t")
ChIP_data[is.na(ChIP_data)] <- 0
colnames(ChIP_data)[1:3] <- colnames(loops)[4:6]
loops <- left_join(loops, ChIP_data)
colnames(loops)[4:7] <- c("DNase_chr","DNase_start","DNase_end", "interaction_score")
loops$relative_posit <- loops$distance > 0
loops$distance <- abs(loops$distance)
write.table(loops, "Documents/Jul2022_MChIPC/fig.3/random forest/full_data_for_forests.txt", quote = F, col.names = T, row.names = F, sep = "\t")
rm(list=ls())
# starting from here, if tha data is ready 
loops <- read.csv("Documents/Jul2022_MChIPC/full_data_for_forests.txt", header = T, sep = "\t")
loops1 <- filter(loops, !(interaction_score %in% sort(loops$interaction_score, decreasing = TRUE)[c(1,2)]))

set.seed(42)

# building distance-only model
data_for_model <- loops1[,c(7:8)]
ind <- sample(3, nrow(data_for_model), replace = TRUE, prob = c(0.33, 0.33, 0.33))
test_1 <- data_for_model[ind==1,]
test_2 <- data_for_model[ind==2,]
test_3 <- data_for_model[ind==3,]
model_1 <- ranger(interaction_score~., data=bind_rows(test_2, test_3), num.trees=100, mtry = 1)
model_2 <- ranger(interaction_score~., data=bind_rows(test_1, test_3), num.trees=100, mtry = 1)
model_3 <- ranger(interaction_score~., data=bind_rows(test_1, test_2), num.trees=100, mtry = 1)
rsq_distance <- c((cor(predict(model_1, test_1)$predictions, test_1$interaction_score))^2,(cor(predict(model_2, test_2)$predictions, test_2$interaction_score))^2,(cor(predict(model_3, test_3)$predictions, test_3$interaction_score))^2)

# building initial model
data_for_model <- loops1[,c(7:13,16,285)]
ind <- sample(3, nrow(data_for_model), replace = TRUE, prob = c(0.33, 0.33, 0.33))
test_1 <- data_for_model[ind==1,]
test_2 <- data_for_model[ind==2,]
test_3 <- data_for_model[ind==3,]
model_1 <- ranger(interaction_score~., data=bind_rows(test_2, test_3), num.trees=100, mtry = 3)
model_2 <- ranger(interaction_score~., data=bind_rows(test_1, test_3), num.trees=100, mtry = 3)
model_3 <- ranger(interaction_score~., data=bind_rows(test_1, test_2), num.trees=100, mtry = 3)
rsq_initial <- c((cor(predict(model_1, test_1)$predictions, test_1$interaction_score))^2,(cor(predict(model_2, test_2)$predictions, test_2$interaction_score))^2,(cor(predict(model_3, test_3)$predictions, test_3$interaction_score))^2)

# step-wise addition of features
set.seed(42)
data_for_model <- loops1[,c(7:13,16,285)]
features <- colnames(loops1)[14:length(colnames(loops1))]
features <- features[c(-3,-272)]
all_models <- tibble(features=features)
step <- tibble(features=features)
rsquared <- c()
for (a in 1:10){
  print(paste0("starting step ",a))
  step$rsquared <- 0
  step$r1 <- 0
  step$r2 <- 0
  step$r3 <- 0
 for (i in features){
  print(paste0("analyzing ",i, match(i, features)))
  data_for_model[i] <- loops1[i]
  ind <- sample(3, nrow(data_for_model), replace = TRUE, prob = c(0.33, 0.33, 0.33))
  test_1 <- data_for_model[ind==1,]
  test_2 <- data_for_model[ind==2,]
  test_3 <- data_for_model[ind==3,]
  model_1 <- ranger(interaction_score~., data=bind_rows(test_2, test_3), num.trees=100, mtry = 3)
  model_2 <- ranger(interaction_score~., data=bind_rows(test_1, test_3), num.trees=100, mtry = 3)
  model_3 <- ranger(interaction_score~., data=bind_rows(test_1, test_2), num.trees=100, mtry = 3)
  step$r1[step$features == i] <- (cor(predict(model_1, test_1)$predictions, test_1$interaction_score))^2
  step$r2[step$features == i] <- (cor(predict(model_2, test_2)$predictions, test_2$interaction_score))^2
  step$r3[step$features == i] <- (cor(predict(model_3, test_3)$predictions, test_3$interaction_score))^2
  step$rsquared[step$features == i] <- mean(c(step$r1[step$features == i], step$r2[step$features == i], step$r3[step$features == i]))
  data_for_model[i] <- NULL}
all_models <- left_join(all_models, step[,1:2], by="features")
step <- arrange(step, desc(rsquared))
rsquared <- append(rsquared,  c(step$r1[1], step$r2[1],step$r3[1]))
data_for_model[step$features[1]] <- loops1[step$features[1]]
features <- setdiff(features, step$features[1])}
# summarizing and writing a df
rsquared <- c(rsq_distance,rsq_initial,rsquared)
write.table(rsquared, "Documents/Jul2022_MChIPC/fig.3/random forest/rsq.txt", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(all_models, "Documents/Jul2022_MChIPC/fig.3/random forest/all_models.txt", quote = F, col.names = F, row.names = F, sep = "\t")
# plot all steps
rsq <- read.csv("Documents/Jul2022_MChIPC/fig.3/random forest/rsq.txt", header = F, sep = " ")
rsq <- tibble(value=as.numeric(rsq[,1]), step = sort(rep(seq(-1,10,1), 3)))
ggplot(rsq) + geom_point(aes(step, value)) + ylim(0,0.45) + scale_x_continuous(breaks=seq(-1,10))

# second step with more details
set.seed(41)
loops <- read.csv("Documents/Jul2022_MChIPC/full_data_for_forests.txt", header = T, sep = "\t")
loops1 <- filter(loops, !(interaction_score %in% sort(loops$interaction_score, decreasing = TRUE)[c(1,2)]))
data_for_model <- loops1[,c(7:13,16,285,22)]
features <- colnames(loops1)[14:length(colnames(loops1))]
features <- features[c(-3,-272, -9)]
all_models <- tibble(features=features)
step <- tibble(features=features)
for (a in 1:1){
  print(paste0("starting step ",a))
  step$rsquared <- 0
  step$r1 <- 0
  step$r2 <- 0
  step$r3 <- 0
  for (i in features){
    print(paste0("analyzing ",i, match(i, features)))
    data_for_model[i] <- loops1[i]
    ind <- sample(3, nrow(data_for_model), replace = TRUE, prob = c(0.33, 0.33, 0.33))
    test_1 <- data_for_model[ind==1,]
    test_2 <- data_for_model[ind==2,]
    test_3 <- data_for_model[ind==3,]
    model_1 <- ranger(interaction_score~., data=bind_rows(test_2, test_3), num.trees=100, mtry = 3)
    model_2 <- ranger(interaction_score~., data=bind_rows(test_1, test_3), num.trees=100, mtry = 3)
    model_3 <- ranger(interaction_score~., data=bind_rows(test_1, test_2), num.trees=100, mtry = 3)
    step$r1[step$features == i] <- (cor(predict(model_1, test_1)$predictions, test_1$interaction_score))^2
    step$r2[step$features == i] <- (cor(predict(model_2, test_2)$predictions, test_2$interaction_score))^2
    step$r3[step$features == i] <- (cor(predict(model_3, test_3)$predictions, test_3$interaction_score))^2
    step$rsquared[step$features == i] <- mean(c(step$r1[step$features == i], step$r2[step$features == i], step$r3[step$features == i]))
    data_for_model[i] <- NULL}
  all_models <- left_join(all_models, step[,1:2], by="features")
  step <- arrange(step, desc(rsquared))}

# plotting the most interesting features at step2
step <- separate(step, features, sep="\\.", into = c(NA, "features", NA, NA))
step$features[step$features=="BRD4"] <- "BRD4.hg19"
step$features[step$features=="cdk8_hg19"] <- "CDK8.hg19"
step$features[step$features=="med_hg19"] <- "MED1.hg19"
bigwig <- read.csv("Documents/Nov2021_MCC/Chips/list_of_beds_and_bigwigs.txt", header = T, sep = "\t")
step <- left_join(step, bigwig[,c(1,3)], by=c("features"="File.accession"))
data_for_plot <- filter(step, target %in% c("EP300","ARID1B","DPF2", "POLR2AphosphoS5","POLR2AphosphoS2", "BRD4","CDK8","YY1","MED1","SMARCE1"))
data_for_plot <- bind_rows(data_for_plot, tibble(features="step1",rsquared=mean(rsq$value[7:9]), r1=rsq$value[7], r2=rsq$value[8], r3=rsq$value[9], target="step1"))
data_for_plot <- melt(data_for_plot, id.vars = "target", measure.vars = c("r1","r2","r3"))
data_for_plot <- data_for_plot %>% mutate(target=factor(target, levels=c("step1","EP300","ARID1B","DPF2","POLR2AphosphoS5", "BRD4", "SMARCE1","CDK8","POLR2AphosphoS2","MED1","YY1")))
write.table(data_for_plot, "Documents/Jul2022_MChIPC/fig.3/random forest/step2.txt", quote = F, col.names = F, row.names = F, sep = "\t")
data_for_plot <- data_for_plot %>% group_by(target) %>% summarise(value=mean(value))
data_for_plot <- arrange(data_for_plot, desc(value)) %>% mutate(target=factor(target, levels=data_for_plot$target))
ggplot(data_for_plot) + geom_point(aes(target, value)) + ylim(0.35,0.45)

# building binary classifier for strong/weak interactions
set.seed(42)
# building distance-only auc model
data_for_model <- loops1[,c(7:8)]
data_for_model$interaction <- data_for_model$interaction_score > median(data_for_model$interaction_score) # using median as a threshold
data_for_model$interaction_score <- NULL
ind <- sample(3, nrow(data_for_model), replace = TRUE, prob = c(0.33, 0.33, 0.33))
test_1 <- data_for_model[ind==1,]
test_2 <- data_for_model[ind==2,]
test_3 <- data_for_model[ind==3,]
model_1 <- ranger(interaction~., data=arrange(bind_rows(test_2, test_3), interaction), num.trees=100, probability = TRUE, mtry = 1)
model_2 <- ranger(interaction~., data=arrange(bind_rows(test_1, test_3), interaction), num.trees=100, probability = TRUE, mtry = 1)
model_3 <- ranger(interaction~., data=arrange(bind_rows(test_1, test_2), interaction), num.trees=100, probability = TRUE, mtry = 1)
auc_distance <- c(roc.curve(scores.class0=predict(model_1, test_1, type="response")$predictions[,2], weights.class0=test_1$interaction)$auc,
                 roc.curve(scores.class0=predict(model_2, test_2, type="response")$predictions[,2], weights.class0=test_2$interaction)$auc,
                 roc.curve(scores.class0=predict(model_3, test_3, type="response")$predictions[,2], weights.class0=test_3$interaction)$auc)

# building initial auc model 
data_for_model <- loops1[,c(7:13,16,285)]
data_for_model$interaction <- data_for_model$interaction_score > median(data_for_model$interaction_score) # using median as a threshold
data_for_model$interaction_score <- NULL
ind <- sample(3, nrow(data_for_model), replace = TRUE, prob = c(0.33, 0.33, 0.33))
test_1 <- data_for_model[ind==1,]
test_2 <- data_for_model[ind==2,]
test_3 <- data_for_model[ind==3,]
model_1 <- ranger(interaction~., data=arrange(bind_rows(test_2, test_3), interaction), num.trees=100, probability = TRUE, mtry = 3)
model_2 <- ranger(interaction~., data=arrange(bind_rows(test_1, test_3), interaction), num.trees=100, probability = TRUE, mtry = 3)
model_3 <- ranger(interaction~., data=arrange(bind_rows(test_1, test_2), interaction), num.trees=100, probability = TRUE, mtry = 3)
auc_initial <- c(roc.curve(scores.class0=predict(model_1, test_1, type="response")$predictions[,2], weights.class0=test_1$interaction)$auc,
                roc.curve(scores.class0=predict(model_2, test_2, type="response")$predictions[,2], weights.class0=test_2$interaction)$auc,
                roc.curve(scores.class0=predict(model_3, test_3, type="response")$predictions[,2], weights.class0=test_3$interaction)$auc)

# step-wise addition of features
data_for_model <- loops1[,c(7:13,16,285)]
data_for_model$interaction <- data_for_model$interaction_score > median(data_for_model$interaction_score) # using median as a threshold
data_for_model$interaction_score <- NULL
features <- colnames(loops1)[14:length(colnames(loops1))]
features <- features[c(-3,-272)]
all_models <- tibble(features=features)
step <- tibble(features=features)
auc <- c()
for (a in 1:1){
  print(paste0("starting step ",a))
  step$auc_roc <- 0
  step$auc_roc1 <- 0
  step$auc_roc2 <- 0
  step$auc_roc3 <- 0
  for (i in features){
    print(paste0("analyzing ",i, match(i, features)))
    data_for_model[i] <- loops1[i]
    ind <- sample(3, nrow(data_for_model), replace = TRUE, prob = c(0.33, 0.33, 0.33))
    test_1 <- data_for_model[ind==1,]
    test_2 <- data_for_model[ind==2,]
    test_3 <- data_for_model[ind==3,]
    model_1 <- ranger(interaction~., data=arrange(bind_rows(test_2, test_3), interaction), num.trees=100, probability = TRUE, mtry = 3)
    model_2 <- ranger(interaction~., data=arrange(bind_rows(test_1, test_3), interaction), num.trees=100, probability = TRUE, mtry = 3)
    model_3 <- ranger(interaction~., data=arrange(bind_rows(test_1, test_2), interaction), num.trees=100, probability = TRUE, mtry = 3)
    step$auc_roc1[step$features == i] <- roc.curve(scores.class0=predict(model_1, test_1, type="response")$predictions[,2], weights.class0=test_1$interaction)$auc
    step$auc_roc2[step$features == i] <- roc.curve(scores.class0=predict(model_2, test_2, type="response")$predictions[,2], weights.class0=test_2$interaction)$auc
    step$auc_roc3[step$features == i] <- roc.curve(scores.class0=predict(model_3, test_3, type="response")$predictions[,2], weights.class0=test_3$interaction)$auc
    data_for_model[i] <- NULL}
  step$auc_roc <- (step$auc_roc1+step$auc_roc2+step$auc_roc3)/3
  all_models <- left_join(all_models, step[,1:2], by="features")
  step <- arrange(step, desc(auc_roc))
  auc <- append(auc,  c(step$auc_roc1[1],step$auc_roc2[1],step$auc_roc3[1]))
  data_for_model[step$features[1]] <- loops1[step$features[1]]
  features <- setdiff(features, step$features[1])}

# summarizing and writing a df
auc <- c(auc_distance,auc_initial,auc)
write.table(all_models, "Documents/Jul2022_MChIPC/fig.3/random forest/all_models_auc.txt", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(auc, "Documents/Jul2022_MChIPC/fig.3/random forest/auc.txt", quote = F, col.names = F, row.names = F, sep = "\t")
# plot all steps
auc <- read.csv("Documents/Jul2022_MChIPC/fig.3/random forest/auc.txt", header = F, sep = " ")
auc <- tibble(value=as.numeric(auc[,1]), step = sort(rep(seq(-1,10,1), 3)))
ggplot(auc) + geom_point(aes(step, value)) + ylim(0.5,0.95) + scale_x_continuous(breaks=seq(-1,10))