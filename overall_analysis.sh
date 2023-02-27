# MChIPC analysis
# biological replicates: rep_0 (3 seq runs), rep_3, rep_4, rep_5
# raw fastq files: MChIPC_rep0_1121 (rep_0), MChIPC_rep0_0122 (rep0), MChIPC_rep0 (rep_0), MChIPC_rep3 (rep_3), MChIPC_rep4 (rep4), MChIPC_rep5 (rep5)
# rep_0 mapping and primary analysis (have to use 'bash' command - not 'sh' to run scripts)
# rep_0_1 (December 2021) 105311061 raw reads, 70140341 within 1kb
../tools/bwa/bwa mem -SP5M -t 8 ../genomes/hg19/BWAIndex/version0.6.0/genome.fa fastqs/MChIPC_rep0_1121_1.fq.gz fastqs/MChIPC_rep0_1121_2.fq.gz | gzip -3 > MChIPC_rep0_1.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../genomes/hg19/hg19.chrom.sizes --output-stats stats_MChIPC_rep0_1_map.txt -o MChIPC_rep0_1.pairsam.gz MChIPC_rep0_1.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o MChIPC_rep0_1.sorted.pairsam.gz MChIPC_rep0_1.pairsam.gz
rm MChIPC_rep0_1.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 5000))' -o MChIPC_rep0_1.cis.pairsam.gz MChIPC_rep0_1.sorted.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 100) and (abs(pos1-pos2) < 200) and (strand1=="+") and (strand2=="-"))' -o ChIP_rep0_1.pairsam.gz MChIPC_rep0_1.sorted.pairsam.gz
rm *.sorted.pairsam.gz
# rep_0_2 (January 2022) 349161883 raw reads, 233365023 within 1 kb
# have to re-name fastq files here
../tools/bwa/bwa mem -SP5M -t 8 ../genomes/hg19/BWAIndex/version0.6.0/genome.fa fastqs/MChIPC_rep0_0122_1.fq.gz fastqs/MChIPC_rep0_0122_2.fq.gz | gzip -3 > MChIPC_rep0_2.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../genomes/hg19/hg19.chrom.sizes --output-stats stats_MChIPC_rep0_2_map.txt -o MChIPC_rep0_2.pairsam.gz MChIPC_rep0_2.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o MChIPC_rep0_2.sorted.pairsam.gz MChIPC_rep0_2.pairsam.gz
rm MChIPC_rep0_2.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 5000))' -o MChIPC_rep0_2.cis.pairsam.gz MChIPC_rep0_2.sorted.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 100) and (abs(pos1-pos2) < 200) and (strand1=="+") and (strand2=="-"))' -o ChIP_rep0_2.pairsam.gz MChIPC_rep0_2.sorted.pairsam.gz
rm *.sorted.pairsam.gz
# rep_0_3 (July 2022) 420867478 raw reads, 312819345 within 1 kb
../tools/bwa/bwa mem -SP5M -t 8 ../genomes/hg19/BWAIndex/version0.6.0/genome.fa fastqs/MChIPC_rep0_1.fq.gz fastqs/MChIPC_rep0_2.fq.gz | gzip -3 > MChIPC_rep0_3.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../genomes/hg19/hg19.chrom.sizes --output-stats stats_MChIPC_rep0_3_map.txt -o MChIPC_rep0_3.pairsam.gz MChIPC_rep0_3.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o MChIPC_rep0_3.sorted.pairsam.gz MChIPC_rep0_3.pairsam.gz
rm MChIPC_rep0_3.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 5000))' -o MChIPC_rep0_3.cis.pairsam.gz MChIPC_rep0_3.sorted.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 100) and (abs(pos1-pos2) < 200) and (strand1=="+") and (strand2=="-"))' -o ChIP_rep0_3.pairsam.gz MChIPC_rep0_3.sorted.pairsam.gz
rm *.sorted.pairsam.gz
# merging ChIP_rep0 technical replicates and performing downstream analysis: downsampling, deduplication & peak calling
pairtools merge --cmd-in gunzip --cmd-out gzip --nproc-in 8 -o ChIP_rep0.pairsam.gz ChIP_rep0_1.pairsam.gz ChIP_rep0_2.pairsam.gz ChIP_rep0_3.pairsam.gz
nreads=$(gunzip -c ChIP_rep0.pairsam.gz | wc -l)
prop=$(bc <<< "scale=2; 20000000/$nreads")
pairtools sample --cmd-in gunzip --cmd-out gzip --nproc-in 8 -o ChIP_rep0.subsampled.pairsam.gz $prop ChIP_rep0.pairsam.gz
pairtools dedup --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-stats stats.ChIP_rep0.dedup.txt -o ChIP_output/rep_0/ChIP_rep0.dedup.pairsam.gz ChIP_rep0.subsampled.pairsam.gz  
pairtools split --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-sam ChIP_output/rep_0/ChIP_rep0.sam.gz ChIP_output/rep_0/ChIP_rep0.dedup.pairsam.gz
# I have to add flag 2 (properly paired) for MACS2 to work normally
samtools view -@ 7 -bh ChIP_output/rep_0/ChIP_rep0.sam.gz --add-flag 2 > ChIP_output/rep_0/ChIP_rep0.bam
samtools sort -@ 7 ChIP_output/rep_0/ChIP_rep0.bam > ChIP_output/rep_0/ChIP_rep0.sorted.bam
samtools index -@ 7 ChIP_output/rep_0/ChIP_rep0.sorted.bam
bamCoverage -b ChIP_output/rep_0/ChIP_rep0.sorted.bam -bs 50 -e -p 7 -o ChIP_output/rep_0/ChIP_rep0.bw
macs2 callpeak -t ChIP_output/rep_0/ChIP_rep0.sorted.bam --outdir ChIP_output/rep_0/macs/ -n ChIP_rep0 -f BAMPE -q 0.0001 --max-gap 1000
# merging MChIPC_rep0 technical replicates and all downstream
pairtools merge --cmd-in gunzip --cmd-out gzip --nproc-in 8 -o MChIPC_rep0.cis.pairsam.gz MChIPC_rep0_1.cis.pairsam.gz MChIPC_rep0_2.cis.pairsam.gz MChIPC_rep0_3.cis.pairsam.gz
pairtools dedup --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-stats stats.MChIPC_rep0.dedup.txt -o MChIPC_output/rep_0/MChIPC_rep0.dedup.pairsam.gz MChIPC_rep0.cis.pairsam.gz
pairtools split --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-sam MChIPC_rep0.sam.gz MChIPC_rep0.dedup.pairsam.gz
samtools view -@ 7 -bh MChIPC_rep0.sam.gz > MChIPC_rep0.bam
samtools sort -@ 7 MChIPC_rep0.bam > MChIPC_rep0.sorted.bam
samtools index -@ 7 MChIPC_rep0.sorted.bam
samtools view -b -P -h -@ 7 --region-file ../../binned_peaks.bed MChIPC_rep0.sorted.bam > MChIPC_rep0.filtered.bam
samtools index -@ 7 MChIPC_rep0.filtered.bam
bamCoverage -b MChIPC_rep0.filtered.bam -bs 250 -o MChIPC_rep0.bw
pairtools split --cmd-in gunzip --cmd-out bgzip --nproc-in 8 --output-pairs MChIPC_rep0.pairs.gz MChIPC_rep0.dedup.pairsam.gz
# separate script for creating 'interactions' bedpe file
file="MChIPC_rep0.pairs.gz"
rep="rep0"
gunzip -c $file | grep -v '^#'| awk '{OFS="\t"}{print $2,$3,$3+1,$4,$5,$5+1}' > contacts.$rep.bedpe
# overlapping contacts with binned peaks
bedtools pairtobed -b ../../binned_peaks.bed -a contacts.$rep.bedpe -bedpe > contacts.peakoverlap.$rep.bedpe
#flipping contacts with baits last and removing coordinates of bait read
awk '{OFS="\t"} {if ($2<=$8) {print $1,$2,$3,$7,$8,$9} else {print $4,$5,$6,$7,$8,$9}}' contacts.peakoverlap.$rep.bedpe > contacts.peakoverlap.flipped.$rep.bedpe
# overlapping OE with bins (do not understand sort completely)
bedtools intersect -wo -a contacts.peakoverlap.flipped.$rep.bedpe -b ../../genomic_bins.bed | sort -k4,4 -k5,8n -k8,8n | cut -f 4,5,6,7,8,9 | uniq -c | awk {'OFS="\t"}{print $2,$3,$4,$5,$6,$7,$1}' > interactions.sum.$rep.txt
rm contacts.$rep.bedpe contacts.peakoverlap.$rep.bedpe contacts.peakoverlap.flipped.$rep.bedpe

# rep_3 mapping and primary analysis 625108983 raw reads, 415260796 within 1 kb
../tools/bwa/bwa mem -SP5M -t 8 ../genomes/hg19/BWAIndex/version0.6.0/genome.fa fastqs/MChIPC_rep3_1.fq.gz fastqs/MChIPC_rep3_2.fq.gz | gzip -3 > MChIPC_rep3.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../genomes/hg19/hg19.chrom.sizes --output-stats stats_MChIPC_rep3_map.txt -o MChIPC_rep3.pairsam.gz MChIPC_rep3.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o MChIPC_rep3.sorted.pairsam.gz MChIPC_rep3.pairsam.gz
rm MChIPC_rep3.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 5000))' -o MChIPC_rep3.cis.pairsam.gz MChIPC_rep3.sorted.pairsam.gz
pairtools dedup --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-stats stats.MChIPC_rep3.dedup.txt -o MChIPC_output/MChIPC_rep3.dedup.pairsam.gz MChIPC_rep3.cis.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 100) and (abs(pos1-pos2) < 200) and (strand1=="+") and (strand2=="-"))' -o ChIP_rep3.pairsam.gz MChIPC_rep3.sorted.pairsam.gz
rm *.sorted.pairsam.gz
nreads=$(gunzip -c ChIP_rep3.pairsam.gz | wc -l)
prop=$(bc <<< "scale=2; 20000000/$nreads")
pairtools sample --cmd-in gunzip --cmd-out gzip --nproc-in 8 -o ChIP_rep3.subsampled.pairsam.gz $prop ChIP_rep3.pairsam.gz
pairtools dedup --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-stats ChIP_output/rep_3/stats.ChIP_rep3.dedup.txt -o ChIP_output/rep_3/ChIP_rep3.dedup.pairsam.gz ChIP_rep3.subsampled.pairsam.gz
pairtools split --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-sam ChIP_output/rep_3/ChIP_rep3.sam.gz ChIP_output/rep_3/ChIP_rep3.dedup.pairsam.gz
# I have to add flag 2 (properly paired) for MACS2 to work normally
samtools view -@ 7 -bh ChIP_output/rep_3/ChIP_rep3.sam.gz --add-flag 2 > ChIP_output/rep_3/ChIP_rep3.bam
samtools sort -@ 7 ChIP_output/rep_3/ChIP_rep3.bam > ChIP_output/rep_3/ChIP_rep3.sorted.bam
samtools index -@ 7 ChIP_output/rep_3/ChIP_rep3.sorted.bam
bamCoverage -b ChIP_output/rep_3/ChIP_rep3.sorted.bam -bs 50 -e -p 7 -o ChIP_output/rep_3/ChIP_rep3.bw
macs2 callpeak -t ChIP_output/rep_3/ChIP_rep3.sorted.bam --outdir ChIP_output/rep_3/macs/ -n ChIP_rep3 -f BAMPE -q 0.0001 --max-gap 1000
# downstream MChIPC_rep3
pairtools dedup --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-stats stats.MChIPC_rep3.dedup.txt -o MChIPC_output/rep_3/MChIPC_rep3.dedup.pairsam.gz MChIPC_rep3.cis.pairsam.gz
pairtools split --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-sam MChIPC_rep3.sam.gz MChIPC_rep3.dedup.pairsam.gz
samtools view -@ 7 -bh MChIPC_rep3.sam.gz > MChIPC_rep3.bam
samtools sort -@ 7 MChIPC_rep3.bam > MChIPC_rep3.sorted.bam
samtools index -@ 7 MChIPC_rep3.sorted.bam
samtools view -b -P -h -@ 7 --region-file ../../binned_peaks.bed MChIPC_rep3.sorted.bam > MChIPC_rep3.filtered.bam
samtools index -@ 7 MChIPC_rep3.filtered.bam
bamCoverage -b MChIPC_rep3.filtered.bam -bs 250 -o MChIPC_rep3.bw
pairtools split --cmd-in gunzip --cmd-out bgzip --nproc-in 8 --output-pairs MChIPC_rep3.pairs.gz MChIPC_rep3.dedup.pairsam.gz
# separate script for creating 'interactions' bedpe file
file="MChIPC_rep3.pairs.gz"
rep="rep3"
gunzip -c $file | grep -v '^#'| awk '{OFS="\t"}{print $2,$3,$3+1,$4,$5,$5+1}' > contacts.$rep.bedpe
# overlapping contacts with binned peaks
bedtools pairtobed -b ../../binned_peaks.bed -a contacts.$rep.bedpe -bedpe > contacts.peakoverlap.$rep.bedpe
#flipping contacts with baits last and removing coordinates of bait read
awk '{OFS="\t"} {if ($2<=$8) {print $1,$2,$3,$7,$8,$9} else {print $4,$5,$6,$7,$8,$9}}' contacts.peakoverlap.$rep.bedpe > contacts.peakoverlap.flipped.$rep.bedpe
# overlapping OE with bins (do not understand sort completely)
bedtools intersect -wo -a contacts.peakoverlap.flipped.$rep.bedpe -b ../../genomic_bins.bed | sort -k4,4 -k5,8n -k8,8n | cut -f 4,5,6,7,8,9 | uniq -c | awk {'OFS="\t"}{print $2,$3,$4,$5,$6,$7,$1}' > interactions.sum.$rep.txt
rm contacts.$rep.bedpe contacts.peakoverlap.$rep.bedpe contacts.peakoverlap.flipped.$rep.bedpe

# rep_4 mapping and primary analysis 546281204 raw reads, 343583426 within 1kb
../tools/bwa/bwa mem -SP5M -t 8 ../genomes/hg19/BWAIndex/version0.6.0/genome.fa fastqs/MChIPC_rep4_1.fq.gz fastqs/MChIPC_rep4_2.fq.gz | gzip -3 > MChIPC_rep4.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../genomes/hg19/hg19.chrom.sizes --output-stats stats_MChIPC_rep4_map.txt -o MChIPC_rep4.pairsam.gz MChIPC_rep4.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o MChIPC_rep4.sorted.pairsam.gz MChIPC_rep4.pairsam.gz
rm MChIPC_rep4.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 5000))' -o MChIPC_rep4.cis.pairsam.gz MChIPC_rep4.sorted.pairsam.gz
pairtools dedup --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-stats stats.MChIPC_rep4.dedup.txt -o MChIPC_output/MChIPC_rep4.dedup.pairsam.gz MChIPC_rep4.cis.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 100) and (abs(pos1-pos2) < 200) and (strand1=="+") and (strand2=="-"))' -o ChIP_rep4.pairsam.gz MChIPC_rep4.sorted.pairsam.gz
rm *.sorted.pairsam.gz
nreads=$(gunzip -c ChIP_rep4.pairsam.gz | wc -l)
prop=$(bc <<< "scale=2; 20000000/$nreads")
pairtools sample --cmd-in gunzip --cmd-out gzip --nproc-in 8 -o ChIP_rep4.subsampled.pairsam.gz $prop ChIP_rep4.pairsam.gz
pairtools dedup --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-stats ChIP_output/rep_4/stats.ChIP_rep4.dedup.txt -o ChIP_output/rep_4/ChIP_rep4.dedup.pairsam.gz ChIP_rep4.subsampled.pairsam.gz
pairtools split --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-sam ChIP_output/rep_4/ChIP_rep4.sam.gz ChIP_output/rep_4/ChIP_rep4.dedup.pairsam.gz
# I have to add flag 2 (properly paired) for MACS2 to work normally
samtools view -@ 7 -bh ChIP_output/rep_4/ChIP_rep4.sam.gz --add-flag 2 > ChIP_output/rep_4/ChIP_rep4.bam
samtools sort -@ 7 ChIP_output/rep_4/ChIP_rep4.bam > ChIP_output/rep_4/ChIP_rep4.sorted.bam
samtools index -@ 7 ChIP_output/rep_4/ChIP_rep4.sorted.bam
bamCoverage -b ChIP_output/rep_4/ChIP_rep4.sorted.bam -bs 50 -e -p 7 -o ChIP_output/rep_4/ChIP_rep4.bw
macs2 callpeak -t ChIP_output/rep_4/ChIP_rep4.sorted.bam --outdir ChIP_output/rep_4/macs/ -n ChIP_rep4 -f BAMPE -q 0.0001 --max-gap 1000
# downstream MChIPC_rep4
pairtools dedup --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-stats stats.MChIPC_rep4.dedup.txt -o MChIPC_output/rep_4/MChIPC_rep4.dedup.pairsam.gz MChIPC_rep4.cis.pairsam.gz
pairtools split --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-sam MChIPC_rep4.sam.gz MChIPC_rep4.dedup.pairsam.gz
samtools view -@ 7 -bh MChIPC_rep4.sam.gz > MChIPC_rep4.bam
samtools sort -@ 7 MChIPC_rep4.bam > MChIPC_rep4.sorted.bam
samtools index -@ 7 MChIPC_rep4.sorted.bam
samtools view -b -P -h -@ 7 --region-file ../../binned_peaks.bed MChIPC_rep4.sorted.bam > MChIPC_rep4.filtered.bam
samtools index -@ 7 MChIPC_rep4.filtered.bam
bamCoverage -b MChIPC_rep4.filtered.bam -bs 250 -o MChIPC_rep4.bw
pairtools split --cmd-in gunzip --cmd-out bgzip --nproc-in 8 --output-pairs MChIPC_rep4.pairs.gz MChIPC_rep4.dedup.pairsam.gz
# separate script for creating 'interactions' bedpe file
file="MChIPC_rep4.pairs.gz"
rep="rep4"
gunzip -c $file | grep -v '^#'| awk '{OFS="\t"}{print $2,$3,$3+1,$4,$5,$5+1}' > contacts.$rep.bedpe
# overlapping contacts with binned peaks
bedtools pairtobed -b ../../binned_peaks.bed -a contacts.$rep.bedpe -bedpe > contacts.peakoverlap.$rep.bedpe
#flipping contacts with baits last and removing coordinates of bait read
awk '{OFS="\t"} {if ($2<=$8) {print $1,$2,$3,$7,$8,$9} else {print $4,$5,$6,$7,$8,$9}}' contacts.peakoverlap.$rep.bedpe > contacts.peakoverlap.flipped.$rep.bedpe
# overlapping OE with bins (do not understand sort completely)
bedtools intersect -wo -a contacts.peakoverlap.flipped.$rep.bedpe -b ../../genomic_bins.bed | sort -k4,4 -k5,8n -k8,8n | cut -f 4,5,6,7,8,9 | uniq -c | awk {'OFS="\t"}{print $2,$3,$4,$5,$6,$7,$1}' > interactions.sum.$rep.txt
rm contacts.$rep.bedpe contacts.peakoverlap.$rep.bedpe contacts.peakoverlap.flipped.$rep.bedpe

# rep_5 mapping and primary analysis
# first, split (to solve HD space issues)
fastqsplitter -i MChIPC_rep5_1.fq.gz -o MChIPC_rep5_1.1.fq.gz -o MChIPC_rep5_2.1.fq.gz
fastqsplitter -i MChIPC_rep5_2.fq.gz -o MChIPC_rep5_1.2.fq.gz -o MChIPC_rep5_2.2.fq.gz
# dealing with rep5_1 pseudo-replicate 747866194 raw reads, 595062674 within 1kb
../tools/bwa/bwa mem -SP5M -t 8 ../genomes/hg19/BWAIndex/version0.6.0/genome.fa fastqs/MChIPC_rep5_1.1.fq.gz fastqs/MChIPC_rep5_1.2.fq.gz | gzip -3 > MChIPC_rep5_1.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../genomes/hg19/hg19.chrom.sizes --output-stats stats_MChIPC_rep5_1_map.txt -o MChIPC_rep5_1.pairsam.gz MChIPC_rep5_1.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o MChIPC_rep5_1.sorted.pairsam.gz MChIPC_rep5_1.pairsam.gz
rm MChIPC_rep5_1.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 5000))' -o MChIPC_rep5_1.cis.pairsam.gz MChIPC_rep5_1.sorted.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 100) and (abs(pos1-pos2) < 200) and (strand1=="+") and (strand2=="-"))' -o ChIP_rep5_1.pairsam.gz MChIPC_rep5_1.sorted.pairsam.gz
rm *.sorted.pairsam.gz
# rep5_2 747866100 raw reads, 595074839 within 1kb
../tools/bwa/bwa mem -SP5M -t 8 ../genomes/hg19/BWAIndex/version0.6.0/genome.fa fastqs/MChIPC_rep5_2.1.fq.gz fastqs/MChIPC_rep5_2.2.fq.gz | gzip -3 > MChIPC_rep5_2.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../genomes/hg19/hg19.chrom.sizes --output-stats stats_MChIPC_rep5_2_map.txt -o MChIPC_rep5_2.pairsam.gz MChIPC_rep5_2.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o MChIPC_rep5_2.sorted.pairsam.gz MChIPC_rep5_2.pairsam.gz
rm MChIPC_rep5_2.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 5000))' -o MChIPC_rep5_2.cis.pairsam.gz MChIPC_rep5_2.sorted.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 100) and (abs(pos1-pos2) < 200) and (strand1=="+") and (strand2=="-"))' -o ChIP_rep5_2.pairsam.gz MChIPC_rep5_2.sorted.pairsam.gz
rm *.sorted.pairsam.gz
# merging pseudoreplicates and performing downstream analysis
pairtools merge --cmd-in gunzip --cmd-out gzip --nproc-in 8 -o ChIP_rep5.pairsam.gz ChIP_rep5_1.pairsam.gz ChIP_rep5_2.pairsam.gz
nreads=$(gunzip -c ChIP_rep5.pairsam.gz | wc -l)
prop=$(bc <<< "scale=2; 20000000/$nreads")
pairtools sample --cmd-in gunzip --cmd-out gzip --nproc-in 8 -o ChIP_rep5.subsampled.pairsam.gz $prop ChIP_rep5.pairsam.gz
pairtools dedup --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-stats stats.ChIP_rep5.dedup.txt -o ChIP_output/rep_5/ChIP_rep5.dedup.pairsam.gz ChIP_rep5.subsampled.pairsam.gz  
pairtools split --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-sam ChIP_output/rep_5/ChIP_rep5.sam.gz ChIP_output/rep_5/ChIP_rep5.dedup.pairsam.gz
# I have to add flag 2 (properly paired) for MACS2 to work normally
samtools view -@ 7 -bh ChIP_output/rep_5/ChIP_rep5.sam.gz --add-flag 2 > ChIP_output/rep_5/ChIP_rep5.bam
samtools sort -@ 7 ChIP_output/rep_5/ChIP_rep5.bam > ChIP_output/rep_5/ChIP_rep5.sorted.bam
samtools index -@ 7 ChIP_output/rep_5/ChIP_rep5.sorted.bam
bamCoverage -b ChIP_output/rep_5/ChIP_rep5.sorted.bam -bs 50 -e -p 7 -o ChIP_output/rep_5/ChIP_rep5.bw
macs2 callpeak -t ChIP_output/rep_5/ChIP_rep5.sorted.bam --outdir ChIP_output/rep_5/macs/ -n ChIP_rep5 -f BAMPE -q 0.0001 --max-gap 1000
# merging MChIPC_rep5 pseudo-replicates and all downstream
pairtools merge --cmd-in gunzip --cmd-out gzip --nproc-in 8 -o MChIPC_rep5.cis.pairsam.gz MChIPC_rep5_1.cis.pairsam.gz MChIPC_rep5_2.cis.pairsam.gz
pairtools dedup --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-stats stats.MChIPC_rep5.dedup.txt -o MChIPC_output/rep_5/MChIPC_rep5.dedup.pairsam.gz MChIPC_rep5.cis.pairsam.gz
pairtools split --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-sam MChIPC_rep5.sam.gz MChIPC_rep5.dedup.pairsam.gz
samtools view -@ 7 -bh MChIPC_rep5.sam.gz > MChIPC_rep5.bam
samtools sort -@ 7 MChIPC_rep5.bam > MChIPC_rep5.sorted.bam
samtools index -@ 7 MChIPC_rep5.sorted.bam
samtools view -b -P -h -@ 7 --region-file ../../binned_peaks.bed MChIPC_rep5.sorted.bam > MChIPC_rep5.filtered.bam
samtools index -@ 7 MChIPC_rep5.filtered.bam
bamCoverage -b MChIPC_rep5.filtered.bam -bs 250 -o MChIPC_rep5.bw
pairtools split --cmd-in gunzip --cmd-out bgzip --nproc-in 8 --output-pairs MChIPC_rep5.pairs.gz MChIPC_rep5.dedup.pairsam.gz

# merging all replicates to create raw and merged bw files
samtools merge -@ 7 -o MChIPC_raw.bam rep_0/MChIPC_rep0.sorted.bam rep_3/MChIPC_rep3.sorted.bam rep_4/MChIPC_rep4.sorted.bam rep_5/MChIPC_rep5.sorted.bam
samtools index -@ 7 MChIPC_raw.bam
bamCoverage -p 7 -b MChIPC_raw.bam -bs 250 -o MChIPC_raw.bw
bamCoverage -p 7 -b MChIPC_raw.bam -bs 100 -o MChIPC_raw.100.bw
samtools merge -@ 7 -o MChIPC_merged.bam rep_0/MChIPC_rep0.filtered.bam rep_3/MChIPC_rep3.filtered.bam rep_4/MChIPC_rep4.filtered.bam rep_5/MChIPC_rep5.filtered.bam
samtools index -@ 7 MChIPC_merged.bam
bamCoverage -p 7 -b MChIPC_merged.bam -bs 250 -o MChIPC_merged.bw
bamCoverage -p 7 -b MChIPC_merged.bam -bs 100 -o MChIPC_merged.100.bw

# separate script for creating 'interactions' bedpe file
file="MChIPC_rep5.pairs.gz"
rep="rep5"
gunzip -c $file | grep -v '^#'| awk '{OFS="\t"}{print $2,$3,$3+1,$4,$5,$5+1}' > contacts.$rep.bedpe
# overlapping contacts with binned peaks
bedtools pairtobed -b ../../binned_peaks.bed -a contacts.$rep.bedpe -bedpe > contacts.peakoverlap.$rep.bedpe
#flipping contacts with baits last and removing coordinates of bait read
awk '{OFS="\t"} {if ($2<=$8) {print $1,$2,$3,$7,$8,$9} else {print $4,$5,$6,$7,$8,$9}}' contacts.peakoverlap.$rep.bedpe > contacts.peakoverlap.flipped.$rep.bedpe
# overlapping OE with bins (do not understand sort completely)
bedtools intersect -wo -a contacts.peakoverlap.flipped.$rep.bedpe -b ../../genomic_bins.bed | sort -k4,4 -k5,8n -k8,8n | cut -f 4,5,6,7,8,9 | uniq -c | awk {'OFS="\t"}{print $2,$3,$4,$5,$6,$7,$1}' > interactions.sum.$rep.txt
rm contacts.$rep.bedpe contacts.peakoverlap.$rep.bedpe contacts.peakoverlap.flipped.$rep.bedpe

# bining the genome 
bedtools makewindows -g ../../../hg19/hg19.chrom.sizes -w 250 > genomic_bins.bed

# finding concensus peaks and creating H3K4me3 mask
bedtools multiinter -i ChIP_rep0_peaks.narrowPeak ChIP_rep3_peaks.narrowPeak ChIP_rep4_peaks.narrowPeak ChIP_rep5_peaks.narrowPeak | awk '$4 >= 3 {print $1"\t"$2"\t"$3}' | bedtools merge -d 1000 -i stdin | bedtools slop -b 750 -i stdin -g ../../../hg19/hg19.txt > concensus_peaks.bed
bedtools multiinter -i ChIP_rep0_peaks.narrowPeak ChIP_rep3_peaks.narrowPeak ChIP_rep4_peaks.narrowPeak ChIP_rep5_peaks.narrowPeak  peak.overlaps.bed | awk '{print $1"\t"$2"\t"$3}' | bedtools merge -i stdin > H3K4me3.mask.bed
# creating bait-file (binned_peaks)
bedtools intersect -wa -a genomic_bins.bed -b concensus_peaks.bed | bedtools merge > binned_peaks.bed


# showing that found peaks are promoters
# first, downloading CAGE files for K562 cells fromm ENCODE (ENCFF019NPZ - CAGE-plus-rep1.bigWig, ENCFF550YWP - CAGE-plus-rep2.bigWig, ENCFF233CVF CAGE-minus-rep1.bigWig, ENCFF570ONH -CAGE-minus-rep2.bigWig) 
# averaging replicates with deeptools:
bigwigCompare -b1 ENCFF019NPZ\ -\ CAGE-plus-rep1.bigWig -b2 ENCFF550YWP\ -\ CAGE-plus-rep2.bigWig --operation mean -p 3 -o CAGE_plus.bw
bigwigCompare -b1 ENCFF233CVF\ CAGE-minus-rep1.bigWig -b2 ENCFF570ONH\ -CAGE-minus-rep2.bigWig --operation mean -p 3 -o CAGE_minus.bw
# creating new DNase without promoter-specific sites with 'creating_DNase_for_CAGE.R' script
# creating matrix and plotting heatmap
computeMatrix reference-point --referencePoint center -S ../../../Nov2021_MCC/Chips/ENCFF352SET\ -\ DNase.bigWig ../../ChIPs/ChIP_mean.bw CAGE_plus.bw CAGE_minus.bw  -R ../../ChIPs/macs/binned_peaks.bed DNase_outside_promoters.bed -p 3 -a 5000 -b 5000 -o matrix_CAGE_250 -bs 250 --missingDataAsZero
plotHeatmap -m matrix_CAGE_250 --colorMap OrRd -o heatmap_CAGE.pdf --zMax 1.2 65 0.001 0.001 --yMax 1.2 65 0.03 0.03 --yMin 0

# estimating coverage for baits (average coverage per 250bp) and all genomic bins for all 4 replicates
samtools bedcov binned_peaks.bed rep_0/ChIP_rep0.sorted.bam > rep_0/peaks_coverage_rep0.bed
awk '{OFS="\t"}{$5=$4*250/($3-$2)}{print}' peaks_coverage_rep0.bed > tmp.bed && mv -f tmp.bed peaks_coverage_rep0.bed
samtools bedcov genomic_bins.bed rep_0/ChIP_rep0.sorted.bam > rep_0/genomic_bins_coverage_rep0.bed
samtools bedcov binned_peaks.bed rep_3/ChIP_rep3.sorted.bam > rep_3/peaks_coverage_rep3.bed
awk '{OFS="\t"}{$5=$4*250/($3-$2)}{print}' peaks_coverage_rep3.bed > tmp.bed && mv -f tmp.bed peaks_coverage_rep3.bed
samtools bedcov genomic_bins.bed rep_3/ChIP_rep3.sorted.bam > rep_3/genomic_bins_coverage_rep3.bed
samtools bedcov binned_peaks.bed rep_4/ChIP_rep4.sorted.bam > rep_4/peaks_coverage_rep4.bed
awk '{OFS="\t"}{$5=$4*250/($3-$2)}{print}' peaks_coverage_rep4.bed > tmp.bed && mv -f tmp.bed peaks_coverage_rep4.bed
samtools bedcov genomic_bins.bed rep_4/ChIP_rep4.sorted.bam > rep_4/genomic_bins_coverage_rep4.bed
samtools bedcov binned_peaks.bed rep_5/ChIP_rep5.sorted.bam > rep_5/peaks_coverage_rep5.bed
awk '{OFS="\t"}{$5=$4*250/($3-$2)}{print}' peaks_coverage_rep5.bed > tmp.bed && mv -f tmp.bed peaks_coverage_rep5.bed
samtools bedcov genomic_bins.bed rep_5/ChIP_rep5.sorted.bam > rep_5/genomic_bins_coverage_rep5.bed

# mapping conventional ChIP-seq (from 11/2021) (now, probably done this with bowtie2 pipeline)
../../tools/bwa/bwa mem -SP5M -t 8 ../../genomes/hg19/BWAIndex/version0.6.0/genome.fa ChIP_H3K4me3_1.fq.gz ChIP_H3K4me3_2.fq.gz | gzip -3 > ChIP_H3K4me3.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../../genomes/hg19/hg19.chrom.sizes --output-stats stats_ChIP_H3K4me3_map.txt -o ChIP_H3K4me3.pairsam.gz ChIP_H3K4me3.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o ChIP_H3K4me3.sorted.pairsam.gz ChIP_H3K4me3.pairsam.gz
rm ChIP_H3K4me3.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 100) and (abs(pos1-pos2) < 200) and (strand1=="+") and (strand2=="-"))' -o ChIP_H3K4me3.pairsam.gz ChIP_H3K4me3.sorted.pairsam.gz
pairtools dedup --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-stats stats.ChIP_H3K4me3.dedup.txt -o ChIP_H3K4me3.dedup.pairsam.gz ChIP_H3K4me3.pairsam.gz  
pairtools split --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-sam ChIP_H3K4me3.sam.gz ChIP_H3K4me3.dedup.pairsam.gz
samtools view -@ 7 -bh ChIP_H3K4me3.sam.gz --add-flag 2 > ChIP_H3K4me3.bam
samtools sort -@ 7 ChIP_H3K4me3.bam > ChIP_H3K4me3.sorted.bam
samtools index -@ 7 ChIP_H3K4me3.sorted.bam
bamCoverage -b ChIP_H3K4me3.sorted.bam -bs 50 -e -p 7 -o ChIP_H3K4me3.bw

# now I can take DHSs and look at correlations between ChIP-seq signal in MChIP-C and conv.ChIP in these sites
# all analysis in script "correlations_ChIP.R"
# correlations: rep0,rep3=0.96 rep0,rep4=0.96 rep0,rep5=0.96 rep0,ChIP=0.90 rep3,rep4=0.97 rep3,rep5=0.95 rep3,ChIP=0.91 rep4,rep5=0.95 rep4,ChIP=0.90 rep5,ChIP=0.89 

# coreelation analysis of MChIP-C biological replicates
# first, created coverage dataframe with deeptools:
multiBigwigSummary bins -bs 250 -p 3 -b MChIPC_rep0.bw MChIPC_rep3.bw MChIPC_rep4.bw MChIPC_rep5.bw --outRawCounts fig.S1/replicates_coverage.txt -o results.npz
# then re-direct data to R script "correlations_MChIP-C.R" to calulate correlations and plot graphs
# corelations:
# separate bins: rep0,rep3=0.68 rep0,rep4=0.65 rep0,rep5=0.76 rep3,rep4=0.64 rep3,rep5=0.74 rep4,rep5=0.69
# merged data: rep0,rep3=0.84 rep0,rep4=0.83 rep0,rep5=0.85 rep3,rep4=0.91 rep3,rep5=0.93 rep4,rep5=0.90 
# also assessing distance decay in replicates "distance_distr_in_replicates.R"

# now everything is ready for MChIPC loop calling (see 'MChIPC_loop_calling.R')


# PLAC-seq re-analysis
# downloading fastq files for both replicates from tps://data.4dnucleome.org/experiment-set-replicates/4DNESWX1J3QU/#raw-file

# rep1
../../tools/bwa/bwa mem -SP5M -t 8 ../../genomes/hg19/BWAIndex/version0.6.0/genome.fa 4DNFIM5BWS5H.fastq.gz 4DNFITMUHEBY.fastq.gz | gzip -3 > PLACseq_rep1_1.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../../genomes/hg19/hg19.chrom.sizes --output-stats stats_PLACseq_rep1_1_map.txt -o PLACseq_rep1_1.pairsam.gz PLACseq_rep1_1.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o PLACseq_rep1_1.sorted.pairsam.gz PLACseq_rep1_1.pairsam.gz
rm PLACseq_rep1_1.pairsam.gz
../../tools/bwa/bwa mem -SP5M -t 8 ../../genomes/hg19/BWAIndex/version0.6.0/genome.fa 4DNFI9TNKOBD.fastq.gz 4DNFI2I5A8ZU.fastq.gz | gzip -3 > PLACseq_rep1_2.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../../genomes/hg19/hg19.chrom.sizes --output-stats stats_PLACseq_rep1_2_map.txt -o PLACseq_rep1_2.pairsam.gz PLACseq_rep1_2.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o PLACseq_rep1_2.sorted.pairsam.gz PLACseq_rep1_2.pairsam.gz
rm PLACseq_rep1_2.pairsam.gz
../../tools/bwa/bwa mem -SP5M -t 8 ../../genomes/hg19/BWAIndex/version0.6.0/genome.fa 4DNFIZEP644W.fastq.gz 4DNFIMWVLMTQ.fastq.gz | gzip -3 > PLACseq_rep1_3.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../../genomes/hg19/hg19.chrom.sizes --output-stats stats_PLACseq_rep1_3_map.txt -o PLACseq_rep1_3.pairsam.gz PLACseq_rep1_3.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o PLACseq_rep1_3.sorted.pairsam.gz PLACseq_rep1_3.pairsam.gz
rm PLACseq_rep1_3.pairsam.gz
pairtools merge --nproc-in 8 --cmd-in gunzip --cmd-out gzip --nproc 8 -o PLACseq_rep1.merged.sorted.pairsam.gz *.sorted.pairsam.gz
pairtools dedup --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-stats stats.PLACseq_rep1.dedup.txt -o PLACseq_rep1.dedup.pairsam.gz PLACseq_rep1_merged.sorted.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 5000))' -o PLACseq_rep1.cis.pairsam.gz PLACseq_rep1.dedup.pairsam.gz
pairtools split --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-sam PLACseq_rep1.cis.sam.gz PLACseq_rep1.cis.pairsam.gz
samtools view -@ 7 -bh PLACseq_rep1.cis.sam.gz > PLACseq_rep1.cis.bam
samtools sort -@ 7 PLACseq_rep1.cis.bam > PLACseq_rep1.cis.sorted.bam
samtools index -@ 7 PLACseq_rep1.cis.sorted.bam
samtools view -b -P -h -@ 7 --region-file ../../MChIPC/binned_peaks.bed PLACseq_rep1.cis.sorted.bam > PLACseq_rep1.cis.filtered.bam
samtools index -@ 7 PLACseq_rep1.cis.filtered.bam
# rep2
../../tools/bwa/bwa mem -SP5M -t 8 ../../genomes/hg19/BWAIndex/version0.6.0/genome.fa 4DNFIP2TS3O8.fastq.gz 4DNFIU98V4FB.fastq.gz | gzip -3 > PLACseq_rep2_1.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../../genomes/hg19/hg19.chrom.sizes --output-stats stats_PLACseq_rep2_1_map.txt -o PLACseq_rep2_1.pairsam.gz PLACseq_rep2_1.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o PLACseq_rep2_1.sorted.pairsam.gz PLACseq_rep2_1.pairsam.gz
rm PLACseq_rep2_1.pairsam.gz
../../tools/bwa/bwa mem -SP5M -t 8 ../../genomes/hg19/BWAIndex/version0.6.0/genome.fa 4DNFI69HRUMH.fastq.gz 4DNFIMTJX2X7.fastq.gz | gzip -3 > PLACseq_rep2_2.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../../genomes/hg19/hg19.chrom.sizes --output-stats stats_PLACseq_rep2_2_map.txt -o PLACseq_rep2_2.pairsam.gz PLACseq_rep2_2.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o PLACseq_rep2_2.sorted.pairsam.gz PLACseq_rep2_2.pairsam.gz
rm PLACseq_rep2_2.pairsam.gz
../../tools/bwa/bwa mem -SP5M -t 8 ../../genomes/hg19/BWAIndex/version0.6.0/genome.fa 4DNFICGRXFU6.fastq.gz 4DNFIJYF2QVF.fastq.gz | gzip -3 > PLACseq_rep2_3.sam.gz
pairtools parse --nproc-in 8 --cmd-in gunzip --cmd-out gzip -c ../../genomes/hg19/hg19.chrom.sizes --output-stats stats_PLACseq_rep2_3_map.txt -o PLACseq_rep2_3.pairsam.gz PLACseq_rep2_3.sam.gz
rm *.sam.gz
pairtools sort --cmd-in gunzip --cmd-out gzip --nproc-in 8 --nproc 8 --nproc-out 8 -o PLACseq_rep2_3.sorted.pairsam.gz PLACseq_rep2_3.pairsam.gz
rm PLACseq_rep2_3.pairsam.gz
pairtools merge --nproc-in 8 --cmd-in gunzip --cmd-out gzip --nproc 8 -o PLACseq_rep2.merged.sorted.pairsam.gz *.sorted.pairsam.gz
pairtools dedup --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-stats stats.PLACseq_rep2.dedup.txt -o PLACseq_rep2.dedup.pairsam.gz PLACseq_rep2.merged.sorted.pairsam.gz
pairtools select --cmd-in gunzip --cmd-out gzip --nproc-in 8 '((chrom1==chrom2) and (abs(pos1 - pos2) > 5000))' -o PLACseq_rep2.cis.pairsam.gz PLACseq_rep2.dedup.pairsam.gz
pairtools split --cmd-in gunzip --cmd-out gzip --nproc-in 8 --output-sam PLACseq_rep2.cis.sam.gz PLACseq_rep2.cis.pairsam.gz
samtools view -@ 7 -bh PLACseq_rep2.cis.sam.gz > PLACseq_rep2.cis.bam
samtools sort -@ 7 PLACseq_rep2.cis.bam > PLACseq_rep2.cis.sorted.bam
samtools index -@ 7 PLACseq_rep2.cis.sorted.bam
samtools view -b -P -h -@ 7 --region-file ../../MChIPC/binned_peaks.bed PLACseq_rep2.cis.sorted.bam > PLACseq_rep2.cis.filtered.bam
samtools index -@ 7 PLACseq_rep2.cis.filtered.bam
# merging replicates https://data.4dnucleome.org/experiment-set-replicates/4DNESWX1J3QU/
samtools merge -@ 7 -o PLACseq_merged.bam rep1/PLACseq_rep1.cis.filtered.bam rep2/PLACseq_rep2.cis.filtered.bam
samtools index -@ 7 PLACseq_merged.bam
bamCoverage -p 7 -b PLACseq_merged.bam -bs 250 -o PLACseq_merged.bw

# selecting individual viewpoints
# list_of_genes - file with viewpoint coordinates
MYC	chr8	128747750	128751500
MYB	chr6	135501750	135505750
TAL1	chr1	47693500	47698500
LINC00853	chr1	47644250	47646750
CMPK1	chr1	47798000	47801250
FOXD2	chr1	47901500	47904750
STIL	chr1	47777500	47780750
GATA1	chrX	48644500	48647750
HBG2	chr11	5274250	5276500
# making MChIPC bw profiles (250bp & 100bp resolution)
cat list_of_genes.txt | while read line; do  promoter=$(echo $line | cut -d " " -f 1); CHR=$(echo $line | cut -d " " -f 2); START=$(echo $line | cut -d " " -f 3); END=$(echo $line | cut -d " " -f 4); basename='MChIPC_merged'; samtools view -bh -P ../$basename.bam $CHR:$START-$END |samtools sort > $basename.$promoter.bam; samtools index $basename.$promoter.bam; bamCoverage -p 7 -b $basename.$promoter.bam -bs 250 -o $basename.$promoter.bw; bamCoverage -p 7 -b $basename.$promoter.bam -bs 100 -o $basename.$promoter.100.bw; done
# making PLACseq bg progiles (DpnII RE resolution)
cat list_of_genes.txt | while read line; do  promoter=$(echo $line | cut -d " " -f 1); CHR=$(echo $line | cut -d " " -f 2); START=$(echo $line | cut -d " " -f 3); END=$(echo $line | cut -d " " -f 4); basename='PLACseq_merged'; samtools view -bh -P ../$basename.bam $CHR:$START-$END |samtools sort > $basename.$promoter.bam; samtools index $basename.$promoter.bam; bedtools coverage -b $basename.$promoter.bam -a ../hg19.DpnII_digest.bed -sorted -counts -g ../../genomes/hg19/BWAIndex/hg19.txt > $basename.$promoter.bedgraph; done

# introducing CTCF and DNase profiles to the analysis: ENCFF000BWF - CTCF bw signal, ENCFF002CEL - peaks; ENCFF352SET - DNase bw signal, ENCFF621ZJY - peaks. Also H3K27: ENCFF010PHG


# selecting peaks of interest with R script 'creating_site_of_interest_files.R' 
# ploting heatmap to show higher sensitivity of MChIPC compared to PLACseq
computeMatrix reference-point --referencePoint center -S ../../../Nov2021_MCC/Chips/ENCFF352SET\ -\ DNase.bigWig ../../../Nov2021_MCC/Chips/ENCFF000BWF\ -\ CTCF.bigWig ../../ChIPs/ChIP_mean.bw ../../MChIPC_merged.bw ../../PLACseq_merged.bw -R CTCF_of_interest.bed DNase_of_interest.bed -p 3 -a 5000 -b 5000 -o matrix -bs 250 --missingDataAsZero 
plotHeatmap -m matrix --colorMap OrRd -o heatmap.pdf --zMax 1 100 10 50 25 --yMax 1 100 10 50 25 --yMin 0
# there are 4 CTCF sites on the profile-plot (one is behind DNase):
chr16	115656	116010	545
chr16	146853	147265	438
chr16	157055	157349	325
chr16	167868	168052	685
# and 4 DNase sites:
chr16	155055	155205	671
chr16	163555	163705	1045
chr16	170015	170165	1524
chr16	214895	215045	2558
	
#Note about ChIP-seq datasets used in analysis:
# 1) 267 bed.gz peak files & 267 corresponfing fold change bigWig files were downloaded directly from ENCODE (hg19)
# 2) 1 ChIP-seq datasets (CHD2) was available from ENCODE only for hg38 reference genome, these files were downloaded and lifted to hg19 with the following commands:
# CrossMap.py bigwig hg38ToHg19.over.chain.gz old.bw new.bw
# ./liftOver old.bed hg38ToHg19.over.chain.gz new.bed unmapped
# 3) 2 ChIP-seq datasets (CDK8 - 56379.bw and MED1) were downloaded from cistrome database and were also lifted to hg19 genome as described before
# 4) finally BRD4 dataset was downloaded from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2635249), here converted Wig to bigWig

# CTCF analysis
# everything is in "CTCF_analysis.R" script

# feature enrichment
# "feature_enrichment.R" + Homer (see the command in the R script)

# hierarcical clustering
# "clustering.R"

# random forest models
# "model_building.R" + Homer and deeptools (see the commands in the R script)

# analysis of loops in verified E-P pairs
# "comparison_with_verified_E-P.R"

