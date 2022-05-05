# OTU/ZOTU pipeline (PEAR + UPARSE/UCHIME)

We compared our original OTU based approach with two ASV based approaches (UNOISE3 and dada2). For dada2 see below

## Pre-processing of raw amplicon data

We started with 244 paired end fastq files, including 96 that are technical replicates (same DNA extract, separately PCR'ed and sequenced).  
Sequencing was done on an Illumina Miseq. We sequenced the V1-V3 region (using the "pA" and "BKL1" primers) of the 16S gene using 2 x 300 bp reads

## merging forward and reverse reads 

```bash
for a in *R1.fastq; do echo ${a}; STEM=$(basename "${a}" R1.fastq); echo ${STEM}; pear -f ${STEM}R1.fastq -r ${STEM}R2.fastq -o ./${STEM}pear -v 10 -m 550 -n 300 -u 0 -j 10 | tee merging_output.txt; done
```

parameters used: minimum overlap (-v)=10, maximum length of assembled seqs (-m)=550, minimum assembled length (-n)=300, maximum proportion of uncalled bases in a read (-u)=0, and number of threads (-j)=10

### labeling the reads

```bash
mkdir ./labeled
out=./labeled

for a in *.fastq; do echo ${a}; STEM=$(basename "${a}" .fastq); echo ${STEM}; sed "-es/^@\(M0\)\(.*\)/@\1\2;barcodelabel=${STEM};/" < ${STEM}.fastq > $out/${STEM}.fastq; done
```

### Quality control

```bash
mkdir maxee2
out=./maxee2

for a in *.fastq;do STEM=$(basename "${a}" .fastq); usearch -fastq_filter ${a} -fastq_maxee 2 -fastq_minlen 450 -fastaout $out/${STEM}.fasta -log $out/${STEM}_maxee2_usearch.log; done
```

### dereplicate reads and sort by size

```bash
cat *.fasta > merged.fasta
vsearch --derep_fulllength merged.fasta --output merged_unique.fasta --sizeout --threads 10
vsearch --sortbysize merged_unique.fasta --output merged_unique_sorted.fasta --minsize 2 --threads 10
```

### cluster OTUs or denoise reads

```bash
usearch -cluster_otus merged_unique_sorted.fasta -otus 16S_all_maxee2_otu97.fasta -relabel OTU_ -log clustered_16S_all_maxee2_otu97.log

#or denoise into ASVs (UNOISE3):

usearch -unoise3 merged_unique_sorted.fasta -zotus 16S_all_maxee2_asv.fasta -tabbedout unoise3_16S_all_maxee2.txt
#alternatively you can use DADA2 instead of UNOISE (we started but did not yet finish comparing UPARSE and UNOISE output to DADA2 output)
```

### Map sample reads to OTUs/ASVs to generate OTU table

```bash
vsearch --usearch_global merged.fasta --db 16S_all_maxee2_otu97.fasta --id 0.97 --notrunclabels --strand plus --otutabout otutab_16S_all_maxee2_97.txt --threads 8

#or UNOISE3 option:

usearch -otutab merged.fasta -zotus 16S_all_maxee2_asv.fasta -otutabout otutab_16S_all_maxee2_asv.txt -notmatched unmatched_16S_all_maxee2_asv.fa -dbmatched ASVs_w_sizes_16S_all_maxee2.fa -log map_reads_to_ASVs_16S_all_maxee2.log -sizeout

# dada2

```bash
# DADA2 sequence processing

# Tutorial: https://benjjneb.github.io/dada2/tutorial.html
# https://benjjneb.github.io/dada2/bigdata_paired.html
# https://benjjneb.github.io/dada2/bigdata.html


rm(list=ls())





############################ test all



#####################

# Filtering script: #

#####################



# File parsing

library(dada2); packageVersion("dada2")



# File parsing




####################
### RUN 1: pool2 ###
####################


path <- "/home/slambrecht/Amplicon/FINAL_RAW_FASTQ_FILES/microbian_16S_final_raw_files/pool2_raw"

list.files(path)

pathF <- path # CHANGE ME to the directory containing your demultiplexed forward-read fastqs

pathR <- path # CHANGE ME ...

filtpathF <- file.path(pathF, "filteredF") # Filtered forward files go into the pathF/filtered/ subdirectory

filtpathR <- file.path(pathR, "filteredR") # ...

fastqFs <- sort(list.files(pathF, pattern="_R1.fastq.gz"))

fastqRs <- sort(list.files(pathR, pattern="_R2.fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")



# # Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS

# default

# filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),

#               rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),

#               truncLen=c(240,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,

#               compress=TRUE, verbose=TRUE, multithread=TRUE)


filtered_out_run1 <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),

                                   rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),

                                   maxEE=c(2,2),truncQ = 2, maxN=0, rm.phix=TRUE,truncLen=c(297,262),

                                   trimLeft = c(20,21),compress=TRUE, verbose=TRUE, multithread=TRUE)



###########################

# sample inference script: #

###########################



library(dada2); packageVersion("dada2")

# File parsing

filtpathF <- paste(path,"/filteredF",sep="") # CHANGE ME to the directory containing your filtered forward fastqs

filtpathR <- paste(path,"/filteredR",sep="") # CHANGE ME ...

filtFs <- list.files(filtpathF, pattern="fastq", full.names = TRUE)

filtRs <- list.files(filtpathR, pattern="fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_R1.fastq.gz"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

sample.namesR <- sapply(strsplit(basename(filtRs), "_R2.fastq.gz"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names

names(filtRs) <- sample.names

set.seed(100)

# Learn forward error rates

errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)

# Learn reverse error rates

errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Sample inference and merger of paired-end reads

mergers <- vector("list", length(sample.names))

names(mergers) <- sample.names

for(sam in sample.names) {

  cat("Processing:", sam, "\n")

  derepF <- derepFastq(filtFs[[sam]])

  ddF <- dada(derepF, err=errF, multithread=TRUE)

  derepR <- derepFastq(filtRs[[sam]])

  ddR <- dada(derepR, err=errR, multithread=TRUE)

  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose=TRUE)

  mergers[[sam]] <- merger

}

rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras

seqtab <- makeSequenceTable(mergers)

saveRDS(seqtab, paste(path,"_seqtab.rds")) # CHANGE ME to where you want sequence table saved



#### COUNT OVERVIEW RUN 1 ###########################

# set function

getN <- function(x) sum(getUniques(x))

# making count overview table

row.names=sample.names
dada2_input=filtered_out_run1[,1]
filtered=filtered_out_run1[,2]
merged=sapply(mergers, getN)
final_perc_reads_retained=round(merged/filtered_out_run1[,1]*100, 1)

#dada_f=sapply(ddF, getN) # doesn't work
#dada_r=sapply(ddR, getN) # doesn't work

summary_pool2 <- data.frame(row.names, dada2_input, filtered, merged, final_perc_reads_retained)

summary_pool2

write.table(summary_pool2, "read-count-tracking_16S_pool2_maxee22_297_262.tsv", quote=FALSE, sep="\t", col.names=NA)




################################
### RUN 2: resequenced pool2 ###
################################



path <- "/home/slambrecht/Amplicon/FINAL_RAW_FASTQ_FILES/microbian_16S_final_raw_files/resequenced_pool2_raw_but_ready"

list.files(path)

pathF <- path # CHANGE ME to the directory containing your demultiplexed forward-read fastqs

pathR <- path # CHANGE ME ...

filtpathF <- file.path(pathF, "filteredF") # Filtered forward files go into the pathF/filtered/ subdirectory

filtpathR <- file.path(pathR, "filteredR") # ...

fastqFs <- sort(list.files(pathF, pattern="_R1.fastq"))

fastqRs <- sort(list.files(pathR, pattern="_R2.fastq"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")



# # Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS

# default

# filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),

#               rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),

#               truncLen=c(240,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,

#               compress=TRUE, verbose=TRUE, multithread=TRUE)



# filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),

#               rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),

#               maxEE=2, truncQ=20, maxN=0, rm.phix=TRUE,

#               compress=TRUE, verbose=TRUE, multithread=TRUE)



filtered_out_run2 <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),

                                   rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),

                                   maxEE=c(2,2),truncQ = 2, maxN=0, rm.phix=TRUE,truncLen=c(297,262),

                                   trimLeft = c(20,21),compress=TRUE, verbose=TRUE, multithread=TRUE)


#####################

# sample inference script: #

#####################



library(dada2); packageVersion("dada2")

# File parsing

filtpathF <- paste(path,"/filteredF",sep="") # CHANGE ME to the directory containing your filtered forward fastqs

filtpathR <- paste(path,"/filteredR",sep="") # CHANGE ME ...

filtFs <- list.files(filtpathF, pattern="fastq", full.names = TRUE)

filtRs <- list.files(filtpathR, pattern="fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_R1.fastq"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

sample.namesR <- sapply(strsplit(basename(filtRs), "_R2.fastq"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names

names(filtRs) <- sample.names

set.seed(100)

# Learn forward error rates

errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)

# Learn reverse error rates

errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Sample inference and merger of paired-end reads

mergers <- vector("list", length(sample.names))

names(mergers) <- sample.names

for(sam in sample.names) {

  cat("Processing:", sam, "\n")

  derepF <- derepFastq(filtFs[[sam]])

  ddF <- dada(derepF, err=errF, multithread=TRUE)

  derepR <- derepFastq(filtRs[[sam]])

  ddR <- dada(derepR, err=errR, multithread=TRUE)

  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose=TRUE)

  mergers[[sam]] <- merger

}

rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras

seqtab <- makeSequenceTable(mergers)

saveRDS(seqtab, paste(path,"_seqtab.rds")) # CHANGE ME to where you want sequence table saved



#### COUNT OVERVIEW RUN 2 ###########################

# set function

getN <- function(x) sum(getUniques(x))

# making count overview table

row.names=sample.names
dada2_input=filtered_out_run2[,1]
filtered=filtered_out_run2[,2]
merged=sapply(mergers, getN)
final_perc_reads_retained=round(merged/filtered_out_run2[,1]*100, 1)

#dada_f=sapply(ddF, getN) # doesn't work
#dada_r=sapply(ddR, getN) # doesn't work

summary_resequenced_pool2 <- data.frame(row.names, dada2_input, filtered, merged, final_perc_reads_retained)

summary_resequenced_pool2

write.table(summary_resequenced_pool2, "read-count-tracking_16S_resequenced_pool2_maxee22_297_262.tsv", quote=FALSE, sep="\t", col.names=NA)






####################
###   plate 3  #####
####################





path <- "/home/slambrecht/Amplicon/FINAL_RAW_FASTQ_FILES/microbian_16S_final_raw_files/plate3_raw_but_ready_correctnames"

list.files(path)

pathF <- path # CHANGE ME to the directory containing your demultiplexed forward-read fastqs

pathR <- path # CHANGE ME ...

filtpathF <- file.path(pathF, "filteredF") # Filtered forward files go into the pathF/filtered/ subdirectory

filtpathR <- file.path(pathR, "filteredR") # ...

fastqFs <- sort(list.files(pathF, pattern="_R1.fastq"))

fastqRs <- sort(list.files(pathR, pattern="_R2.fastq"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")



# # Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS

# default

# filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),

#               rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),

#               truncLen=c(240,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,

#               compress=TRUE, verbose=TRUE, multithread=TRUE)



# filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),

#               rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),

#               maxEE=2, truncQ=20, maxN=0, rm.phix=TRUE,

#               compress=TRUE, verbose=TRUE, multithread=TRUE)

#



filtered_out_run3 <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),

                                   rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),

                                   maxEE=c(2,2),truncQ = 2, maxN=0, rm.phix=TRUE,truncLen=c(297,262),

                                   trimLeft = c(20,21),compress=TRUE, verbose=TRUE, multithread=TRUE)





#####################

# sample inference script: #

#####################



library(dada2); packageVersion("dada2")

# File parsing

filtpathF <- paste(path,"/filteredF",sep="") # CHANGE ME to the directory containing your filtered forward fastqs

filtpathR <- paste(path,"/filteredR",sep="") # CHANGE ME ...

filtFs <- list.files(filtpathF, pattern="fastq", full.names = TRUE)

filtRs <- list.files(filtpathR, pattern="fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_R1.fastq"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

sample.namesR <- sapply(strsplit(basename(filtRs), "_R2.fastq"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names

names(filtRs) <- sample.names

set.seed(100)

# Learn forward error rates

errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)

# Learn reverse error rates

errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Sample inference and merger of paired-end reads

mergers <- vector("list", length(sample.names))

names(mergers) <- sample.names

for(sam in sample.names) {

  cat("Processing:", sam, "\n")

  derepF <- derepFastq(filtFs[[sam]])

  ddF <- dada(derepF, err=errF, multithread=TRUE)

  derepR <- derepFastq(filtRs[[sam]])

  ddR <- dada(derepR, err=errR, multithread=TRUE)

  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose=TRUE)

  mergers[[sam]] <- merger

}

rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras

seqtab <- makeSequenceTable(mergers)

saveRDS(seqtab, paste(path,"_seqtab.rds")) # CHANGE ME to where you want sequence table saved




#### COUNT OVERVIEW RUN 3 ###########################

# set function

getN <- function(x) sum(getUniques(x))

# making count overview table

row.names=sample.names
dada2_input=filtered_out_run3[,1]
filtered=filtered_out_run3[,2]
merged=sapply(mergers, getN)
final_perc_reads_retained=round(merged/filtered_out_run3[,1]*100, 1)

#dada_f=sapply(ddF, getN) # doesn't work
#dada_r=sapply(ddR, getN) # doesn't work

summary_plate3 <- data.frame(row.names, dada2_input, filtered, merged, final_perc_reads_retained)

summary_plate3

write.table(summary_plate3, "read-count-tracking_16S_plate3_maxee22_297_262.tsv", quote=FALSE, sep="\t", col.names=NA)







####################
### 2019 samples ###
####################



path <- "/home/slambrecht/Amplicon/FINAL_RAW_FASTQ_FILES/microbian_16S_final_raw_files/2019samples_16S"

list.files(path)

pathF <- path # CHANGE ME to the directory containing your demultiplexed forward-read fastqs

pathR <- path # CHANGE ME ...

filtpathF <- file.path(pathF, "filteredF") # Filtered forward files go into the pathF/filtered/ subdirectory

filtpathR <- file.path(pathR, "filteredR") # ...

fastqFs <- sort(list.files(pathF, pattern="_R1.fastq"))

fastqRs <- sort(list.files(pathR, pattern="_R2.fastq"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")



# # Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS

# default

# filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),

#               rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),

#               truncLen=c(240,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,

#               compress=TRUE, verbose=TRUE, multithread=TRUE)



# filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),

#               rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),

#               maxEE=2, truncQ=20, maxN=0, rm.phix=TRUE,

#               compress=TRUE, verbose=TRUE, multithread=TRUE)


filtered_out_run4 <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),

                                   rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),

                                   maxEE=c(2,2),truncQ = 2, maxN=0, rm.phix=TRUE,truncLen=c(297,262),

                                   trimLeft = c(20,21),compress=TRUE, verbose=TRUE, multithread=TRUE)



#####################

# sample inference script: #

#####################



library(dada2); packageVersion("dada2")

# File parsing

filtpathF <- paste(path,"/filteredF",sep="") # CHANGE ME to the directory containing your filtered forward fastqs

filtpathR <- paste(path,"/filteredR",sep="") # CHANGE ME ...

filtFs <- list.files(filtpathF, pattern="fastq", full.names = TRUE)

filtRs <- list.files(filtpathR, pattern="fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_R1.fastq"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

sample.namesR <- sapply(strsplit(basename(filtRs), "_R2.fastq"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names

names(filtRs) <- sample.names

set.seed(100)

# Learn forward error rates

errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)

# Learn reverse error rates

errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Sample inference and merger of paired-end reads

mergers <- vector("list", length(sample.names))

names(mergers) <- sample.names

for(sam in sample.names) {

  cat("Processing:", sam, "\n")

  derepF <- derepFastq(filtFs[[sam]])

  ddF <- dada(derepF, err=errF, multithread=TRUE)

  derepR <- derepFastq(filtRs[[sam]])

  ddR <- dada(derepR, err=errR, multithread=TRUE)

  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose=TRUE)

  mergers[[sam]] <- merger

}

rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras

seqtab <- makeSequenceTable(mergers)

saveRDS(seqtab, paste(path,"_seqtab.rds")) # CHANGE ME to where you want sequence table saved





#### COUNT OVERVIEW RUN 4 ###########################

# set function

getN <- function(x) sum(getUniques(x))

# making count overview table

row.names=sample.names
dada2_input=filtered_out_run4[,1]
filtered=filtered_out_run4[,2]
merged=sapply(mergers, getN)
final_perc_reads_retained=round(merged/filtered_out_run4[,1]*100, 1)

#dada_f=sapply(ddF, getN) # doesn't work
#dada_r=sapply(ddR, getN) # doesn't work

summary_2019samples_16S <- data.frame(row.names, dada2_input, filtered, merged, final_perc_reads_retained)

summary_2019samples_16S

write.table(summary_2019samples_16S, "read-count-tracking_16S_2019samples_maxee22_297_262.tsv", quote=FALSE, sep="\t", col.names=NA)






#####################

# chimera/taxonomy script: #

#####################



library(dada2); packageVersion("dada2")

# Merge multiple runs (if necessary)

st1 <- readRDS("/home/slambrecht/Amplicon/FINAL_RAW_FASTQ_FILES/microbian_16S_final_raw_files/pool2_raw_seqtab.rds")

st2 <- readRDS("/home/slambrecht/Amplicon/FINAL_RAW_FASTQ_FILES/microbian_16S_final_raw_files/resequenced_pool2_raw_but_ready_seqtab.rds")

st3 <- readRDS("/home/slambrecht/Amplicon/FINAL_RAW_FASTQ_FILES/microbian_16S_final_raw_files/plate3_raw_but_ready_correctnames_seqtab.rds")

st4 <- readRDS("/home/slambrecht/Amplicon/FINAL_RAW_FASTQ_FILES/microbian_16S_final_raw_files/2019samples_16S_seqtab.rds")

st.all <- mergeSequenceTables(st1,st2,st3, st4) #st2,

# Remove chimeras

seqtab <- removeBimeraDenovo(st.all, method="consensus", verbose=T, multithread=TRUE)

getwd()



#########################################

#Track counts after chimera removal#

#########################################


getN <- function(x) sum(getUniques(x))
chimera_count <- cbind(rowSums(st.all), rowSums(seqtab))
colnames(chimera_count) <- c("tabled", "nonchim")
#rownames(chimera_count) <- sample.names
head(chimera_count)

#####


write.csv(seqtab,"dada2_multi16S_maxee22_trunclen297-262_all.csv")


# Assign taxonomy

tax <- assignTaxonomy(seqtab, "/home/refsets/ribosomal_DB/combinations_for_dada2/silva/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

# Write to disk

saveRDS(seqtab, "/home/slambrecht/Amplicon/FINAL_RAW_FASTQ_FILES/microbian_16S_final_raw_files/seqtab_final.rds") # CHANGE ME to where you want sequence table saved

saveRDS(tax, "/home/slambrecht/Amplicon/FINAL_RAW_FASTQ_FILES/microbian_16S_final_raw_files/tax_final.rds") # CHANGE ME ...



# OTU/ASV table

# replace colnames by OTU_

ncol(seqtab)

otu.vec <- as.vector(NULL)

for (otu in c(1:ncol(seqtab))){

  otu.count <- paste("OTU_",otu, sep="")

  otu.vec <- c(otu.vec,otu.count)

}



length(otu.vec)



seqtab.nochim.otu <- seqtab

colnames(seqtab.nochim.otu) <- otu.vec

seqtab.nochim.otu[c(1:10), c(1:10)]



write.csv(seqtab.nochim.otu, "microbian_16S_dada2_otutab_maxee22_trunclen297-262_all.csv")



# create fasta file

fasta.out <- NULL

for (otu in c(1:length(otu.vec))){

  print(otu)

  fasta.intermediate<- paste(">OTU_",otu,"\n",colnames(seqtab)[otu],sep="")

  fasta.out <-c(fasta.out,fasta.intermediate)

}



head(fasta.out)

write.table(fasta.out,"microbian_16S_dada2_otutab_maxee22_trunclen297_262_all.fasta", row.names = FALSE, quote = FALSE, sep="\t",col.names = FALSE)


##Transpose otu table


library(data.table)


##Transpose otu table


library(data.table)

data_almost <- read.csv("microbian_16S_dada2_otutab_maxee22_trunclen297-262_all.csv", header = T, sep = ",", stringsAsFactors = F)

data_almost <- read.csv("microbian_16S_dada2_otutab_maxee22_trunclen297_262_final.csv", header = T, sep = ",", stringsAsFactors = F)

tax_mothur <- read.csv("16S_dada2_262.csv", header = T, sep = ";", stringsAsFactors = F)


View(head(data_almost,100))

View(head(data_almost2,100))

View(head(tax_mothur,100))

dim(data_almost2)

dim(tax_mothur)

otutab_16S_maxee22_297_262_incl_tax <- cbind(data_almost2,tax_mothur)

View(head(otutab_16S_maxee22_297_262_incl_tax,100))

write.csv(otutab_16S_maxee22_297_262_incl_tax, "microbian_16S_dada2_otutab_maxee22_trunclen297_262_incl_tax.csv")



#transpose data frame
df_t <- transpose(data_almost)

#redefine row and column names
rownames(df_t) <- colnames(data_almost)
colnames(df_t) <- rownames(data_almost)

names(df_t) <- df_t[1,]
df_t <- df_t[-1,]

rownames(df_t) <- rownames(df_t)
colnames(df_t) <- colnames(df_t)

df_t <- df_t[,order(colnames(df_t))]

df_t[c(1:10), c(1:10)]

data_almost2 <- data_almost[,-1]
# rownames(data_almost2) <- data_almost[,1]

# data_almost2[c(1:10), c(1:10)]

View(head(df_t,100))

# View(head(data_almost2,100))


dat <- as.data.frame(sapply(df_t, as.numeric)) #<- sapply is here

View(head(dat,100))

data_almost2 <- dat

samples_above_5000reads <- data_almost2[, colSums(data_almost2) > 5000]
samples_below_5000reads <- data_almost2[, colSums(data_almost2) < 5000]

samples_above_4000reads <- data_almost2[, colSums(data_almost2) > 4000]
samples_below_4000reads <- data_almost2[, colSums(data_almost2) < 4000]

samples_above_3000reads <- data_almost2[, colSums(data_almost2) > 3000]
samples_below_3000reads <- data_almost2[, colSums(data_almost2) < 3000]


colSums(samples_below_5000reads) > 5397
colSums(samples_below_5000reads)

dim(data_almost2)
dim(samples_above_5000reads)
dim(samples_below_5000reads)

dim(data_almost2)
dim(samples_above_4000reads)
dim(samples_below_4000reads)

dim(data_almost2)
dim(samples_above_3000reads)
dim(samples_below_3000reads)


# View(head(samples_below_5000reads,100))


write.csv(dat, "microbian_16S_dada2_otutab_maxee22_trunclen297_262_final.csv")


# write.csv(df_t, "microbian_16S_dada2_otutab_maxee22_trunclen297_262_final.csv")

taxa.print <- tax # Removing sequence rownames for display only
rownames(taxa.print) <- NULL

# head(taxa.print)




```


