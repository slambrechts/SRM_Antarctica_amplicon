```bash
rm(list = ls())
```

# load packages

```bash
require(dplyr)
require(plyr)
require(RAM)
library(phyloseq)
library(vegan)
library(Biostrings) # To read fasta file
library(reshape2)
library(ape) # to read tree file
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggbiplot)
library(devtools)
library(corrplot)
```


# load otutable 
```bash
data_almost <- read.csv("otutab_16S_all_maxee2_97.csv", header = T, sep = ",", stringsAsFactors = F)
#26282 OTUs, 257 samples
```
# change wrong names in the otu table
```bash
oldnames = c("DV_SF2_CON_1_A_R", "DV_SF2_CON_2_A", "DV_SF2_CON_2_A_18", "DV_SF2_CON_2_B", "DV_SF2_CON_2_B_18", "DV_SF2_CON_2_C", "DV_SF2_CON_2_C_18", "DV_SF2_CON_3_A", "DV_SF2_CON_3_A_18", "DV_SF2_CON_3_B", "DV_SF2_CON_3_B_18", "DV_SF2_CON_3_C", "DV_SF2_CON_3_C_18", "OTC15_CON_A_19", "UT_OTC_N_A_18", "UT_OTC_N_B_18", "UT_OTC_N_C_18")
newnames = c("DV_SF2_CON_1_A_18_R", "DV_SF2_2_A", "DV_SF2_2_A_18", "DV_SF2_2_B", "DV_SF2_2_B_18", "DV_SF2_2_C", "DV_SF2_2_C_18", "DV_SF2_3_A", "DV_SF2_3_A_18", "DV_SF2_3_B", "DV_SF2_3_B_18", "DV_SF2_3_C", "DV_SF2_3_C_18", "OTC15_16_17_CON_A_19", "UT_OTC2_A_18", "UT_OTC2_B_18", "UT_OTC2_C_18")

data_almost <- data_almost %>% 
  rename_with(~ newnames[which(all_of(oldnames == .x))], .cols = all_of(oldnames))
```
# subset otutable

### remove chloroplasts, mitochondria, eukaryota, "unknown" otus
(no Archaea present in the starting otu table)
```bash
chloroplasts <- data_almost[grepl('Chloroplast',data_almost$taxonomy),]
data_almost2 <- data_almost[!grepl('Chloroplast',data_almost$taxonomy),]
#26229 OTUs, 256 samples

mithocondria <- data_almost2[grepl('Mitochondria',data_almost2$taxonomy),]
data_almost5 <- data_almost2[!grepl('Mitochondria',data_almost2$taxonomy),]
#26227 OTUs

#remove unknown
unknown <- data_almost5[grepl('k__unknown',data_almost5$taxonomy),]
data_almost6 <- data_almost5[!grepl('k__unknown',data_almost5$taxonomy),]
#26220 otus, 256 samples

#remove Eukaryota
eukaryota <- data_almost6[grepl('Eukaryota',data_almost6$taxonomy),]
data_almost3 <- data_almost6[!grepl('Eukaryota',data_almost6$taxonomy),]
#26216 otus, 256 samples
```

### remove OTCs and other odd samples
```bash
data_almost4 <- select(data_almost3, - c("UT_OTC2_A_18", "UT_OTC2_B_18", "UT_OTC2_C_18", "DV_SF_1_B_RR_18", "DV_SF_1_B_R_18", "DV_SF_2_B_18", "UT_1_18_R", "AU_7_18", "C2", "C3"))
data_almost4

data <- select(data_almost4, - contains(c("OTC6", "OTC5", "control", "MOCK")))
#26216 OTUs, 231 samples

head(data) #to check if they were successfully removed
```

### remove samples without all the bedrock data
```bash
data_almost5 <- select(data, - c("PB_N_2", "PB_N_2_18", "PB_N_2_R_18", "PB_N_5_18", "PB_N_5_N", "PB_N_4_18", "PB_N_4_18_plate3", "PB_N_4_N", "PT_1_B_19", "UT_5", "UT_5_18", "UT_8", "UT_8_18", "UT_OTC11_B", "UT_OTC11_B_18", "UT_OTC11_B_R", "UT_OTC11_CON_A_18", "UT_OTC11_C_18"))


row.names(data_almost5) <- data_almost5[, 1] # make the first column as the row names
dim(data_almost5)
#26216, 211 
data <- data_almost5[, -1]
dim(data)
#26216 otus, 212 samples
```
### order the columns by alphabetic order
```bash
data <- data[,order(colnames(data))]

head(data) # check column names 

grep("taxonomy", colnames(data)) # find the "taxonomy" column

data<-data[,c(1:170, 172:212, 171)] #place the taxonomy column at the end of the table
head(data) # check
dim(data)
#26216 otus, 212 samples
```
# load metadata table
```bash
metadata_original <- read.csv("metadata_16S_all.csv", header = T, sep = ",", stringsAsFactors = T)
```
### subset metadata table
The idea is now to subset the metadata table without the aquatic samples, store the colnames as a vector and use it to finally subset the otutable in order to have the exact same set of samples (in the same order)

# remove aquatic samples
```bash
metadata_aq <- metadata_original %>% filter(habitat != "aquatic") #filter aquatic samples
dim(metadata_aq) #check dimension
```
# store the aquatic sample names 
```bash
aquatic_samples <- metadata_original %>% filter(habitat == "aquatic") 
row.names(aquatic_samples) <- aquatic_samples[, 1]
aquatic_samples2 <- row.names(aquatic_samples)
```
### final substet of the otutable

# remove the aquatic samples
```bash 
data_aq <- data %>% select(- all_of(aquatic_samples2))
dim(data_aq)
#26216 otus, 197 samples
head(data_aq)
```

# remove otus with 0 or 1 reads in the whole dataset 
```bash 
data_aq <- data_aq[rowSums(data_aq[!names(data_aq) %in% "taxonomy"])>1,] #select all columns but the taxonomy one
dim(data_aq)
#24193 otus, 197 samples

rowSums(data_aq[!names(data_aq) %in% "taxonomy"]) == 1 #check if all singletones were removed
```

### final subset of the metadata table: only retain the same samples as in the otutable
```bash
samples_labels <- colnames(data_aq[, !names(data_aq) %in% "taxonomy"])
samples_labels2 <- paste0(samples_labels, collapse = '$|')

metadata_almost <- metadata_original%>% filter(grepl(samples_labels2,SampleID))
dim(metadata_almost)
metadata_almost$SampleID[!(metadata_almost$SampleID %in% samples_labels2)]
row.names(metadata_almost) <- metadata_almost[, 1]
metadata <- metadata_almost[,-1]
meta_SampleID <- rownames(metadata)
```

# check for mismatches in both tables 
We use the SampleID column and the just created variable "samples_labels"
```bash
samples_labels[!(samples_labels %in% meta_SampleID)]
meta_SampleID[!(meta_SampleID %in% samples_labels)]
```
# both tables with the same sample names order (alphabetically)
```bash
metadata <- metadata[order(rownames(metadata)),]
data_aq2 <- data_aq[ , order(names(data_aq[, !names(data_aq2) %in% "taxonomy"]))]
data_aq2$taxonomy <- data_aq$taxonomy #add the taxonomy column
```

### remove replicates in both tables
```bash
metadata$SampleID <- rownames(metadata) 
metadata$sum_reads <- colSums(data_aq2[, !names(data_aq2) %in% "taxonomy"])

metadata_no_repl <- metadata %>% 
  group_by(sample_ID) %>%         #sample_ID is a column with the same sample name for every replicates of one same sample
  top_n(1, sum_reads)             #only keep samples that have the highest number of reads out of all the replicates for that same sample
```

# final metadata table
```bash
metadata_no_repl <- as.data.frame(metadata_no_repl) 
rownames(metadata_no_repl) <- metadata_no_repl$SampleID #with the rows called with the samplename
```
# create variable with the same sample names as in the metatable + "taxonomy", so that only the right columns will be kept in the otutable
```bash
samples_no_repl <- c(rownames(metadata_no_repl), "taxonomy") 
```
# final otu table
```bash
data_no_repl <- data_aq2 %>% select(all_of(samples_no_repl))
dim(data_no_repl)
#24193 otus, 104 samples
```

# remove otus with 0 or 1 reads in the whole dataset 
```bash 
data_no_repl <- data_no_repl[rowSums(data_no_repl[!names(data_no_repl) %in% "taxonomy"])>1,]
dim(data_no_repl)
#21543 otus, 104 samples

rowSums(data_no_repl[!names(data_no_repl) %in% "taxonomy",]) == 1 #just to check
```
# create .csv files for both tables 
```bash
write.csv(data_no_repl, file = "otutab_16S_all_maxee2_97_no_replicates.csv")
write.csv(metadata_no_repl, file = "meta_no_replicates.csv")

```

### Rarefaction

# keep samples which minimum reads number is 1000, before rarefaction
```bash 
samples_above_1000reads <- data_no_repl[, colSums(data_no_repl[!names(data_no_repl) %in% "taxonomy"]) > 1000]
samples_below_1000reads <- data_no_repl[, colSums(data_no_repl[!names(data_no_repl) %in% "taxonomy"]) < 1000]
colSums(samples_below_1000reads) > 1000 
colSums(samples_below_1000reads)        
colSums(samples_above_1000reads[!names(samples_above_1000reads) %in% "taxonomy"])
colSums(samples_above_1000reads[!names(samples_above_1000reads) %in% "taxonomy"]) < 1000
dim(samples_above_1000reads)
#24193 otus, 192 samples
```
# keep samples which minimum reads number is 2683, before rarefaction
```bash
samples_above_2683reads <- data_no_repl[, colSums(data_no_repl[!names(data_no_repl) %in% "taxonomy"]) > 2683]
samples_below_2683reads <- data_no_repl[, colSums(data_no_repl[!names(data_no_repl) %in% "taxonomy"]) < 2683]
colSums(samples_below_2683reads) > 2683
colSums(samples_below_2683reads)
dim(samples_above_2683reads)
# 21543 otus, 99 samples
colSums(samples_above_2683reads[!names(samples_above_2683reads) %in% "taxonomy"]) < 2683
colSums(samples_above_2683reads[!names(samples_above_2683reads) %in% "taxonomy"])
```

# keep samples which minimum reads number is 5000, before rarefaction
```bash
samples_above_5397reads <- data_no_repl[, colSums(data_no_repl[!names(data_no_repl) %in% "taxonomy"]) > 5397]
samples_below_5397reads <- data_no_repl[, colSums(data_no_repl[!names(data_no_repl) %in% "taxonomy"]) < 5397]
colSums(samples_below_5397reads) > 5397
colSums(samples_below_5397reads)
dim(samples_above_5397reads)
# 21543 otus    99 samples
colSums(samples_above_5397reads[, !names(samples_above_5397reads) %in% "taxonomy"]) < 5397
colSums(samples_above_5397reads[, !names(samples_above_5397reads) %in% "taxonomy"])
```
# rarefaction 
```bash
data.l <- list(data = samples_above_5397reads)

data_rar <- as.data.frame(OTU.rarefy(data.l))

colnames(data_rar) <- colnames(samples_above_5397reads)
```
# CLR transformation
```bash
library(microbiome)

data_rar_clr <- as.data.frame(microbiome::transform(data_rar[, !names(data_rar) %in% "taxonomy"], "clr"))  
data_rar_clr$taxonomy <- data_rar$taxonomy

write.csv(data_rar_clr, "otutab_CLR_rarefied5397.csv")

data_rar_clr$taxonomy <- data_rar$taxonomy
```
