# Pre-processing of raw amplicon data

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
