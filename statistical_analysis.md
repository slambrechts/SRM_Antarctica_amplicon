# statistical analysis

Insert description  

##  Alpha diversity

richness & Shannon and Simpson indices

calculate this on the matrices resulting from multiple rarefactions, with our smallest selected sample size (2683 reads?) as the maximum depth Collate results and average to obtain a single representative value for each sample.

```bash
insert code here
```

### Prepare OTU table for downstream analysis

normalize the non-rarefied OTU table using cumulative-sum scaling (CSS)
Calculate Bray-Curtis dissimilarity distance of OTU table using phyloseq or vegan

```bash
insert code here
```

### Calculate correlations between and investigate distribution of environmental variables

```bash
library(corrplot)

setwd(dir = "C:/Amplicon new/16S_all_maxee2_uparse97/distribution and correlations metadata/")

metadata <- read.csv("metadata_16S_numeric.csv", check.names = FALSE, sep = ";")

Correlation <- cor(metadata, use="complete.obs") # using complete.obs

corrplot(Correlation, method = "number") # Display the Pearson correlation coefficient (Pearson is the default method)
```

### Incorporate distance between sampling sites

```bash
insert code here
```
### Variation partitioning analysis

```bash
insert code here
```

### Indicator taxa analysis

```bash
insert code here
```
