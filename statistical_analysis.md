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

Correlations:

```bash
library(corrplot)

setwd(dir = "C:/Amplicon new/16S_all_maxee2_uparse97/distribution and correlations metadata/")

metadata <- read.csv("metadata_16S_numeric.csv", check.names = FALSE, sep = ";")

Correlation <- cor(metadata, use="complete.obs") # using complete.obs

corrplot(Correlation, method = "number") # Display the Pearson correlation coefficient (Pearson is the default method)
```

Distribution:

```bash
library(ggpubr)
library(moments)

as.numeric(x = metadata$dem) ##does not work work
as.numeric(x = metadata$aspect) ##does not work
as.numeric(x = metadata$slope) ##does not work

skewness(metadata$pH_dest, na.rm = TRUE)
skewness(metadata$Soil_dry_weight, na.rm = TRUE)
skewness(metadata$Soil_dry_way_II, na.rm = TRUE)
skewness(metadata$TOC, na.rm = TRUE)
skewness(metadata$NH4, na.rm = TRUE)
skewness(metadata$conductivity, na.rm = TRUE)
skewness(metadata$TN, na.rm = TRUE)
skewness(metadata$TP, na.rm = TRUE)
skewness(metadata$NO3, na.rm = TRUE)
skewness(metadata$PO4, na.rm = TRUE)

kurtosis(metadata$pH_dest, na.rm = TRUE)
kurtosis(metadata$Soil_dry_weight, na.rm = TRUE)
kurtosis(metadata$Soil_dry_way_II, na.rm = TRUE)
kurtosis(metadata$TOC, na.rm = TRUE)
kurtosis(metadata$NH4, na.rm = TRUE)
kurtosis(metadata$conductivity, na.rm = TRUE)
kurtosis(metadata$TN, na.rm = TRUE)
kurtosis(metadata$TP, na.rm = TRUE)
kurtosis(metadata$NO3, na.rm = TRUE)
kurtosis(metadata$PO4, na.rm = TRUE)

metadata$conductivity <- log10(metadata$conductivity)

hist(metadata$pH_dest, col='steelblue')
hist(metadata$pH_KCl, col='steelblue')
hist(metadata$conductivity, col='steelblue')
hist(metadata$TOC, col='steelblue')
hist(metadata$dry_weight, col='steelblue')
hist(metadata$dry_weight_2, col='steelblue')
hist(metadata$`N-NH4`, col='steelblue')
hist(metadata$N_NO3, col='steelblue')
hist(metadata$TN, col='steelblue')
hist(metadata$TP, col='steelblue')
hist(metadata$dry_weight_2, col='steelblue')
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
