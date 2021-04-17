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

```

### Calculate correlations between and investigate distribution of environmental variables

Correlations:

```bash
library(corrplot)

setwd(dir = "C:/Amplicon new/16S_all_maxee2_uparse97/distribution and correlations metadata/")

metadata <- read.csv("metadata_16S_numeric.csv", check.names = FALSE, sep = ";")

Correlation <- cor(metadata, use="complete.obs") # using complete.obs

corrplot(Correlation, method = "number") # Display the Pearson correlation coefficient (Pearson is the default method)

# Multicorrelation

vifstep(metadata[, c(22, 23, 32:44)],th=10)

vifstep(metadata_no_repl[, c(22, 23, 32, 34, 35, 37:44)],th=10) #removing pH_K_Cl (because less accurate - see Josef's email), and Soil_dry_weight_I

vifstep(metadata_no_repl[, c(22, 23, 32, 34, 35, 38:44)],th=10) #without N_NH4
```

Distribution:

```bash
library(ggpubr)
library(moments)

hist(metadata$pH_dest, col='steelblue')
hist(metadata$conductivity, col='steelblue')
hist(metadata$TOC, col='steelblue')
hist(metadata$Soil_dry_weight, col='steelblue')
hist(metadata$Soil_dry_way_II, col='steelblue')
hist(metadata$N.NH4, col='steelblue')
hist(metadata$N_NO3, col='steelblue')
hist(metadata$TN_g_kg, col='steelblue')
hist(metadata$TP_g_kg, col='steelblue')
hist(metadata$dem, col='steelblue')
hist(metadata$aspect, col='steelblue')
hist(metadata$slope, col='steelblue')

skewness(metadata$pH_dest, na.rm = TRUE)
skewness(metadata$Soil_dry_weight, na.rm = TRUE)
skewness(metadata$Soil_dry_way_II, na.rm = TRUE)
skewness(metadata$TOC, na.rm = TRUE)
skewness(metadata$N.NH4, na.rm = TRUE)
skewness(metadata$conductivity, na.rm = TRUE)
skewness(metadata$TN_g_kg, na.rm = TRUE)
skewness(metadata$TP_g_kg, na.rm = TRUE)
skewness(metadata$N_NO3, na.rm = TRUE)
skewness(metadata$P_PO4, na.rm = TRUE)
skewness(metadata$dem, na.rm = TRUE)
skewness(metadata$aspect, na.rm = TRUE)
skewness(metadata$slope, na.rm = TRUE)

kurtosis(metadata$pH_dest, na.rm = TRUE)
kurtosis(metadata$Soil_dry_weight, na.rm = TRUE)
kurtosis(metadata$Soil_dry_way_II, na.rm = TRUE)
kurtosis(metadata$TOC, na.rm = TRUE)
kurtosis(metadata$N.NH4, na.rm = TRUE)
kurtosis(metadata$conductivity, na.rm = TRUE)
kurtosis(metadata$TN_g_kg, na.rm = TRUE)
kurtosis(metadata$TP_g_kg, na.rm = TRUE)
kurtosis(metadata$N_NO3, na.rm = TRUE)
kurtosis(metadata$P_PO4, na.rm = TRUE)
kurtosis(metadata$dem, na.rm = TRUE)
kurtosis(metadata$aspect, na.rm = TRUE)
kurtosis(metadata$slope, na.rm = TRUE)


#################### Transformations performed ###########################
metadata$conductivity <- log10(metadata_MA$conductivity)
metadata$TOC <- log10(metadata_MA$TOC)
metadata$TP_g_kg <- log10(metadata_MA$TP_g_kg)
metadata$N_NO3 <- log10(metadata_MA$N_NO3)
metadata$N.NH4 <- log10(metadata_MA$N.NH4)
metadata$TN_g_kg <- log10(metadata_MA$TN_g_kg)
metadata$P_PO4 <- log10(metadata_MA$P_PO4)
metadata$Soil_dry_weight <- log10(max(metadata_MA$Soil_dry_weight+1) - metadata_MA$Soil_dry_weight)
metadata$Soil_dry_way_II <- log10(max(metadata_MA$Soil_dry_way_II+1) - metadata_MA$Soil_dry_way_II)
metadata$slope <- log10(metadata_MA$slope)
metadata$aspect <- sqrt(metadata_MA$aspect)
## no transformation carried out for pH and dem

##Check distribution again
hist...
skewness...
kurtosis...
...
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
NbClust(data = metadata$TOC, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 5, method = "kmeans")

metadata_scaled <- cbind(metadata)

metadata_scaled$TOC_scaled <- scale(metadata$TOC)

metadata_scaled$conductivity <- scale(metadata$conductivity)
metadata_scaled$TP_g_kg <- scale(metadata$TP_g_kg)
metadata_scaled$N_NO3 <- scale(metadata$N_NO3)
metadata_scaled$N.NH4 <- scale(metadata$N.NH4)
metadata_scaled$TN_g_kg <- scale(metadata$TN_g_kg)
metadata_scaled$P_PO4 <- scale(metadata$P_PO4)
metadata_scaled$Soil_dry_weight <- scale(metadata$Soil_dry_weight)
metadata_scaled$Soil_dry_way_II <- scale(metadata$Soil_dry_way_II)
metadata_scaled$slope <- scale(metadata$slope)
metadata_scaled$pH_dest <- scale(metadata$pH_dest)
metadata_scaled$pH_K_Cl <- scale(metadata$pH_K_Cl)
metadata_scaled$dem <- scale(metadata$dem)

NbClust(data = metadata_scaled$TOC, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 12, method = "kmeans")
NbClust(data = metadata_scaled$TOC_scaled, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 12, method = "kmeans")

NbClust(data = metadata_scaled$TOC_scaled, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans")

set.seed(123)
km.res <- kmeans(metadata_scaled$TOC_scaled, 5, nstart = 25)
print(km.res)

metadata_scaled_new <- cbind(metadata_scaled, TOC_cat = km.res$cluster)

TOC <- subset(metadata_scaled_new, select=c(TOC,TOC_scaled,TOC_cat))


#rename the clusters in TOC_cat

metadata_scaled_new$TOC_cat[metadata_scaled_new$TOC_cat == 2] = "9-11"
metadata_scaled_new$TOC_cat[metadata_scaled_new$TOC_cat == 4] = "4-6"
metadata_scaled_new$TOC_cat[metadata_scaled_new$TOC_cat == 1] = "2.1-3.5"
metadata_scaled_new$TOC_cat[metadata_scaled_new$TOC_cat == 3] = "1-2"
metadata_scaled_new$TOC_cat[metadata_scaled_new$TOC_cat == 5] = "< 1"
```
