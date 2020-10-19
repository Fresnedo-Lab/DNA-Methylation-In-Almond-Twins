## Merging Methylation and Gene Expression

Create a visualization of percent methylation within the gene-associated DMRs merged with expression of those genes (log2 fold change) based on proximity (upstream, overlapping, and downstream) or based on methylation context (CG, CHG, and CHH). Heatmaps show the log2 fold change for the gene (comparing BF to no-BF) and percent methylation difference (comparing BF to no-BF) in each twin pair. The DMRs are subsetted either by their proximity to a gene or by the methylation context.

**Software:**
R v. 3.6.3

Load required packages into R.


Read in the csv file containing results from the gene expression analysis with the gene ID, the basemean, the log2foldchange (comparing BF to no-BF), the standard error of the log2foldchange, the p-value, and the adjusted p-value.

```Gene_Exp <- read.csv("condition_BF_results.csv", header = TRUE)```

Read in the files containing the percent methylation data for each DMR for each twin pair based on methylation context.

```CG_TP1 <- read.csv("../../Smoothing CG_DMRs_Pair1_0.01_DSS.csv", header = TRUE)```
```CG_TP2 <- read.csv("../../Smoothing CG_DMRs_Pair2_0.01_DSS.csv", header = TRUE)```
```CHG_TP1 <- read.csv("../../Smoothing CHG_DMRs_Pair1_0.01_DSS.csv", header = TRUE)```
```CHG_TP2 <- read.csv("../../Smoothing CHG_DMRs_Pair2_0.01_DSS.csv", header = TRUE)```
```CHH_TP1 <- read.csv("../../Smoothing CHH_DMRs_Pair1_0.01_DSS.csv", header = TRUE)```
```CHH_TP2 <- read.csv("../../Smoothing CHH_DMRs_Pair2_0.01_DSS.csv", header = TRUE)```

**Merge Methylation and Expression Data Based on Proximity**

Read in the files containing the genomic coordinates for DMRs from twin pair 1 and twin pair 2 located either upstream (w/in 2,000 bp), downstream (w/in 2,000 bp), or overlapping the same gene, as well as genomic coordinates for the gene.

```Upstream <- read.csv("Upstream_Subset.csv", header = TRUE)```
```Overlapping <- read.csv("Overlapping_Subset.csv", header = TRUE)```
```Downstream <- read.csv("Downstream_Subset.csv", header = TRUE)```

Merge upstream DMRs with percent methylation Data from DMR files for each twin pair for each methylation context.

Twin pair 1 - CG context
``Merged_Upstream_CG <- merge(Upstream, CG_TP1[,c(2:4,7:8)], by.x = c("Chr_1", "start_1", "end_1"), by.y = c("chr", "start", "end"))``
``colnames(Merged_Upstream_CG)[c(18,19)] <- c("PerMeth_1a_BF", "PerMeth_1b_noBF")``

Twin pair 2 - CG context
```Merged_Upstream_CG <- merge(Merged_Upstream_CG, CG_TP2[,c(2:4,7:8)], by.x = c("Chr_2", "start_2", "end_2"), by.y = c("chr", "start", "end"))```
```colnames(Merged_Upstream_CG)[c(20,21)] <- c("PerMeth_2a_BF", "PerMeth_2b_noBF")```

Twin pair 1 - CHG context
```Merged_Upstream_CHG <- merge(Upstream, CHG_TP1[,c(2:4,7:8)], by.x = c("Chr_1", "start_1", "end_1"), by.y = c("chr", "start", "end"))```
```colnames(Merged_Upstream_CHG)[c(18,19)] <- c("PerMeth_1a_BF", "PerMeth_1b_noBF")```

Twin pair 2 - CHG context
```Merged_Upstream_CHG <- merge(Merged_Upstream_CHG, CHG_TP2[,c(2:4,7:8)], by.x = c("Chr_2", "start_2", "end_2"), by.y = c("chr", "start", "end"))```
```colnames(Merged_Upstream_CHG)[c(20,21)] <- c("PerMeth_2a_BF", "PerMeth_2b_noBF")```

Twin pair 1 - CHH context
```Merged_Upstream_CHH <- merge(Upstream, CHH_TP1[,c(2:4,7:8)], by.x = c("Chr_1", "start_1", "end_1"), by.y = c("chr", "start", "end"))```
```colnames(Merged_Upstream_CHH)[c(18,19)] <- c("PerMeth_1a_BF", "PerMeth_1b_noBF")```

Twin pair 2 - CHH context
```Merged_Upstream_CHH <- merge(Merged_Upstream_CHH, CHH_TP2[,c(2:4,7:8)], by.x = c("Chr_2", "start_2", "end_2"), by.y = c("chr", "start", "end"))```
```colnames(Merged_Upstream_CHH)[c(20,21)] <- c("PerMeth_2a_BF", "PerMeth_2b_noBF")```

Merge percent methylation infomration for all methylation contexts upstream of a gene
```Merged_Upstream <- rbind(Merged_Upstream_CG, Merged_Upstream_CHG, Merged_Upstream_CHH)```

**Repeat the above merging for both the downstream and overlapping proximitiy groups**
