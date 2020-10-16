## Merging Methylation and Gene Expression

Create a visualization of percent methylation within the gene-associated DMRs merged with expression of those genes (log2 fold change) based on proximity (upstream, overlapping, and downstream) or based on methylation context (CG, CHG, and CHH). Heatmaps show the log2 fold change for the gene (comparing BF to no-BF) and percent methylation difference (comparing BF to no-BF) in each twin pair. The DMRs are subsetted either by their proximity to a gene or by the methylation context.

**Software:**
R v. 3.6.3

Load required packages into R.


**Merge Methylation and Expression Data Based on Proximity**

Read in the csv file containing results from the gene expression analysis with the gene ID, the basemean, the log2foldchange (comparing BF to no-BF), the standard error of the log2foldchange, the p-value, and the adjusted p-value.

```Gene_Exp <- read.csv("condition_BF_results.csv", header = TRUE)```
