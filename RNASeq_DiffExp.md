## RNASeq Data Analysis with DESeq2

**Software:**
R v. 3.6.3

Load required packages into R.

```library(DESeq2)```
```library(ggplot2)```

Using the gene counts text files generated for each Stukey library in the RNASeq_Alignment, create a master text file containing a column for geneID and a column for each library containing counts for each cooresponding gene. Read this text file into R.

```StukeyCts <- read.csv("gene_counts_all.csv", row.names = "gene_id")```

Upload an annotation file containing information on the study design. Headers for this study include "SampleName" containing the individual library names, "Condition" indicating whether the particular library comes from a Stukey individual with or without budfailure (i.e. "no-BF" or "BF"), and "Twinpair" indicating which twin pair the particular library is from.

```StukeyAnno <- read.csv("design.csv", header=TRUE, sep=",", row.names = "SampleName")```

In the annotation file, make the "Twinpair" variable a factor.

```StukeyAnno$TwinPair <- as.factor(StukeyAnno$TwinPair)```

Create a DESeqDataSet object using the DESeqDataSetFromMatrix method.

```dds <- DESeqDataSetFromMatrix(countData = StukeyCts, colData = StukeyAnno, design = ~Condition+TwinPair)```

Filter the data to keep only the gene IDs that have counts greater than 10.

```keep <- rowSums(counts(dds)) >= 10```
```dds <- dds[keep,]```

Create a "run" variable that indicates which original sample each technical replicate comes from. In this case, we are simply labeling them 1, 2, 3, or 4 and classifying the Stukey libraries using these values.

```dds$run <- paste0("run",c(1,1,2,2,3,3,4,4))```

Collapse the technical replicates for each Stukey individual.

```ddsColl <- collapseReplicates(dds, dds$run)```

Perform differential expression analysis on the collapsed technical replicate DESeqDataSet object using the DESeq method.

```ddsColl <- DESeq(ddsColl)```

Prepare a DESeqResults object using the results method and a contrast of the bud failure and no bud failure condition.

```res<- results(ddsColl, contrast=c("Condition", "BF", "no-BF"))```

Visualize a summary of the results from the differential expression analysis comparing gene expression in bud failure and no bud failure twins.

```summary(res)```

Add shrunken log2 fold changes and standard errors to the the DESeqDataSet using the shrinkage estimator "ashr" (adaptive shrinkage estimator) which creates a new results object.

```res<- lfcShrink(ddsColl, contrast=c("Condition", "BF", "no-BF"), type="ashr")```

Generate a plot displaying the relationship between the mean of normalized counts for each gene and the log2 fold change. Genes with a log2 fold change greater than or less than a specific cutoff value will be highlighted.

```DESeq2::plotMA(res, ylim=c(-2,3), cex = 0.7)```

Order the results object by p-value.

```resOrdered <- res[order(res$pvalue),]```

Write a csv of ordered results containing each gene passing the initial count threshold (see above), the basemean, the log2foldchange (comparing BF to no-BF), the standard error of the log2foldchange, the p-value, and the adjusted p-value.

```write.csv(as.data.frame(resOrdered), file="condition_BF_results.csv")```

Subset the results object and create a results object only containing those genes with a significant log2 fold change based on a defined adjust p-value cutoff.

```resSig <- subset(resOrdered, padj < 0.1)```

Write a csv containing only the genes passing the significance filter.

```write.csv(as.data.frame(resSig), file = "condition_BF_sigresults.csv")```

Perform a variance stabilizing transformation using the DESeqDataSet generated following collapse of the technical replicates. This transformation is being performed blind of the sample design to perform quality assurance on the data. This command produces a DESeqTransform object.

```vsd<- varianceStabilizingTransformation(ddsColl, blind=T)```

Create a datarame using plotPCA and the previously created DESeqTransform object. Use the top 500 genes based on row variance for the principal components.

```pcaData<- plotPCA(vsd, ntop=500, intgroup=c("Condition", "TwinPair"), returnData = T)```

Create an object containing the percent variance explained by the first two principal compoenents.

```percentVar<- round(100 * attr(pcaData, "percentVar"))```

Plot the PCA using ggplot2.

```ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=TwinPair)) + geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))```
