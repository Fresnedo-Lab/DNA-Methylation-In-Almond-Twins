## RNASeq Data Analysis with DESeq2

**Software:**
R v. 3.6.3

Load required packages into R.

```library(DESeq2)```

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
