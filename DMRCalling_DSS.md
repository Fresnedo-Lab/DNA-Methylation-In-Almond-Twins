## DMR calling with DSS

**Software:**
R v. 3.6.3

Load required packages into R.

```library(DSS)```
```library(bsseq)```
```library(ggplot2)```

Create a path variable with the path to the directory where the text files containing the methylation calls for each methylation context for each twin are stored.

```path <- "path/to/methylation/calls/for/each/context"```

Read in the methylation call file for the first methylation context (in this case CG) for each twin.

```Stukey1a_CG <- read.table(file.path(path, "CGStukey1a_DSS.txt"), header = FALSE)```

```Stukey1b_CG <- read.table(file.path(path, "CGStukey1b_DSS.txt"), header = FALSE)```

```Stukey2a_CG <- read.table(file.path(path, "CGStukey2a_DSS.txt"), header = FALSE)```

```Stukey2b_CG <- read.table(file.path(path, "CGStukey2b_DSS.txt"), header = FALSE)```

Change the column names in the data frames created above to work with DSS>

```colnames(Stukey1a_CG) <- c("chr", "pos", "N", "X")```
```colnames(Stukey1b_CG) <- c("chr", "pos", "N", "X")```
```colnames(Stukey2a_CG) <- c("chr", "pos", "N", "X")```
```colnames(Stukey2b_CG) <- c("chr", "pos", "N", "X")```

Create a bsseq object for Stukey twin pair 1 and Stukey twin pair 2.

```BSobj1 = makeBSseqData( list(Stukey1a_CG, Stukey1b_CG), c("1a","1b"))```

```BSobj2 = makeBSseqData( list(Stukey2a_CG, Stukey2b_CG), c("2a","2b"))```

Create a variable representing the number of cores available for analysis. If using the Ohio Supercomputer Center, one node represents 28 cores.

```mParam = MulticoreParam(workers=28, progressbar=TRUE)```

Using DSS, call differentially methylated loci (DML) comparing twins in twin pair 1 and twins in twin pair 2. This compares the bud failure twin to the no bud failure twin in each twin pair.

```dmlTest1 = DMLtest(BSobj1, group1="1a", group2="1b", smoothing = TRUE, BPPARAM=mParam)```
```dmlTest2 = DMLtest(BSobj2, group1="2a", group2="2b", smoothing = TRUE, BPPARAM=mParam)```

Using DSS, call differentially methylated regions (DMRs) comparing twins in twin pair 1 and twins in twin pair 2. DMRs represent regions of the almond genome with significantly different levels of methylation (alpha = 0.01) when comparing the bud failure twin to the no bud failure twin in each twin pair.

```dmrsPair1_0.01 = callDMR(dmlTest1, p.threshold = 0.01)```

```dmrsPair2_0.01 = callDMR(dmlTest2, p.threshold = 0.01)```

Write the identified DMRs in each twin pair to a csv file to use for downstream analysis.

```write.csv(dmrsPair1_0.01, "Smoothing_CG_DMRs_Pair1_0.01_DSS.csv")```
```write.csv(dmrsPair2_0.01, "Smoothing_CG_DMRs_Pair2_0.01_DSS.csv")```

Create histograms of DMR lengths in basepairs using ggplot2 for DMRs in all methylation contexts for both twin pairs. First create a csv file with one column that is the length of the DMRs from each csv file generated above (i.e. Smoothing_CG_DMRs_Pair1_0.01_DSS.csv). Calculate the mean, median, and mode for the length of all DMRs.

Read in the csv file with one column for length

```DMR_Hist <- read.csv("DMR_Length_DSS_Smoothing.csv")```
```str(DMR_Hist)```
```summary(DMR_Hist)```
```pdf("Histogram DMR Length DSS.pdf")```
```ggplot(DMR_Hist, aes(x=Length)) + ```
  ```geom_histogram(binwidth = 10, color = "black", fill = "slateblue") +```
  ```xlab("DMR Length (bp)")+```
  ```ylab("Count")```
```dev.off()```

Creates a function to estimate the mode for DMR lengths

```estimate_mode <- function(x) {```
  ```d <- density(x)```
  ```d$x[which.max(d$y)]```
```}```
```estimate_mode(DMR_Hist$Length)```
