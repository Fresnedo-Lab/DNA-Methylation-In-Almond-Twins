## Create Heatmaps and Venn Diagrams of DMRs for Stukey Twinpairs

**Software:**
R v. 3.6.3

Load required packages into R.
```library(grid)                                        library(VennDiagram)                                   library(wesanderson)                              library(ComplexHeatmap)                           library(circlize)                                        library(cowplot)```

Read in the csv file generated when calling DMRs using DSS for each methylation context and twin pair.

```CG_TP1 <- read.csv("../../Smoothing CG_DMRs_Pair1_0.01_DSS.csv", header = TRUE)```
```CG_TP2 <- read.csv("../../Smoothing CG_DMRs_Pair2_0.01_DSS.csv", header = TRUE)```

Read in the csv file containing the DMRs associated with the same gene in the same proximitiy class between the two twinpairs in each methylation context.

```CG <- read.csv("ObsCG_All.csv", header = TRUE)```

Create variables of the count of DMRs in each twin pair and the number of shared DMRs using count.

```count.Stukey1.CG <- nrow(CG_TP1)```
```count.Stukey2.CG <- nrow(CG_TP2)```
```count.Shared.CG <- nrow(CG)```

Create a Venn diagram displaying the number of DMRs in the CG context for each twin pair with an overlap region showing the number of shared DMRs in this context.

```V1 <- draw.pairwise.venn(count.Stukey1.CG, count.Stukey2.CG, count.Shared.CG, category = c("Stukey1", "Stukey2"), lty = rep("blank", 1), fill = wes_palette("Royal1",2), alpha = rep(0.5, 2), cat.pos = c(350,10), cat.dist = rep(0.045, 2), cex = c(1,1,1), rotation.degree = 180, scaled = FALSE, cat.cex = c(1,1))```

**Repeat this process for both the CHH and CHG context creating objects V2 and V3**

Merge shared CG DMRs with mean methylation data from the CG DMR file for twinpair 1.

```Merged_CG <- merge(CG, CG_TP1[,c(2:4,7:8)], by.x = c("Scaffold_TP1", "Start_TP1", "End_TP1"), by.y = c("chr", "start", "end"))```

Change the name of the columns in the newly created dataframe.

```colnames(Merged_CG)[c(11,12)] <- c("PerMeth_1a_BF", "PerMeth_1b_noBF")```

Merge dataframe created above with mean methylation data from the CG DMR file for twinpair 1.

```Merged_CG <- merge(Merged_CG, CG_TP2[,c(2:4,7:8)], by.x = c("Scaffold_TP2", "Start_TP2", "End_TP2"), by.y = c("chr", "start", "end"))```

Change the name of the new columsn added to the dataframe.

```colnames(Merged_CG)[c(13,14)] <- c("PerMeth_2a_BF", "PerMeth_2b_noBF")```

Calculate percent methylation from the mean methylation values presented for each DMR for each twin.

```Merged_CG$PerMeth_1a_BF <- Merged_CG$PerMeth_1a_BF*100```
```Merged_CG$PerMeth_1b_noBF <- Merged_CG$PerMeth_1b_noBF*100```
```Merged_CG$PerMeth_2a_BF <- Merged_CG$PerMeth_2a_BF*100```
```Merged_CG$PerMeth_2b_noBF <- Merged_CG$PerMeth_2b_noBF*100```

Create a matrix using the merged dataframe and subset that matrix to contain only the percent methylation values for each DMR for each twin.

```Merged_Mat_CG <- as.matrix(Merged_CG[, c(1,11:14)])```
```MethMatCG <- Merged_Mat_CG[,c(2:5)]```

Convert the subsetted matrix to a numeric matrix and recreate the matrix.

```MethMatCG <- mapply(MethMatCG, FUN=as.numeric)```
```MethMatCG <- matrix(data=MethMatCG, ncol=4, nrow=82)```

Change the column and row names of the matrix.

```colnames(MethMatCG) <- c("Stukey1a", "Stukey1b", "Stukey2a", "Stukey2b")```
```rownames(MethMatCG) <- Merged_CG$GeneName```

Create a variable storing information to make the colors of the heatmap with a min value of 0, a midpoint value of 50, and a max value of 100.

```col_fun = colorRamp2(c(0, 50, 100), c("red", "white", "blue"))```

Create a plot generating a heatmap displaying percent methylation in the DMRs associated with the same gene (shared) between twinpairs. Percent methylation is shown for each individual twin.

```Ht1 <- Heatmap(MethMatCG, name = "% CG DNA Methylation", show_column_dend = FALSE, col = col_fun, column_order = c("Stukey1a", "Stukey1b", "Stukey2a", "Stukey2b"), show_row_names = FALSE, show_column_names = TRUE)```

**Repeat this process for both the CHH and CHG context creating objects Ht2 and Ht3**

Create a pdf image with heatmaps and Venn diagrams for each methylation context.

```Ht1 <- grid.grabExpr(draw(Ht1))```
```Ht2 <- grid.grabExpr(draw(Ht2))```
```Ht3 <- grid.grabExpr(draw(Ht3))```

```pdf("CombinedHeatmapVennDiagram.pdf")```
```ggdraw() +```
  ```draw_plot(Ht1, x = 0, y = 0.65, width = 0.6, height = 0.3) +```
  ```draw_plot(Ht2, x = 0, y = 0.35, width = 0.6, height = 0.3) +```
  ```draw_plot(Ht3, x = 0, y = 0.05, width = 0.6, height = 0.3) +```
  ```draw_plot(V1, x = 0.65, y = 0.65, width = 0.3, height = 0.3) +```
  ```draw_plot(V2, x = 0.65, y = 0.35, width = 0.3, height = 0.3) +```
  ```draw_plot(V3, x = 0.65, y = 0.05, width = 0.3, height = 0.3) +```
  ```draw_plot_label(label = c("a", "b", "c"), size = 15, x = c(0,0,0), y=c(1,0.7,0.4))```
```dev.off()```
