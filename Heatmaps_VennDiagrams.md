## Create Heatmaps and Venn Diagrams of DMRs for Stukey Twinpairs

**Software:**
R v. 3.6.3

Load required packages into R.
```library(grid)```
```library(VennDiagram)```
```library(wesanderson)```

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

Create a pdf image with heatmaps and Venn diagrams for each methylation context.

```pdf("CombinedHeatmapVennDiagram2.pdf")```
```ggdraw() +```
  ```draw_plot(Ht1, x = 0, y = 0.65, width = 0.6, height = 0.3) +```
  ```draw_plot(Ht2, x = 0, y = 0.35, width = 0.6, height = 0.3) +```
  ```draw_plot(Ht3, x = 0, y = 0.05, width = 0.6, height = 0.3) +```
  ```draw_plot(V1, x = 0.65, y = 0.65, width = 0.3, height = 0.3) +```
  ```draw_plot(V2, x = 0.65, y = 0.35, width = 0.3, height = 0.3) +```
  ```draw_plot(V3, x = 0.65, y = 0.05, width = 0.3, height = 0.3) +```
  ```draw_plot_label(label = c("a", "b", "c"), size = 15, x = c(0,0,0), y=c(1,0.7,0.4))```
```dev.off()```
