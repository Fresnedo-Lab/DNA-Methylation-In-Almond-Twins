## Create Heatmaps and Venn Diagrams of DMRs for Stukey Twinpairs

Read in the csv file generated when calling DMRs using DSS for each methylation context and twin pair.

```CG_TP1 <- read.csv("../../Smoothing CG_DMRs_Pair1_0.01_DSS.csv", header = TRUE)```
```CG_TP2 <- read.csv("../../Smoothing CG_DMRs_Pair2_0.01_DSS.csv", header = TRUE)```

Read in the csv file containing the DMRs associated with the same gene in the same proximitiy class between the two twinpairs in each methylation context.

```CG <- read.csv("ObsCG_All.csv", header = TRUE)```

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
