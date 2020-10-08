## Generate Chromosome Level Methylation Profiles

Create plots showing the percent total methylation in all contexts across each of the 8 almond chromosomes for each Stukey twin using R and the CX reports genreated for each methylation context and twin in Bismark.

**Software:**
R v. 3.6.3

Read in the CX report for the CG context for each Stukey twin.

```CX_1a<- read.delim("CHGStukey1a_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)```

```CX_1b<- read.delim("CHGStukey1b_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)```

```CX_2a<- read.delim("CHGStukey2a_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)```

```CX_2b<- read.delim("CHGStukey2b_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)```

Change the column names in each dataframe.

```colnames(CX_1a)<- c("Scaffold","pos","meth","unmeth")```
```colnames(CX_2a)<- c("Scaffold","pos","meth","unmeth")```
```colnames(CX_1b)<- c("Scaffold","pos","meth","unmeth")```
```colnames(CX_2b)<- c("Scaffold","pos","meth","unmeth")```

Calculate the total number of counts for each cytosine (methylated and unmethylated counts).

```CX_1a$total<- CX_1a$meth + CX_1a$unmeth```
```CX_1b$total<- CX_1b$meth + CX_1b$unmeth```
```CX_2a$total<- CX_2a$meth + CX_2a$unmeth```
```CX_2b$total<- CX_2b$meth + CX_2b$unmeth```


agg_1a<- aggregate(. ~ Scaffold + pos, data=CX_1a, sum)
agg_1b<- aggregate(. ~ Scaffold + pos, data=CX_1b, sum)
agg_2a<- aggregate(. ~ Scaffold + pos, data=CX_2a, sum)
agg_1b<- aggregate(. ~ Scaffold + pos, data=CX_2b, sum)
