## Generate Chromosome Level Methylation Profiles

Create plots showing the percent total methylation in all contexts across each of the 8 almond chromosomes for each Stukey twin using R and the CX reports genreated for each methylation context and twin in Bismark.

**Software:**
R v. 3.6.3

Read in the Bismark CX report for the CG context for each Stukey twin.

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

Create a new variable in the aggregate data frames representing percent methylation (count methylated/total count)

```agg_1a$pm<- agg_1a$meth / agg_1a$total```
```agg_1b$pm<- agg_1b$meth / agg_1b$total```
```agg_2a$pm<- agg_2a$meth / agg_2a$total```
```agg_2b$pm<- agg_2b$meth / agg_2b$total```

Write the aggreated data frame to a text file for each twin pair.

```write.txt(agg_1a, "CG_1a_agg.txt", sep="\t")```
```write.txt(agg_1b, "CG_1b_agg.txt", sep="\t")```
```write.txt(agg_2a, "CG_2a_agg.txt", sep="\t")```
```write.txt(agg_2b, "CG_2b_agg.txt", sep="\t")```

**Repeat the above for the CHG and CHH methylation contexts**
