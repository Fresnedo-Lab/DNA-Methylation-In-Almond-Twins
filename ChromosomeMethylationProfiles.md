## Generate Chromosome Level Methylation Profiles

Create plots showing the percent total methylation in all contexts across each of the 8 almond chromosomes for each Stukey twin using R and the CX reports genreated for each methylation context and twin in Bismark.

**Software:**
R v. 3.6.3

# This script aggregates tech reps for use in later analyses

Define column classes
``cc<- rep('NULL', 7)``
``cc[1]<- 'factor'``
``cc[c(2,4,5)]<- 'integer'``

Read in CX reports generated from Bismark methylation extractor
``CX_1a<- read.delim("CGStukey1a_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)``
``CX_1b<- read.delim("CGStukey1b_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)``
``CX_2a<- read.delim("CGStukey2a_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)``
``CX_2b<- read.delim("CGStukey2b_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)``

Rename columns in data frames
``colnames(CX_1a)<- c("Scaffold","pos","meth","unmeth")``
``colnames(CX_2a)<- c("Scaffold","pos","meth","unmeth")``
``colnames(CX_1b)<- c("Scaffold","pos","meth","unmeth")``
``colnames(CX_2b)<- c("Scaffold","pos","meth","unmeth")``

Aggregate tech reps. This may take awhile.
``agg_1a<- aggregate(. ~ Scaffold + pos, data=CX_1a, sum)``
``agg_1b<- aggregate(. ~ Scaffold + pos, data=CX_1b, sum)``
``agg_2a<- aggregate(. ~ Scaffold + pos, data=CX_2a, sum)``
``agg_2b<- aggregate(. ~ Scaffold + pos, data=CX_2b, sum)``

Write to table.
``write.table(agg_1a, "agg_1a.txt", sep='\t', row.names=F, quote=F)``
``write.table(agg_1b, "agg_1b.txt", sep='\t', row.names=F, quote=F)``
``write.table(agg_2a, "agg_2a.txt", sep='\t', row.names=F, quote=F)``
``write.table(agg_2b, "agg_2b.txt", sep='\t', row.names=F, quote=F)``


Repeat above for CHG context

``CX_1a<- read.delim("CHGStukey1a_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)``
``CX_1b<- read.delim("CHGStukey1b_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)``
``CX_2a<- read.delim("CHGStukey2a_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)``
``CX_2b<- read.delim("CHGStukey2b_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)``

``colnames(CX_1a)<- c("Scaffold","pos","meth","unmeth")``
``colnames(CX_2a)<- c("Scaffold","pos","meth","unmeth")``
``colnames(CX_1b)<- c("Scaffold","pos","meth","unmeth")``
``colnames(CX_2b)<- c("Scaffold","pos","meth","unmeth")``

``agg_1a<- aggregate(. ~ Scaffold + pos, data=CX_1a, sum)``
``agg_1b<- aggregate(. ~ Scaffold + pos, data=CX_1b, sum)``
``agg_2a<- aggregate(. ~ Scaffold + pos, data=CX_2a, sum)``
``agg_2b<- aggregate(. ~ Scaffold + pos, data=CX_2b, sum)``

``write.table(agg_1a, "agg_1a.txt", sep='\t', row.names=F, quote=F)``
``write.table(agg_1b, "agg_1b.txt", sep='\t', row.names=F, quote=F)``
``write.table(agg_2a, "agg_2a.txt", sep='\t', row.names=F, quote=F)``
``write.table(agg_2b, "agg_2b.txt", sep='\t', row.names=F, quote=F)``


CHH context

``CX_1a<- read.delim("CHHStukey1a_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)``
``CX_1b<- read.delim("CHHStukey1b_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)``
``CX_2a<- read.delim("CHHStukey2a_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)``
``CX_2b<- read.delim("CHHStukey2b_CXreport.txt", header=F, sep="\t", dec=".", colClasses=cc)``

``colnames(CX_1a)<- c("Scaffold","pos","meth","unmeth")``
``colnames(CX_2a)<- c("Scaffold","pos","meth","unmeth")``
``colnames(CX_1b)<- c("Scaffold","pos","meth","unmeth")``
``colnames(CX_2b)<- c("Scaffold","pos","meth","unmeth")``

``agg_1a<- aggregate(. ~ Scaffold + pos, data=CX_1a, sum)``
``agg_1b<- aggregate(. ~ Scaffold + pos, data=CX_1b, sum)``
``agg_2a<- aggregate(. ~ Scaffold + pos, data=CX_2a, sum)``
``agg_2b<- aggregate(. ~ Scaffold + pos, data=CX_2b, sum)``

``write.table(agg_1a, "agg_1a.txt", sep='\t', row.names=F, quote=F)``
``write.table(agg_1b, "agg_1b.txt", sep='\t', row.names=F, quote=F)``
``write.table(agg_2a, "agg_2a.txt", sep='\t', row.names=F, quote=F)``
``write.table(agg_2b, "agg_2b.txt", sep='\t', row.names=F, quote=F)``
