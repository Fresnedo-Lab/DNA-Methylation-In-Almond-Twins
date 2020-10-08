## RNASeq Data Alignment

**Software:**

FASTQC v. 0.11.7\
STAR v. 2.6.0a\
Picard v. 2.3.0\
Subread v. 1.5.0-p2\


Perform fastqc analysis on all RNAseq data files for each Stukey twin before alignment. Prepare the Nonpareil almond genome prior to use in alignment.

```/usr/local/star/2.6.0a/bin/Linux_x86_64/STAR --runThreadN 28 --runMode genomeGenerate --genomeDir /path/to/genome/directory --genomeFastaFiles /path/to/genome/fasta/file --sjdbGTFfile /path/to/gff/file --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 --genomeSAindexNbases 12```

Perform alignment using almond reference genome and each technical replicate RNASeq data file for each Stukey twin.

```/usr/local/star/2.6.0a/bin/Linux_x86_64/STAR --runThreadN 28 --genomeDir /path/to/genome/directory --readFilesIn /path/to/forward/read/fastq /path/to/reverse/read/fastq --readFilesCommand zcat --outFileNamePrefix /path/to/output/file --alignIntronMin 60 --alignIntronMax 6000 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outSAMtype BAM Unsorted SortedByCoordinate```

Mark duplicate reads using picard and the bam files generated above. This generates a new bam file to use for additional analyses.

```java -jar $PICARD MarkDuplicates I=/path/to/soted/bam O=/path/to/output/file M=/path/to/output/metrics/file```

Generate a feature counts file for each technical replicate using the gff file for the Nonpareil almond reference genome. For this particular analysis, counts for gene features are generated.

```featureCounts -T 12 -a /path/to/reference/gff/file -t gene -f -p -M -g ID -o /path/to/output/file /path/to/sorted/bam/file```
