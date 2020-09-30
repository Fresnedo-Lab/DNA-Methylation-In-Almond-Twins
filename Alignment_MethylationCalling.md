## Alignment & Methylation Calling

Bash script used to align concatenated Stukey reads to the almond genome and call methylated cytosines in each methylation context.

**Software:**   
FastQC v. 0.11.7\
TrimGalore v. 0.6.0\
Cutadapt v. 2.1\
Bismark v. 0.21.0

Perform FastQC analysis before and after trimming the fastq files. Trimmed single end files for each sample using TrimGalore with additional 5' and 3' clip for each read.

```trim_galore --three_prime_clip_R1 10 ----clip_R1 10 R1.fastq.gz```

Align paired end data following trimming to almond 'Nonpareil' genome using Bismark. Before aligning, perform genome preparation for reference genome of interest.

```bismark_genome_preparation --path_to_aligner /usr/bin/bowtie2/ --verbose /path/to/reference/genome```

```bismark -n 1 /path/to/Bismark/genome/files --non_directional Stukey1.fastq.gz -o Nonpareil_Bismark_Alignment```

Once all files are aligned to reference genome, perform deduplication of each bam file.

```deduplicate_bismark -s Nonpareil_Bismark_Alignment Stukey1_bismark_bt2.bam -o Nonpareil_Bismark_Alignment/DeduplicatedFiles```

Extract methylation calls producing reports for each methylation context (CG, CHG, & CHH).

```bismark_methylation_extractor Stukey2_bismark_bt2.deduplicated.bam -s --comprehensive --samtools_path /usr/local/samtools/1.9/bin --output ExtractedMethylationFiles_BedReport --parallel 9 --bedGraph --CX_context --scaffolds --cytosine_report --genome_folder /path/to/Bismark/genome/files```
