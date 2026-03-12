# Differential Exon Usage Analysis Pipeline on RNA-seq data (HISAT2 to exon-level counting to statistical analyses) 
DEXSeq tests whether the proportion of reads mapping to a specific exon relative to the rest of the gene changes between conditions.


## Overview

This project performs **RNA-seq analysis to detect differential exon usage (DEU)** between **control** and **knockdown (KD)** human aortic smooth muscle cell (HAoSMC) samples.

The analysis is performed in three major stages:

1. **Read preprocessing and alignment**
2. **Generation of exon-level count data for DEXSeq**
3. **Statistical testing for differential exon usage in R**

The workflow combines several bioinformatics tools and programming environments:

| Tool                   | Purpose                                          |
| ---------------------- | ------------------------------------------------ |
| FastQC                 | Assess quality of raw sequencing reads           |
| Trimmomatic            | Remove sequencing adapters and low-quality bases |
| HISAT2                 | Align RNA-seq reads to the human genome          |
| SAMtools               | Process BAM alignment files                      |
| MultiQC                | Summarize quality control metrics                |
| HTSeq / DEXSeq scripts | Generate exon-level counts                       |
| R (DEXSeq package)     | Perform statistical analysis of exon usage       |

The goal of this analysis is to identify **changes in exon usage between control and knockdown conditions**, which may indicate **alternative splicing or transcript isoform changes**.

---

# Biological Objective

RNA-seq experiments measure the abundance of RNA molecules in biological samples.
However, many genes produce **multiple transcript isoforms** due to **alternative splicing**.

Alternative splicing can cause:

* exon inclusion
* exon skipping
* alternative splice sites
* transcript isoform switching

DEXSeq detects whether **individual exons within a gene are used differently between conditions**.

This allows us to identify **regulatory effects on RNA splicing** caused by the knockdown experiment.

---

# Overall Analysis Workflow

The analysis proceeds in the following steps:

```
Raw FASTQ reads
        ↓
FastQC quality assessment
        ↓
Adapter trimming (Trimmomatic)
        ↓
Genome alignment (HISAT2)
        ↓
Sorted BAM files (SAMtools)
        ↓
DEXSeq exon counting
        ↓
DEXSeq statistical analysis in R
        ↓
Differential exon usage results
```

---

# Directory Structure

The project uses the following directory structure.

```
Scripts/
│
├── Final_working_hisat2_Ensembl112_haosmc.sh
├── HAoSMCs_KD_DEXSeq_Analysis.sh
│
├── Samplesheet_HAoSMCs.csv
│
├── HAoSMC_HISAT2_Trimmomatic_Output_Ensembl112/
├── trimmed_data/
├── qc_reports/
├── logs/
├── multiqc_reports/
│
└── DEXSeq_HAoSMCs_KD_112/
      ├── annotation/
      ├── counts/
      └── logs/
```

---

# Step 1 — RNA-seq Alignment Pipeline

Script used:

```
Final_working_hisat2_Ensembl112_haosmc.sh
```

This script processes raw FASTQ files and produces **sorted BAM alignment files**.

## Input Files

The script requires:

### Reference genome

```
/home/s2451842/GTF/Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

### Gene annotation

```
/home/s2451842/GTF/Homo_sapiens.GRCh38.112.gtf
```

### Sample sheet

```
Samplesheet_HAoSMCs.csv
```

Example format:

```
sample,fastq_1,fastq_2,strandedness
CTR1_Z20_T,file_R1.fastq.gz,file_R2.fastq.gz,reverse
KD1_Z20_T,file_R1.fastq.gz,file_R2.fastq.gz,reverse
```

---

## Step 1.1 — FastQC Quality Control

FastQC is used to evaluate the quality of the raw sequencing reads.

The script runs:

```
fastqc -t 8 -o qc_reports fastq_1 fastq_2
```

Outputs include:

```
qc_reports/*.html
qc_reports/*.zip
```

These reports contain information about:

* sequence quality scores
* adapter contamination
* GC content
* duplication levels

---

## Step 1.2 — Adapter Trimming

Sequencing adapters and low-quality bases are removed using **Trimmomatic**.

Command used:

```
trimmomatic PE
```

Important parameters:

```
ILLUMINACLIP:TruSeq3-PE.fa
SLIDINGWINDOW:4:20
MINLEN:50
```

Outputs:

```
trimmed_data/sample_1P.fastq.gz
trimmed_data/sample_2P.fastq.gz
```

Only **paired trimmed reads (P)** are used for alignment.

---

## Step 1.3 — HISAT2 Genome Index Construction

Before alignment, a genome index is built.

HISAT2 uses the reference genome and gene annotation to construct a splice-aware index.

Commands used:

```
hisat2_extract_splice_sites.py
hisat2_extract_exons.py
hisat2-build
```

This improves alignment across splice junctions.

---

## Step 1.4 — RNA-seq Alignment

Trimmed reads are aligned to the reference genome using **HISAT2**.

Important parameters:

```
--rna-strandness RF
--no-mixed
--no-discordant
--sensitive
```

The aligner outputs BAM files via SAMtools.

---

## Step 1.5 — Sorting and Indexing

SAMtools is used to process alignment files.

Steps performed:

```
samtools view
samtools sort
samtools index
```

Final output files:

```
sample.sorted.bam
sample.sorted.bam.bai
```

These BAM files contain mapped RNA-seq reads.

---

## Step 1.6 — MultiQC Summary

MultiQC aggregates all QC reports.

Command:

```
multiqc qc_reports alignment logs
```

Output:

```
multiqc_reports/multiqc_report.html
```

This provides a summary of:

* FastQC metrics
* alignment statistics
* trimming results

---

# Step 2 — DEXSeq Exon Counting

Script used:

```
HAoSMCs_KD_DEXSeq_Analysis.sh
```

This script prepares data for differential exon usage analysis.

---

## Step 2.1 — Flattening the Annotation

DEXSeq requires a **flattened annotation file**.

This converts the original GTF into **non-overlapping exon bins**.

Command used:

```
dexseq_prepare_annotation.py
```

Output:

```
Homo_sapiens.GRCh38.112.DEXSeq.gff
```

Each exon is assigned a unique identifier.

---

## Step 2.2 — Generating Exon Counts

DEXSeq counts reads overlapping each exon bin.

Command used:

```
dexseq_count.py
```

Important parameters:

```
-p yes
-r pos
-s reverse
-f bam
```

Outputs:

```
sample.dexseq_counts.txt
```

Each file contains counts per exon bin.

Example entry:

```
ENSG000001234:001    42
```

---

# Step 3 — Differential Exon Usage Analysis in R

The final analysis is performed in **R using the DEXSeq package**.

---

## Step 3.1 — Load Required Packages

Packages used:

```
DEXSeq
BiocParallel
GenomicRanges
IRanges
rtracklayer
```

These packages perform statistical modeling and genomic data handling.

---

## Step 3.2 — Sample Table

Samples are assigned to experimental groups.

```
control
KD
```

Example:

```
CTR1_Z20_T   control
CTR2_Z30_T   control
KD1_Z20_T    KD
KD2_Z30_T    KD
```

---

## Step 3.3 — Load Exon Count Files

DEXSeq count files are read and converted into a matrix.

Rows represent:

```
exon bins
```

Columns represent:

```
samples
```

Example matrix:

```
exon_bin   CTR1   CTR2   KD1   KD2
ENSG:001    34     40     10    12
ENSG:002    50     47     49    51
```

---

## Step 3.4 — Create DEXSeq Dataset

The dataset object is constructed:

```
DEXSeqDataSet()
```

Design formula:

```
~ sample + exon + condition:exon
```

This model tests whether exon usage differs between conditions.

---

## Step 3.5 — Statistical Analysis

DEXSeq performs several statistical steps.

1. **Estimate size factors**

Normalizes library sizes.

2. **Estimate dispersions**

Models variability in exon counts.

3. **Test for differential exon usage**

Identifies exons whose usage changes between conditions.

4. **Estimate exon fold changes**

Calculates effect size.

---

# Step 4 — Output Results

Several outputs are generated.

---

## Full DEXSeq Dataset

```
HAoSMC_DEXSeq_dxd.rds
```

Contains the full DEXSeq object.

---

## Results Table

```
HAoSMC_DEXSeq_results_annotated.csv
```

Columns include:

| Column         | Description           |
| -------------- | --------------------- |
| geneID         | Gene identifier       |
| exonID         | Exon bin              |
| chr            | Chromosome            |
| start          | Start coordinate      |
| end            | End coordinate        |
| log2FoldChange | Exon usage change     |
| pvalue         | Raw p-value           |
| padj           | FDR-corrected p-value |

---

## Significant Results

```
HAoSMC_DEXSeq_significant_FDR0.05.csv
```

Contains exons with:

```
padj < 0.05
```

These represent **statistically significant differential exon usage events**.

---

# Step 5 — Visualization

DEXSeq generates diagnostic plots.

---

## Dispersion Plot

```
HAoSMC_DEXSeq_dispersion_plot.pdf
```

Shows dispersion estimates used in statistical modeling.

---

## MA Plot

```
HAoSMC_DEXSeq_MAplot.pdf
```

Displays log fold change versus mean expression.

---

## Top Gene Plot

```
HAoSMC_DEXSeq_top_gene_plot.pdf
```

Shows exon usage differences for the most significant gene.

---

# Step 6 — Region of Interest (ROI) Analysis

The script also extracts results for a **specific genomic region**.

Example:

```
chr12
start = ROI_start
end = ROI_end
```

Exons overlapping this region are identified.

Outputs include:

```
ROI_chr12_wholelocus_results.csv
ROI_wholelocus_significant_FDR0.05.csv
ROI_chr12_wholelocus_counts.csv
```

These files allow focused investigation of a locus of interest.

---

# Step 7 — Gene-level Visualization for ROI

For each gene overlapping the ROI:

DEXSeq plots are generated:

```
plotDEXSeq()
```

Outputs:

```
ROI_chr12_wholelocus_gene_plots/*.pdf
```

These plots show exon usage differences across conditions.

---

# Final Outcome

This workflow produces:

* genome-aligned RNA-seq reads
* exon-level count tables
* statistical tests for differential exon usage
* annotated genomic results
* visualization of exon usage patterns

The results allow identification of **genes whose splicing patterns change in response to knockdown**, providing insight into **RNA processing and regulatory mechanisms**.

---

# Author

Maryam Uhabila Usman
RNA-seq Differential Exon Usage Analysis
HAoSMC Knockdown Study
