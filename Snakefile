import os
import csv

configfile: "config.yaml"

# -----------------------------------
# Read sample sheet
# -----------------------------------
SAMPLES = {}
with open(config["samplesheet"]) as f:
    reader = csv.DictReader(f)
    for row in reader:
        SAMPLES[row["sample"]] = {
            "fastq_1": row["fastq_1"],
            "fastq_2": row["fastq_2"],
            "strandedness": row.get("strandedness", "reverse")
        }

SAMPLE_NAMES = list(SAMPLES.keys())

# -----------------------------------
# Variables
# -----------------------------------
FASTA = config["fasta"]
GTF = config["gtf"]
THREADS = config["threads"]
ADAPTERS = config["adapters"]
INDEX_PREFIX = config["index_prefix"]

TRIMDIR = config["trimdir"]
QC_DIR = config["qc_dir"]
ALIGN_DIR = config["align_dir"]
LOG_DIR = config["log_dir"]
MULTIQC_DIR = config["multiqc_dir"]

DEXSEQ_PROJECT_DIR = config["dexseq_project_dir"]
DEXSEQ_ANNOT_DIR = config["dexseq_annot_dir"]
DEXSEQ_COUNT_DIR = config["dexseq_count_dir"]
DEXSEQ_LOG_DIR = config["dexseq_log_dir"]

FLATTENED_GFF = f"{DEXSEQ_ANNOT_DIR}/Homo_sapiens.GRCh38.112.DEXSeq.gff"

# -----------------------------------
# Final target
# -----------------------------------
rule all:
    input:
        expand(f"{ALIGN_DIR}/{{sample}}.sorted.bam", sample=SAMPLE_NAMES),
        expand(f"{ALIGN_DIR}/{{sample}}.sorted.bam.bai", sample=SAMPLE_NAMES),
        expand(f"{DEXSEQ_COUNT_DIR}/{{sample}}.dexseq_counts.txt", sample=SAMPLE_NAMES),
        f"{MULTIQC_DIR}/multiqc_report.html",
        FLATTENED_GFF

# -----------------------------------
# Download adapters
# -----------------------------------
rule download_adapters:
    output:
        ADAPTERS
    shell:
        """
        mkdir -p adapters
        wget -O {output} https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa
        """

# -----------------------------------
# FastQC
# -----------------------------------
rule fastqc:
    input:
        r1=lambda wc: SAMPLES[wc.sample]["fastq_1"],
        r2=lambda wc: SAMPLES[wc.sample]["fastq_2"]
    output:
        html1=f"{QC_DIR}/{{sample}}_R1_fastqc.html",
        zip1=f"{QC_DIR}/{{sample}}_R1_fastqc.zip",
        html2=f"{QC_DIR}/{{sample}}_R2_fastqc.html",
        zip2=f"{QC_DIR}/{{sample}}_R2_fastqc.zip"
    threads: THREADS
    params:
        outdir=QC_DIR
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2}
        """

# -----------------------------------
# Trimmomatic
# -----------------------------------
rule trimmomatic:
    input:
        r1=lambda wc: SAMPLES[wc.sample]["fastq_1"],
        r2=lambda wc: SAMPLES[wc.sample]["fastq_2"],
        adapters=ADAPTERS
    output:
        p1=f"{TRIMDIR}/{{sample}}_1P.fastq.gz",
        u1=f"{TRIMDIR}/{{sample}}_1U.fastq.gz",
        p2=f"{TRIMDIR}/{{sample}}_2P.fastq.gz",
        u2=f"{TRIMDIR}/{{sample}}_2U.fastq.gz"
    log:
        f"{LOG_DIR}/{{sample}}.trimmomatic.log"
    threads: THREADS
    shell:
        """
        mkdir -p {TRIMDIR} {LOG_DIR}
        trimmomatic PE -threads {threads} \
          {input.r1} {input.r2} \
          {output.p1} {output.u1} \
          {output.p2} {output.u2} \
          ILLUMINACLIP:{input.adapters}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 \
          > {log} 2>&1
        """

# -----------------------------------
# Extract splice sites
# -----------------------------------
rule splice_sites:
    input:
        GTF
    output:
        "hisat2_index/splicesites.txt"
    shell:
        """
        mkdir -p hisat2_index
        hisat2_extract_splice_sites.py {input} > {output}
        """

# -----------------------------------
# Extract exons
# -----------------------------------
rule exons:
    input:
        GTF
    output:
        "hisat2_index/exons.txt"
    shell:
        """
        mkdir -p hisat2_index
        hisat2_extract_exons.py {input} > {output}
        """

# -----------------------------------
# Build HISAT2 index
# -----------------------------------
rule build_hisat2_index:
    input:
        fasta=FASTA,
        ss="hisat2_index/splicesites.txt",
        exons="hisat2_index/exons.txt"
    output:
        expand("hisat2_index/Homo_sapiens.GRCh38_direct_index.{i}.ht2", i=range(1, 9))
    log:
        f"{LOG_DIR}/hisat2_build.log"
    threads: THREADS
    shell:
        """
        mkdir -p hisat2_index {LOG_DIR}
        hisat2-build --ss {input.ss} --exon {input.exons} \
          {input.fasta} {INDEX_PREFIX} > {log} 2>&1
        """

# -----------------------------------
# Align
# -----------------------------------
rule hisat2_align:
    input:
        idx=expand("hisat2_index/Homo_sapiens.GRCh38_direct_index.{i}.ht2", i=range(1, 9)),
        r1=f"{TRIMDIR}/{{sample}}_1P.fastq.gz",
        r2=f"{TRIMDIR}/{{sample}}_2P.fastq.gz"
    output:
        bam=temp(f"{ALIGN_DIR}/{{sample}}.bam"),
        summary=f"{ALIGN_DIR}/{{sample}}.hisat2.summary.log"
    threads: THREADS
    params:
        strandness=config["hisat2_rna_strandness"]
    shell:
        """
        mkdir -p {ALIGN_DIR}
        hisat2 -x {INDEX_PREFIX} \
          -1 {input.r1} \
          -2 {input.r2} \
          --rna-strandness {params.strandness} \
          --summary-file {output.summary} \
          --threads {threads} \
          --rg-id {wildcards.sample} --rg SM:{wildcards.sample} \
          --no-mixed --no-discordant --sensitive -I 1 -X 1000 \
          | samtools view -@ {threads} -bS - > {output.bam}
        """

# -----------------------------------
# Sort BAM
# -----------------------------------
rule sort_bam:
    input:
        f"{ALIGN_DIR}/{{sample}}.bam"
    output:
        f"{ALIGN_DIR}/{{sample}}.sorted.bam"
    threads: THREADS
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """

# -----------------------------------
# Index BAM
# -----------------------------------
rule index_bam:
    input:
        f"{ALIGN_DIR}/{{sample}}.sorted.bam"
    output:
        f"{ALIGN_DIR}/{{sample}}.sorted.bam.bai"
    shell:
        """
        samtools index {input}
        """

# -----------------------------------
# MultiQC
# -----------------------------------
rule multiqc:
    input:
        expand(f"{ALIGN_DIR}/{{sample}}.hisat2.summary.log", sample=SAMPLE_NAMES),
        expand(f"{LOG_DIR}/{{sample}}.trimmomatic.log", sample=SAMPLE_NAMES)
    output:
        f"{MULTIQC_DIR}/multiqc_report.html"
    shell:
        """
        mkdir -p {MULTIQC_DIR}
        multiqc -o {MULTIQC_DIR} {QC_DIR} {ALIGN_DIR} {LOG_DIR}
        """

# -----------------------------------
# Find DEXSeq scripts
# -----------------------------------
rule dexseq_prepare_annotation:
    input:
        gtf=GTF
    output:
        FLATTENED_GFF
    log:
        stdout=f"{DEXSEQ_LOG_DIR}/dexseq_prepare_annotation.stdout.log",
        stderr=f"{DEXSEQ_LOG_DIR}/dexseq_prepare_annotation.stderr.log"
    shell:
        r"""
        mkdir -p {DEXSEQ_ANNOT_DIR} {DEXSEQ_COUNT_DIR} {DEXSEQ_LOG_DIR}

        DEXSEQ_PY=$(Rscript -e 'cat(system.file("python_scripts", package="DEXSeq"))')
        PREP_SCRIPT="${{DEXSEQ_PY}}/dexseq_prepare_annotation.py"

        python "$PREP_SCRIPT" \
          {input.gtf} \
          {output} \
          > {log.stdout} 2> {log.stderr}
        """

# -----------------------------------
# DEXSeq counting
# -----------------------------------
rule dexseq_count:
    input:
        bam=f"{ALIGN_DIR}/{{sample}}.sorted.bam",
        bai=f"{ALIGN_DIR}/{{sample}}.sorted.bam.bai",
        gff=FLATTENED_GFF
    output:
        f"{DEXSEQ_COUNT_DIR}/{{sample}}.dexseq_counts.txt"
    log:
        stdout=f"{DEXSEQ_LOG_DIR}/{{sample}}.dexseq_count.stdout.log",
        stderr=f"{DEXSEQ_LOG_DIR}/{{sample}}.dexseq_count.stderr.log"
    params:
        stranded=config["dexseq_stranded"]
    shell:
        r"""
        samtools quickcheck {input.bam}

        DEXSEQ_PY=$(Rscript -e 'cat(system.file("python_scripts", package="DEXSeq"))')
        COUNT_SCRIPT="${{DEXSEQ_PY}}/dexseq_count.py"

        python "$COUNT_SCRIPT" \
          -p yes \
          -r pos \
          -s {params.stranded} \
          -f bam \
          {input.gff} \
          {input.bam} \
          {output} \
          > {log.stdout} 2> {log.stderr}
        """
