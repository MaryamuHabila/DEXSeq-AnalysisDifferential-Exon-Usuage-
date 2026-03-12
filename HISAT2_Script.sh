#!/bin/bash
set -euo pipefail

# ----------------------------
# Paths
# ----------------------------
FASTA="/home/s2451842/GTF/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF="/home/s2451842/GTF/Homo_sapiens.GRCh38.112.gtf"
INDEX_PREFIX="Homo_sapiens.GRCh38_direct_index"
THREADS=8
SAMPLESHEET="Samplesheet_HAoSMCs.csv"
ADAPTERS="adapters/TruSeq3-PE.fa"

WORKDIR="$(pwd)"
OUTPUT_DIR="${WORKDIR}/HAoSMC_HISAT2_Trimmomatic_Output_Ensembl112"
TRIMDIR="${WORKDIR}/trimmed_data"
QC_DIR="${WORKDIR}/qc_reports"
LOG_DIR="${WORKDIR}/logs"
MQC_DIR="${WORKDIR}/multiqc_reports"

mkdir -p "$OUTPUT_DIR" "$TRIMDIR" "$QC_DIR" "$LOG_DIR" "$MQC_DIR" adapters

# ----------------------------
# Download adapters if missing
# ----------------------------
if [ ! -f "$ADAPTERS" ]; then
  echo "Downloading adapters..."
  wget -P adapters/ https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa
fi

# ----------------------------
# FastQC
# ----------------------------
echo "Running FastQC..."
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r sample fastq_1 fastq_2 strandedness
do
  fastqc -t "$THREADS" -o "$QC_DIR" "$fastq_1" "$fastq_2"
done

# ----------------------------
# Trimmomatic (SKIP if trimmed exists)
# ----------------------------
echo "Checking trimming status..."
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r sample fastq_1 fastq_2 strandedness
do
  out1="${TRIMDIR}/${sample}_1P.fastq.gz"
  out2="${TRIMDIR}/${sample}_2P.fastq.gz"

  if [[ -s "$out1" && -s "$out2" ]]; then
    echo "✔ Trimmed reads exist for $sample — skipping trimming"
    continue
  fi

  echo "✂ Trimming $sample"
  trimmomatic PE -threads "$THREADS" \
    "$fastq_1" "$fastq_2" \
    "${TRIMDIR}/${sample}_1P.fastq.gz" "${TRIMDIR}/${sample}_1U.fastq.gz" \
    "${TRIMDIR}/${sample}_2P.fastq.gz" "${TRIMDIR}/${sample}_2U.fastq.gz" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 \
    > "${LOG_DIR}/${sample}.trimmomatic.log" 2>&1

  test -s "$out1"
  test -s "$out2"
done

# ----------------------------
# Build HISAT2 index (skip if exists)
# ----------------------------
if [ ! -f "${INDEX_PREFIX}.1.ht2" ]; then
  echo "Building HISAT2 index..."
  hisat2_extract_splice_sites.py "$GTF" > splicesites.txt
  hisat2_extract_exons.py "$GTF" > exons.txt
  hisat2-build --ss splicesites.txt --exon exons.txt \
    "$FASTA" "$INDEX_PREFIX" 2>&1 | tee "${LOG_DIR}/hisat2_build.log"
fi

# ----------------------------
# HISAT2 Alignment (skip if BAM exists)
# ----------------------------
echo "Starting alignment..."
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r sample fastq_1 fastq_2 strandedness
do
  bam="${OUTPUT_DIR}/${sample}.sorted.bam"

  if [[ -s "$bam" ]]; then
    echo "✔ BAM exists for $sample — skipping alignment"
    continue
  fi

  echo "🔹 Aligning $sample"
  hisat2 -x "$INDEX_PREFIX" \
    -1 "${TRIMDIR}/${sample}_1P.fastq.gz" \
    -2 "${TRIMDIR}/${sample}_2P.fastq.gz" \
    --rna-strandness RF \
    --summary-file "${OUTPUT_DIR}/${sample}.hisat2.summary.log" \
    --threads "$THREADS" \
    --rg-id "$sample" --rg SM:"$sample" \
    --no-mixed --no-discordant --sensitive -I 1 -X 1000 \
    | samtools view -@ "$THREADS" -bS - > "${OUTPUT_DIR}/${sample}.bam"

  test -s "${OUTPUT_DIR}/${sample}.bam"

  samtools sort -@ "$THREADS" \
    -o "$bam" "${OUTPUT_DIR}/${sample}.bam"
  samtools index "$bam"
done

# ----------------------------
# MultiQC
# ----------------------------
multiqc -o "$MQC_DIR" "$QC_DIR" "$OUTPUT_DIR" "$LOG_DIR"

echo "Pipeline completed successfully."
