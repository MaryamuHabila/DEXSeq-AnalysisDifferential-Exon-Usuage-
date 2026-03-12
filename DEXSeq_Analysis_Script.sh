#!/bin/bash
set -euo pipefail

# ======================================
# Activate conda environment
# ======================================

eval "$(conda shell.bash hook)"
conda activate base

# ======================================
# DEXSeq exon count generation
# ======================================

# Paths
GTF="/home/s2451842/GTF/Homo_sapiens.GRCh38.112.gtf"
BAM_DIR="/home/s2451842/projects/HAoSMC_KD_150PE_GENCODE/HISAT2_Hardcoded_Analyses2/HISAT2_112_Analyses/HAoSMC_HISAT2_Trimmomatic_Output_Ensembl112"
PROJECT_DIR="/home/s2451842/projects/Monzino_Analyses/HISAT2_HTSeq_Analyses/DEXSeq_HAoSMCs_KD_112"

ANNOT_DIR="${PROJECT_DIR}/annotation"
COUNT_DIR="${PROJECT_DIR}/counts"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$ANNOT_DIR" "$COUNT_DIR" "$LOG_DIR"

FLATTENED_GFF="${ANNOT_DIR}/Homo_sapiens.GRCh38.112.DEXSeq.gff"

echo "======================================"
echo "DEXSeq exon counting pipeline"
echo "Project directory: $PROJECT_DIR"
echo "GTF: $GTF"
echo "BAM directory: $BAM_DIR"
echo "Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"
echo "Python used: $(which python)"
echo "======================================"

# ======================================
# Check required software
# ======================================

command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not installed"; exit 1; }
command -v python >/dev/null 2>&1 || { echo "ERROR: python not installed"; exit 1; }
command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not installed"; exit 1; }

# Check HTSeq python module
python -c "import HTSeq" >/dev/null 2>&1 || {
    echo "ERROR: HTSeq Python module not found in this environment"
    echo "Install with: conda install -c bioconda htseq"
    exit 1
}

echo "HTSeq detected successfully"

# ======================================
# Locate DEXSeq python scripts
# ======================================

DEXSEQ_PY=$(Rscript -e 'cat(system.file("python_scripts", package="DEXSeq"))')

if [[ -z "$DEXSEQ_PY" || ! -d "$DEXSEQ_PY" ]]; then
    echo "ERROR: Could not find DEXSeq python_scripts directory"
    echo "Install DEXSeq in R with:"
    echo "BiocManager::install(\"DEXSeq\")"
    exit 1
fi

PREP_SCRIPT="${DEXSEQ_PY}/dexseq_prepare_annotation.py"
COUNT_SCRIPT="${DEXSEQ_PY}/dexseq_count.py"

[[ -f "$PREP_SCRIPT" ]] || { echo "ERROR: Missing $PREP_SCRIPT"; exit 1; }
[[ -f "$COUNT_SCRIPT" ]] || { echo "ERROR: Missing $COUNT_SCRIPT"; exit 1; }

echo "DEXSeq scripts found in:"
echo "$DEXSEQ_PY"

# ======================================
# Create flattened annotation
# ======================================

if [[ ! -s "$FLATTENED_GFF" ]]; then

    echo "--------------------------------------"
    echo "Creating flattened annotation"
    echo "--------------------------------------"

    python "$PREP_SCRIPT" \
        "$GTF" \
        "$FLATTENED_GFF" \
        > "${LOG_DIR}/dexseq_prepare_annotation.stdout.log" \
        2> "${LOG_DIR}/dexseq_prepare_annotation.stderr.log"

    echo "Flattened annotation created:"
    echo "$FLATTENED_GFF"

else

    echo "Flattened annotation already exists"
    echo "$FLATTENED_GFF"

fi

# ======================================
# Generate exon counts
# ======================================

shopt -s nullglob
BAM_FILES=( "$BAM_DIR"/*.sorted.bam )

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No .sorted.bam files found in:"
    echo "$BAM_DIR"
    exit 1
fi

echo "--------------------------------------"
echo "Found ${#BAM_FILES[@]} sorted BAM files"
echo "--------------------------------------"

for BAM in "${BAM_FILES[@]}"
do

    SAMPLE=$(basename "$BAM" .sorted.bam)
    COUNT_FILE="${COUNT_DIR}/${SAMPLE}.dexseq_counts.txt"

    echo ""
    echo "--------------------------------------"
    echo "Processing sample: $SAMPLE"
    echo "--------------------------------------"

    samtools quickcheck "$BAM" || {
        echo "ERROR: BAM file corrupted:"
        echo "$BAM"
        exit 1
    }

    if [[ -s "$COUNT_FILE" ]]; then
        echo "Count file already exists, skipping:"
        echo "$COUNT_FILE"
        continue
    fi

    python "$COUNT_SCRIPT" \
        -p yes \
        -r pos \
        -s reverse \
        -f bam \
        "$FLATTENED_GFF" \
        "$BAM" \
        "$COUNT_FILE" \
        > "${LOG_DIR}/${SAMPLE}.dexseq_count.stdout.log" \
        2> "${LOG_DIR}/${SAMPLE}.dexseq_count.stderr.log"

    [[ -s "$COUNT_FILE" ]] || {
        echo "ERROR: Failed to generate count file:"
        echo "$COUNT_FILE"
        exit 1
    }

    echo "Counts created:"
    echo "$COUNT_FILE"

done

echo ""
echo "======================================"
echo "DEXSeq counting complete"
echo "Output directory:"
echo "$COUNT_DIR"
echo "======================================"
