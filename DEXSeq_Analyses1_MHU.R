# =====================================
# HAoSMC KD 
# =====================================
suppressPackageStartupMessages({
  library(DEXSeq)
  library(BiocParallel)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

# =====================================
# Paths (Mac)
# =====================================

project_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_20th_Feb/DEXSeq/HAoSMC_Ensembl112_/DEXSeq_HAoSMCs_KD_112/counts"
count_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_20th_Feb/DEXSeq/HAoSMC_Ensembl112_/DEXSeq_HAoSMCs_KD_112/counts" 
flattened_file <- "/Users/maryamuhabilausman/Desktop/Analysis_20th_Feb/DEXSeq/HAoSMC_Ensembl112_/DEXSeq_HAoSMCs_KD_112/annotation/Homo_sapiens.GRCh38.112.DEXSeq.gff"
results_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_20th_Feb/DEXSeq/HAoSMC_Ensembl112_"

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# =====================================
# Sample table
# =====================================

sampleTable <- data.frame(
  row.names = c("CTR1_Z20_T","CTR2_Z30_T","Z350Ctr","Z950Ctr","KD1_Z20_T","KD2_Z30_T","Z650KD","Z1250KD"),
  condition = c("control","control","control","control",
                "KD","KD","KD","KD"),
  stringsAsFactors = FALSE
)

sampleTable$condition <- factor(sampleTable$condition, levels = c("control", "KD"))

# =====================================
# Input files
# =====================================

countFiles <- file.path(count_dir, paste0(rownames(sampleTable), ".dexseq_counts.txt"))

missing_files <- countFiles[!file.exists(countFiles)]
if (length(missing_files) > 0) {
  stop("Missing count files:\n", paste(missing_files, collapse = "\n"))
}

if (!file.exists(flattened_file)) {
  stop("Flattened DEXSeq GFF not found:\n", flattened_file)
}

# =====================================
# Clean count files
# =====================================

clean_dir <- file.path(project_dir, "counts_clean")
dir.create(clean_dir, showWarnings = FALSE)

clean_dexseq_count_file <- function(infile, outfile) {
  x <- readLines(infile)
  x <- x[!grepl("^_", x)]
  x <- x[nzchar(x)]
  writeLines(x, outfile)
}

cleanFiles <- file.path(clean_dir, basename(countFiles))

for (i in seq_along(countFiles)) {
  clean_dexseq_count_file(countFiles[i], cleanFiles[i])
}

# =====================================
# Read cleaned count files manually
# =====================================

read_dexseq_counts <- function(f) {
  read.delim(
    f,
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    col.names = c("count_id", "count")
  )
}

count_list <- lapply(cleanFiles, read_dexseq_counts)

ref_ids <- count_list[[1]]$count_id
same_ids <- vapply(count_list, function(x) identical(x$count_id, ref_ids), logical(1))

if (!all(same_ids)) {
  stop("Count files do not have identical exon-bin IDs/order:\n",
       paste(basename(cleanFiles)[!same_ids], collapse = "\n"))
}

count_matrix <- do.call(cbind, lapply(count_list, `[[`, "count"))
mode(count_matrix) <- "integer"
rownames(count_matrix) <- ref_ids
colnames(count_matrix) <- rownames(sampleTable)

# Parse IDs from count file
# Example raw ID:  "ENSG00000293600+ENSG00000176087":"009"
count_groupID <- sub(':"[^"]+"$', "", ref_ids)
count_groupID <- gsub('^"|"$', "", count_groupID)

count_featureID <- sub('^.*:"([^"]+)"$', "\\1", ref_ids)

count_key <- paste(count_groupID, count_featureID, sep = ":")

# =====================================
# Import flattened GFF and build feature ranges
# =====================================

gff <- rtracklayer::import(flattened_file)

# Keep only exonic parts
if ("type" %in% names(mcols(gff))) {
  exonic_parts <- gff[mcols(gff)$type == "exonic_part"]
} else {
  stop("Could not find 'type' column in flattened GFF metadata.")
}

gff_cols <- names(mcols(exonic_parts))

if (!("gene_id" %in% gff_cols)) {
  stop("Flattened GFF does not contain a gene_id column.")
}

# exonic part number column name can vary slightly
part_col <- NULL
for (nm in c("exonic_part_number", "exon_number", "exonicPart")) {
  if (nm %in% gff_cols) {
    part_col <- nm
    break
  }
}

if (is.null(part_col)) {
  stop("Could not find exonic part number column in flattened GFF. Available columns:\n",
       paste(gff_cols, collapse = ", "))
}

gff_groupID <- as.character(mcols(exonic_parts)$gene_id)
gff_featureID <- sprintf("%03d", as.integer(mcols(exonic_parts)[[part_col]]))
gff_key <- paste(gff_groupID, gff_featureID, sep = ":")

match_idx <- match(count_key, gff_key)

if (any(is.na(match_idx))) {
  bad <- count_key[is.na(match_idx)]
  stop(
    "Some count IDs could not be matched to the flattened GFF.\nFirst few unmatched IDs:\n",
    paste(head(bad, 10), collapse = "\n")
  )
}

featureRanges <- exonic_parts[match_idx]

# =====================================
# Build DEXSeq dataset directly
# =====================================

dxd <- DEXSeqDataSet(
  countData = count_matrix,
  sampleData = sampleTable,
  design = ~ sample + exon + condition:exon,
  featureID = count_featureID,
  groupID = count_groupID,
  featureRanges = featureRanges
)

# =====================================
# Run full DEXSeq analysis
# =====================================

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")

dxr <- DEXSeqResults(dxd)

saveRDS(dxd, file = file.path(results_dir, "HAoSMC_DEXSeq_dxd.rds"))
saveRDS(dxr, file = file.path(results_dir, "HAoSMC_DEXSeq_dxr.rds"))

# =====================================
# Full annotated results
# =====================================

# =====================================
# Full annotated results
# =====================================

res_df <- as.data.frame(dxr)

rr <- rowRanges(dxd)

if (nrow(res_df) != length(rr)) {
  stop("Number of result rows does not match number of genomic ranges.")
}

coord_df <- data.frame(
  featureID = count_featureID,
  groupID = count_groupID,
  exonbase = count_featureID,
  chr = as.character(seqnames(rr)),
  start = start(rr),
  end = end(rr),
  strand = as.character(strand(rr)),
  stringsAsFactors = FALSE
)

# bind coordinates + DEXSeq statistics by row order
res_annot <- cbind(
  coord_df,
  res_df[, setdiff(colnames(res_df), c("featureID", "groupID")), drop = FALSE]
)

write.csv(
  res_annot,
  file = file.path(results_dir, "HAoSMC_DEXSeq_results_annotated.csv"),
  row.names = FALSE
)

write.csv(
  subset(res_annot, !is.na(padj) & padj < 0.05),
  file = file.path(results_dir, "HAoSMC_DEXSeq_significant_FDR0.05.csv"),
  row.names = FALSE
)
# =====================================
# Global plots
# =====================================

pdf(file.path(results_dir, "HAoSMC_DEXSeq_dispersion_plot.pdf"))
plotDispEsts(dxd)
dev.off()

pdf(file.path(results_dir, "HAoSMC_DEXSeq_MAplot.pdf"))
plotMA(dxr)
dev.off()

res_ordered <- res_annot[!is.na(res_annot$padj), ]
res_ordered <- res_ordered[order(res_ordered$padj), ]

if (nrow(res_ordered) > 0) {
  top_gene <- as.character(res_ordered$groupID[1])
  pdf(file.path(results_dir, "HAoSMC_DEXSeq_top_gene_plot.pdf"))
  plotDEXSeq(dxr, top_gene, legend = TRUE)
  dev.off()
}


# ====ROI USED =================================
# ROI extraction and plotting
# =====================================

library(GenomicRanges)

rr <- rowRanges(dxd)

# detect chromosome naming style automatically
chr_names <- unique(as.character(seqnames(rr)))

roi_chr <- if ("chr12" %in% chr_names) "chr12" else "12"


roi_start <- xxxxxx
roi_end <- xxxxxx


roi_gr <- GRanges(
  seqnames = roi_chr,
  ranges = IRanges(start = roi_start, end = roi_end)
)

hits <- findOverlaps(rr, roi_gr)

if (length(hits) == 0) {
  stop("No DEXSeq exon bins overlap the ROI.")
}

# -------------------------------------
# Build full annotated results
# -------------------------------------

res_df <- as.data.frame(dxr)

coord_df <- data.frame(
  featureID = count_featureID,
  groupID = count_groupID,
  exonbase = count_featureID,
  chr = as.character(seqnames(rr)),
  start = start(rr),
  end = end(rr),
  strand = as.character(strand(rr)),
  stringsAsFactors = FALSE
)

full_res <- cbind(
  coord_df,
  res_df[, setdiff(colnames(res_df), c("featureID","groupID")), drop = FALSE]
)

roi_res <- full_res[queryHits(hits), , drop = FALSE]

# -------------------------------------
# Save ROI result tables
# -------------------------------------

write.csv(
  roi_res,
  file = file.path(results_dir, "ROI_chr12_wholelocus_results.csv"),
  row.names = FALSE
)

roi_sig <- subset(roi_res, !is.na(padj) & padj < 0.05)

write.csv(
  roi_sig,
  file = file.path(results_dir, "ROI_wholelocus_significant_FDR0.05.csv"),
  row.names = FALSE
)

# -------------------------------------
# Save counts for ROI bins
# -------------------------------------

roi_counts <- counts(dxd)[queryHits(hits), , drop = FALSE]

roi_counts_df <- cbind(
  roi_res[, c("featureID","groupID","chr","start","end","strand","exonbase")],
  as.data.frame(roi_counts)
)

write.csv(
  roi_counts_df,
  file = file.path(results_dir, "ROI_chr12_wholelocus_counts.csv"),
  row.names = FALSE
)

# -------------------------------------
# Plot DEXSeq exon usage for genes
# -------------------------------------

roi_genes <- unique(roi_res$groupID)

plot_dir <- file.path(results_dir, "ROI_chr12_wholelocus_gene_plots")
dir.create(plot_dir, showWarnings = FALSE)

for (g in roi_genes) {
  
  pdf(file.path(plot_dir, paste0(g, "_DEXSeqwholelocus_plot.pdf")))
  try(plotDEXSeq(dxr, g, legend = TRUE), silent = TRUE)
  dev.off()
  
}

cat("\nROI extraction complete\n")
cat("Genes overlapping ROI:\n")
print(roi_genes)

cat("\nResults saved in:\n")
cat(results_dir, "\n")

















