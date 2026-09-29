# ============================================================================================
# Script:  DRIMSeq_input_generation.R
# Purpose: Creating the exact input format for DRIMSeq using isoquant output files and GTF
# Author:  Valeria Famà
# ============================================================================================

#######################################
###### import paths and files #########
#######################################
library(data.table)
BASEDIR <- "/projects/CGS_shared/vfama/BRIGHT_PROJECT/"

gtf_path <- paste0(BASEDIR,"/filtered_corrected_assembly.gtf")  # reference GTF used for quantification
output_path  <- paste0(BASEDIR,"/DIU") 

quantif_files_path <- "/projects/CGS_shared/cugolini/BC_atlas/analysis/all_cell_lines/dRNA/quantification_IsoQuant/all_cell_lines/IsoQuant_BRIGHT"

##########################################
### Assign Tx IDs to Gene IDs using GTF###
##########################################

gtf <- read.delim(gtf_path, header = FALSE, comment.char = "#", quote = "",
                  stringsAsFactors = FALSE, fill = TRUE)

gtf <- gtf[gtf$V3 %in% c("transcript", "exon"), ]
if (any(gtf$V3 == "transcript")) gtf <- gtf[gtf$V3 == "transcript", ]  # prefer transcript lines if present

has_ids <- grepl('transcript_id "', gtf$V9, fixed = TRUE) & grepl('gene_id "', gtf$V9, fixed = TRUE)
gtf <- gtf[has_ids, ]

tx_id   <- sub('.*transcript_id "([^"]+)".*', "\\1", gtf$V9)
gene_id <- sub('.*gene_id "([^"]+)".*',       "\\1", gtf$V9)

valid <- !is.na(tx_id) & !is.na(gene_id)

mapping <- unique(data.frame(
  transcript_id = tx_id[valid],
  gene_id = gene_id[valid]
))

conflicts <- mapping$transcript_id[
  duplicated(mapping$transcript_id) |
    duplicated(mapping$transcript_id, fromLast = TRUE)
]

if (length(conflicts)) {
  stop(
    "Some transcripts map to multiple genes: ",
    paste(head(unique(conflicts), 10), collapse = ", ")
  )
}

cat(length(unique(mapping$transcript_id)), "transcripts mapped to genes\n")

##################################################
###### Merge Transcript counts with GeneIDs ######
##################################################

quantif_files_path_all <- list.dirs(path = quantif_files_path, recursive = TRUE, full.names = TRUE)

cell_lines <- c("BT483", "MCF7", "T47D", "MDAMB231", "SUM159", "BT549")

genera_input_DIU <- function(cline, quantif_files_path_all, mapping, BASEDIR) {

  # cartelle di questa cell line, condizioni DMSO e STORM
  pattern <- paste0("^", cline, "_(DMSO|STM)$")
  cline_dirs <- quantif_files_path_all[grepl(pattern, basename(quantif_files_path_all))]

  if (length(cline_dirs) != 2) {
    stop("Attese 2 cartelle (DMSO+STM) per ", cline, ", trovate: ",
        paste(basename(cline_dirs), collapse = ", "))
  }

  tx_counts_files <- file.path(cline_dirs, paste0(basename(cline_dirs), ".transcript_grouped_counts.tsv"))
  tx_counts <- setNames(lapply(tx_counts_files, read.delim, check.names = FALSE), basename(cline_dirs))

  merged_counts <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), tx_counts)
  merged_counts[is.na(merged_counts)] <- 0
  colnames(merged_counts) <- c("feature_id", setdiff(colnames(merged_counts), "gene_id"))

  merged_counts_txID_geneID <- merge(merged_counts, mapping, by.x = "feature_id", by.y = "transcript_id")
  sample_names <- setdiff(colnames(merged_counts_txID_geneID), c("gene_id", "feature_id"))
  merged_counts_txID_geneID <- merged_counts_txID_geneID[, c("gene_id", "feature_id", sample_names)]

  dir.create(file.path(BASEDIR, "DIU", cline), recursive = TRUE, showWarnings = FALSE)
  write.table(merged_counts_txID_geneID,
              file.path(BASEDIR, "DIU", cline, paste0(cline, "_DMSO_STM_Tx_counts.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # metadata: gruppo = DMSO/STORM, estratto dal nome campione
  metadata <- data.frame(
    sample_id = sample_names,
    group = sub("^.*_(DMSO|STM)_[0-9]+$", "\\1", sample_names)
  )
  write.table(metadata,
              file.path(BASEDIR, "DIU", cline, paste0(cline, "_DRIMSEQ_metadata.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  invisible(list(counts = merged_counts_txID_geneID, metadata = metadata))
}

for (cline in cell_lines) {
  genera_input_DIU(cline, quantif_files_path_all, mapping, BASEDIR)
}







 