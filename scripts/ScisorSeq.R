#scisorseq DIU:
#!/usr/bin/env Rscript
#
# scisorseqr_output_production.R 
#
# Two outputs are written:
#   1. <out_prefix>.isoform_ranks.tsv          Gene  Isoform  n  iso.id
#   2. <out_prefix>.isoform_sample_counts.tsv  GeneID::isoID  sample  count
## ============================ CONFIG ===========================
library("data.table")
library("rtracklayer")
library("scisorseqr")
#Define output dir
out_prefix <- ""
#Import files
BASEDIR <- "/projects/CGS_shared/vfama/BRIGHT_PROJECT/"
gtf_path <- paste0(BASEDIR,"filtered_corrected_assembly.gtf")          # reference GTF used for quantification (.gz ok)
linear_counts_names <- list.files(path="/projects/CGS_shared/cugolini/BC_atlas/analysis/all_cell_lines/dRNA/quantification_IsoQuant/all_cell_lines/IsoQuant_BRIGHT/",pattern=".*transcript.*linear",recursive=TRUE)
linear_counts <- lapply(linear_counts_names, function(f) {
  d <- fread(paste0("/projects/CGS_shared/cugolini/BC_atlas/analysis/all_cell_lines/dRNA/quantification_IsoQuant/all_cell_lines/IsoQuant_BRIGHT/",f))
  setnames(d, sub("^#", "", names(d)))
  d
})
names(linear_counts) <- sub("/.*", "", linear_counts_names)
#Define classes
Luminal <- c("BT483", "MCF7", "T47D")
Basal  <- c("SUM159", "MDAMB231","BT549")
groups    <- unlist(lapply(linear_counts, function(i) unique(i$group_id)))
cell_line <- sub("_.*", "", groups)                  
treatment <- sub("^[^_]+_([^_]+)_.*", "\\1", groups) # DMSO / STM
sub_types <- list(Basal = Basal, Luminal = Luminal)
config <- do.call(rbind, lapply(names(sub_types), function(s) {
  sel <- cell_line %in% sub_types[[s]]
  data.frame(
    c1 = paste0(s, "DMSO"),
    g1 = paste(sort(groups[sel & treatment == "DMSO"]), collapse = ","),
    c2 = paste0(s, "STM"),
    g2 = paste(sort(groups[sel & treatment == "STM"]),  collapse = ",")
  )
}))

config_file <- file.path(BASEDIR, "config_DMSO_vs_STM_subtype.tsv")
write.table(config, config_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

keep_zero      <- FALSE   # keep isoforms whose summed count is 0
min_gene_count <- 0       # drop genes whose total summed count is below this
min_iso_count  <- 0       # drop isoforms whose summed count is below this
max_iso        <- 0       # keep only the top N isoforms per gene (0 = keep all)
dense          <- FALSE   # in file 2, also emit explicit 0 rows for absent isoform/sample pairs
as_float       <- FALSE   # write counts with 2 decimals instead of rounding to integers
on_missing     <- "drop"  # transcripts absent from the GTF: "drop" | "error" | "keep"
no_header      <- FALSE   # omit header lines
## ============================================================================
##  1. Read the GTF and build a transcript -> gene lookup
## ============================================================================
gtf <- import(gtf_path)
gtf <- gtf[gtf$type == "transcript" & !is.na(gtf$transcript_id) & !is.na(gtf$gene_id)]
t2g <- setNames(gtf$transcript_id, gtf$gene_id)
cat(length(t2g), "transcripts mapped to genes\n")
## ============================================================================
##  2. Aggregate all samples
## ============================================================================
long <- do.call(rbind, linear_counts)
long$count <- as.numeric(long$count)
long <- aggregate(count ~ feature_id + group_id, data = long, FUN = sum) # collapse any (tx, group_id) pair that appeared more than once
## ============================================================================
##  4. Sum counts per transcript, attach gene id, filter
## ============================================================================
tot <- aggregate(count ~ feature_id, data = long, FUN = sum)
names(tot) <- c("tx", "n")
if (!keep_zero)       tot <- tot[tot$n > 0, ]
if (min_iso_count > 0) tot <- tot[tot$n >= min_iso_count, ]
tot$gene <- names(t2g)[match(tot$tx, t2g)]
## ============================================================================
##  5. Rank isoforms within each gene (0 = most abundant), then filter further
## ============================================================================
tot <- tot[order(tot$gene, -tot$n, tot$tx), ]         # sort first ...
tot$iso.id <- ave(tot$n, tot$gene, FUN = seq_along) - 1  # ... then number within each gene
if (min_gene_count > 0) {
  gene_totals <- ave(tot$n, tot$gene, FUN = sum)
  tot <- tot[gene_totals >= min_gene_count, ]
}
if (max_iso > 0) tot <- tot[tot$iso.id < max_iso, ]
tot$pos <- seq_len(nrow(tot))              # remembers the final gene/rank order for file 2
tot$key <- paste0(tot$gene, "::", tot$iso.id)
## ============================================================================
##  6. Write file 1: Gene  Isoform  n  iso.id
## ============================================================================
n_out <- if (as_float) sprintf("%.2f", tot$n) else as.character(round(tot$n))
file1 <- data.frame(Gene = tot$gene, Isoform = tot$tx, n = n_out, iso.id = tot$iso.id)
dir.create(out_prefix, recursive = TRUE, showWarnings = FALSE)
write.table(file1, paste0(out_prefix, "isoform_ranks.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = !no_header)
## ============================================================================
##  7. Write file 2: GeneID::isoID  sample  count
## ============================================================================
if (dense) {
  grid <- expand.grid(feature_id = tot$tx, group_id = groups, stringsAsFactors = FALSE)
  grid <- merge(grid, long, by = c("feature_id", "group_id"), all.x = TRUE)
  grid$count[is.na(grid$count)] <- 0
} else {
  grid <- long[long$feature_id %in% tot$tx, ]
}
grid$pos   <- tot$pos[match(grid$feature_id, tot$tx)]
grid$key   <- tot$key[match(grid$feature_id, tot$tx)]
grid$lorder <- match(grid$group_id, groups)
grid <- grid[order(grid$pos, grid$lorder), ]
count_out <- if (as_float) sprintf("%.2f", grid$count) else as.character(round(grid$count))
file2 <- data.frame(iso = grid$key, sample = grid$group_id, count = count_out)
write.table(file2, paste0(out_prefix, "NumIsoPerCluster"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
## ============================================================================
cat("[done]", nrow(tot), "isoforms across", length(unique(tot$gene)), "genes and", length(labels), "samples\n")
cat("[done] wrote", paste0(out_prefix, "isoform_ranks.tsv"), "\n")
cat("[done] wrote", paste0(out_prefix, "NumIsoPerCluster"), "\n")
## ---- write the config file for DiffSplicingAnalysis -------------------------
## group names should not contain underscores
config <- do.call(rbind, lapply(names(sub_types), function(s) {
  sel <- cell_line %in% sub_types[[s]]
  data.frame(
    c1 = paste0(s, "DMSO"),
    g1 = paste(sort(groups[sel & treatment == "DMSO"]), collapse = ","),
    c2 = paste0(s, "STM"),
    g2 = paste(sort(groups[sel & treatment == "STM"]),  collapse = ",")
  )
}))

config_file <- file.path(BASEDIR, "config_DMSO_vs_STM_subtype.tsv")
write.table(config, config_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

#####################
# devtools::install_github('noush-joglekar/scisorseqr',build_vignettes = FALSE)
DiffSplicingAnalysis(config_file,
                     numIsoforms = 10, minNumReads = 25,
                     typeOfTest = 'Iso', numThreads = 4)
writeLines(unique(c(config$c1, config$c2)), "groups.txt")
triHeatmap("TreeTraversal_Iso/", "groups.txt",
           typeOfTest = "Iso", outName = "Heatmap_DMSO_vs_STM_subtype")
 