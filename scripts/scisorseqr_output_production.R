#!/usr/bin/env Rscript
#
# scisorseqr_output_production_DMSO_vs_STORM.R
#
# Stessi output della versione originale (isoform_ranks.tsv / NumIsoPerCluster),
# generati una sola volta su TUTTI i campioni. Il confronto DMSO vs STORM viene
# poi eseguito separatamente per ciascuna cell line (6 analisi indipendenti).

## ============================ CONFIG ===========================

BASEDIR <- "/projects/CGS_shared/vfama/BRIGHT_PROJECT/"

gtf_path       <- paste0(BASEDIR,"filtered_corrected_assembly.gtf")
manifest_path  <- paste0(BASEDIR,"DIU/path_sample_bright.tsv")  # deve includere anche i campioni STORM

out_prefix     <- paste0(BASEDIR,"DIU/IsoQuantOutput/")
work_dir       <- paste0(BASEDIR,"DIU/")

## campioni DMSO/STM per cell line (group names senza underscore per scisorseqr)
## ADATTA i nomi/numero di repliche ai tuoi campioni reali
sample_map <- list(
  BT483 = list(
    DMSO  = c("BT483_DMSO_1","BT483_DMSO_2","BT483_DMSO_3","BT483_DMSO_4"),
    STM = c("BT483_STM_1","BT483_STM_2","BT483_STM_3","BT483_STM_4")
  ),
  MCF7 = list(
    DMSO  = c("MCF7_DMSO_1","MCF7_DMSO_2","MCF7_DMSO_3","MCF7_DMSO_4"),
    STM = c("MCF7_STM_1","MCF7_STM_2","MCF7_STM_3","MCF7_STM_4")
  ),
  T47D = list(
    DMSO  = c("T47D_DMSO_1","T47D_DMSO_2","T47D_DMSO_4"),
    STM = c("T47D_STM_1","T47D_STM_2","T47D_STM_4")
  ),
  BT549 = list(
    DMSO  = c("BT549_DMSO_1","BT549_DMSO_2","BT549_DMSO_3","BT549_DMSO_4"),
    STM = c("BT549_STM_1","BT549_STM_2","BT549_STM_3","BT549_STM_4")
  ),
  MDAMB231 = list(
    DMSO  = c("MDAMB231_DMSO_1","MDAMB231_DMSO_2","MDAMB231_DMSO_3","MDAMB231_DMSO_4"),
    STM = c("MDAMB231_STM_1","MDAMB231_STM_2","MDAMB231_STM_3","MDAMB231_STM_4")
  ),
  SUM159 = list(
    DMSO  = c("SUM159_DMSO_1","SUM159_DMSO_2","SUM159_DMSO_3","SUM159_DMSO_4"),
    STM = c("SUM159_STM_1","SUM159_STM_2","SUM159_STM_3","SUM159_STM_4")
  )
)

keep_zero      <- FALSE
min_gene_count <- 0
min_iso_count  <- 0
max_iso        <- 0
dense          <- FALSE
as_float       <- FALSE
on_missing     <- "drop"
no_header      <- FALSE

## ============================================================================
##  1. GTF -> transcript -> gene lookup   (invariato)
## ============================================================================

gtf <- read.delim(gtf_path, header = FALSE, comment.char = "#", quote = "",
                  stringsAsFactors = FALSE, fill = TRUE)
gtf <- gtf[gtf$V3 %in% c("transcript", "exon"), ]
if (any(gtf$V3 == "transcript")) gtf <- gtf[gtf$V3 == "transcript", ]

has_ids <- grepl('transcript_id "', gtf$V9, fixed = TRUE) & grepl('gene_id "', gtf$V9, fixed = TRUE)
gtf <- gtf[has_ids, ]

tx_id   <- sub('.*transcript_id "([^"]+)".*', "\\1", gtf$V9)
gene_id <- sub('.*gene_id "([^"]+)".*',       "\\1", gtf$V9)

keep <- !duplicated(tx_id)
t2g <- setNames(gene_id[keep], tx_id[keep])
cat(length(t2g), "transcripts mapped to genes\n")

## ============================================================================
##  2. Manifest   (invariato — deve contenere anche i campioni STM)
## ============================================================================

manifest_lines <- readLines(manifest_path)
manifest_lines <- trimws(manifest_lines)
manifest_lines <- manifest_lines[nzchar(manifest_lines) & !startsWith(manifest_lines, "#")]
manifest_split <- strsplit(manifest_lines, "[\t ]+")

manifest <- data.frame(
  path   = sapply(manifest_split, function(x) x[1]),
  sample = sapply(manifest_split, function(x) x[2]),
  label  = sapply(manifest_split, function(x) if (length(x) > 2) x[3] else NA),
  stringsAsFactors = FALSE
)

## ============================================================================
##  3. Lettura conteggi per ogni (file, sample) del manifest   (invariato)
## ============================================================================

long_parts <- list()
labels <- c()

for (p in unique(manifest$path)) {

  d <- read.delim(p, header = FALSE, comment.char = "#", quote = "",
                  stringsAsFactors = FALSE, fill = TRUE)
  colnames(d)[1:3] <- c("tx", "sample", "count")

  rows <- manifest[manifest$path == p, ]

  for (i in seq_len(nrow(rows))) {
    samp  <- rows$sample[i]
    label <- rows$label[i]

    if (samp == "*") {
      these_samples <- unique(d$sample)
      these_labels  <- these_samples
    } else {
      these_samples <- samp
      these_labels  <- if (is.na(label) || label == "") samp else label
    }

    for (j in seq_along(these_samples)) {
      s   <- these_samples[j]
      lab <- these_labels[j]

      if (lab %in% labels) stop("duplicate output sample label '", lab, "' (from ", p,
                                "). Use the manifest's 3rd column to rename it.")

      sub <- d[d$sample == s, c("tx", "count")]
      if (nrow(sub) == 0) {
        cat("[warn] sample '", s, "' produced 0 rows from ", p, "\n", sep = "")
        next
      }
      sub$label <- lab
      long_parts[[length(long_parts) + 1]] <- sub
      labels <- c(labels, lab)
      cat("[read]", p, ":", s, "->", lab, "(", nrow(sub), "rows )\n")
    }
  }
}

if (length(long_parts) == 0) stop("no counts were read -- check the sample names in the manifest")

long <- do.call(rbind, long_parts)
long$count <- as.numeric(long$count)
long <- aggregate(count ~ tx + label, data = long, FUN = sum)

## ============================================================================
##  4. Somma per trascritto, gene id, filtri, ranking   (invariato)
## ============================================================================

tot <- aggregate(count ~ tx, data = long, FUN = sum)
names(tot) <- c("tx", "n")

if (!keep_zero)       tot <- tot[tot$n > 0, ]
if (min_iso_count > 0) tot <- tot[tot$n >= min_iso_count, ]

tot$gene <- unname(t2g[tot$tx])
missing_tx <- is.na(tot$gene)

if (any(missing_tx)) {
  cat("[warn]", sum(missing_tx), "transcripts were not found in the GTF (action:", on_missing, "); e.g.",
      paste(head(tot$tx[missing_tx], 5), collapse = ", "), "\n")
  if (on_missing == "error") stop("transcripts missing from the GTF")
  if (on_missing == "drop")  tot <- tot[!missing_tx, ]
  if (on_missing == "keep")  tot$gene[missing_tx] <- tot$tx[missing_tx]
}

tot <- tot[order(tot$gene, -tot$n, tot$tx), ]
tot$iso.id <- ave(tot$n, tot$gene, FUN = seq_along) - 1

if (min_gene_count > 0) {
  gene_totals <- ave(tot$n, tot$gene, FUN = sum)
  tot <- tot[gene_totals >= min_gene_count, ]
}
if (max_iso > 0) tot <- tot[tot$iso.id < max_iso, ]

tot$pos <- seq_len(nrow(tot))
tot$key <- paste0(tot$gene, "::", tot$iso.id)

## ============================================================================
##  5. File 1: Gene  Isoform  n  iso.id   (invariato)
## ============================================================================

n_out <- if (as_float) sprintf("%.2f", tot$n) else as.character(round(tot$n))
file1 <- data.frame(Gene = tot$gene, Isoform = tot$tx, n = n_out, iso.id = tot$iso.id)

dir.create(out_prefix, recursive = TRUE, showWarnings = FALSE)
write.table(file1, paste0(out_prefix, "isoform_ranks.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = !no_header)

## ============================================================================
##  6. File 2: GeneID::isoID  sample  count   (invariato)
## ============================================================================

if (dense) {
  grid <- expand.grid(tx = tot$tx, label = labels, stringsAsFactors = FALSE)
  grid <- merge(grid, long, by = c("tx", "label"), all.x = TRUE)
  grid$count[is.na(grid$count)] <- 0
} else {
  grid <- long[long$tx %in% tot$tx, ]
}

grid$pos   <- tot$pos[match(grid$tx, tot$tx)]
grid$key   <- tot$key[match(grid$tx, tot$tx)]
grid$lorder <- match(grid$label, labels)
grid <- grid[order(grid$pos, grid$lorder), ]

count_out <- if (as_float) sprintf("%.2f", grid$count) else as.character(round(grid$count))
file2 <- data.frame(iso = grid$key, sample = grid$label, count = count_out)

write.table(file2, paste0(out_prefix, "NumIsoPerCluster"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("[done]", nrow(tot), "isoforms across", length(unique(tot$gene)), "genes and", length(labels), "samples\n")

## ============================================================================
##  7. Config + analisi DMSO vs STM, per ciascuna cell line
## ============================================================================

library(scisorseqr)

for (cline in names(sample_map)) {

  cline_dir <- file.path(work_dir, cline)
  dir.create(cline_dir, recursive = TRUE, showWarnings = FALSE)

  dmso_samples  <- sample_map[[cline]]$DMSO
  storm_samples <- sample_map[[cline]]$STM

  ## verifica che i campioni siano tra le label generate dal manifest
  missing_lab <- setdiff(c(dmso_samples, storm_samples), labels)
  if (length(missing_lab)) {
    stop("Campioni non trovati tra le label lette dal manifest per ", cline, ": ",
        paste(missing_lab, collapse = ", "))
  }

  ## config.tsv con un solo confronto: DMSO vs STORM per questa linea
  config_line <- paste("DMSO", paste(dmso_samples, collapse = ","),
                       "STM", paste(storm_samples, collapse = ","),
                       sep = "\t")
  config_out_cline <- file.path(cline_dir, "config.tsv")
  writeLines(config_line, config_out_cline)
  cat("[done] wrote config for", cline, "->", config_out_cline, "\n")

  ## symlink a IsoQuantOutput condiviso, cosi' DiffSplicingAnalysis lo trova
  ## relativo alla working dir della cell line
  link_path <- file.path(cline_dir, "IsoQuantOutput")
  if (!file.exists(link_path)) {
    file.symlink(normalizePath(out_prefix, mustWork = TRUE), link_path)
  }

  setwd(cline_dir)

  DiffSplicingAnalysis(config_out_cline,
                       numIsoforms = 10, minNumReads = 25,
                       typeOfTest = 'Iso', numThreads = 4)

  writeLines(c("DMSO", "STM"), file.path(cline_dir, "groups.txt"))

  triHeatmap("TreeTraversal_Iso/", "groups.txt", typeOfTest = "Iso",
            outName = paste0("Heatmap_", cline, "_DMSO_vs_STM"))

  setwd(work_dir)
}

cat("[done] DIU DMSO vs STM completata per tutte le cell line\n")