base  <- "/projects/CGS_shared/vfama/BRIGHT_PROJECT/DIU/"
lines <- c("BT483", "MCF7", "T47D", "BT549", "MDAMB231", "SUM159")
ranks <- read.table(paste0(base, "IsoQuantOutput/isoform_ranks.tsv"), header = TRUE)
k <- ranks[, c("Gene", "iso.id", "Isoform")]

switch_line <- function(line) {
  d <- paste0(base, line, "/TreeTraversal_Iso/DMSO_STM_10/")

  # geni significativi, con le due isoforme che cambiano di più
  res <- read.table(paste0(d, "DMSO_STM_25_results.csv"), sep = "\t", header = TRUE)
  need_res <- c("Gene", "FDR", "dPI", "maxDeltaPI_ix1", "maxDeltaPI_ix2")
  if (!all(need_res %in% colnames(res)))
    stop(line, ": colonne mancanti in results. Presenti: ", paste(colnames(res), collapse = ", "))
  sig <- res[res$FDR <= 0.05 & abs(res$dPI) >= 0.1, need_res]

  # pi1, pi2, delta per isoforma, dal Robj
  e <- new.env()
  load(paste0(d, "DMSO_STM_10X2_Atleast25reads.Robj"), envir = e)
  pdf <- e$processedDF
  need <- c("Gene", "IsoID", "pi1", "pi2", "delta")
  if (!all(need %in% colnames(pdf)))
    stop(line, ": colonne mancanti in processedDF. Presenti: ", paste(colnames(pdf), collapse = ", "))
  pdf <- pdf[, need]

  # formato lungo: una riga per isoforma (ix1 e ix2) per gene
  long <- rbind(
    data.frame(sig[, c("Gene", "FDR", "dPI")], rank_cambio = "ix1", iso.id = sig$maxDeltaPI_ix1),
    data.frame(sig[, c("Gene", "FDR", "dPI")], rank_cambio = "ix2", iso.id = sig$maxDeltaPI_ix2)
  )
  long <- merge(long, k,   by = c("Gene", "iso.id"), all.x = TRUE)
  long <- merge(long, pdf, by.x = c("Gene", "iso.id"), by.y = c("Gene", "IsoID"), all.x = TRUE)

  long$cell_line <- line
  # delta = pi1 - pi2, pi1 = DMSO -> delta > 0 = isoforma che cala con STM
  long$direzione <- ifelse(long$delta > 0, "giu_con_STM", "su_con_STM")
  long <- long[order(long$Gene, long$rank_cambio), ]
  long[, c("cell_line", "Gene", "rank_cambio", "Isoform", "iso.id",
           "pi1", "pi2", "delta", "direzione", "FDR", "dPI")]
}

switch_all <- lapply(lines, switch_line)
names(switch_all) <- lines
saveRDS(switch_all, file = "Scisorseq_ENST.rds")

