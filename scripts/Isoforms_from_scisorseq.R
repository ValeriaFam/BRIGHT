base  <- "/projects/CGS_shared/vfama/BRIGHT_PROJECT/DIU/"
lines <- c("BT483", "MCF7", "T47D", "BT549", "MDAMB231", "SUM159")
ranks <- read.table(paste0(base, "IsoQuantOutput/isoform_ranks.tsv"), header = TRUE)
k <- ranks[, c("Gene", "iso.id", "Isoform")]

switch_line <- function(line, thr = 0.1) {
  d <- paste0(base, line, "/TreeTraversal_Iso/DMSO_STM_10/")

  # geni significativi (stessi criteri dello script di Fatemeh)
  res <- read.table(paste0(d, "DMSO_STM_25_results.csv"), sep = "\t", header = TRUE)
  sig <- res[res$FDR <= 0.05 & abs(res$dPI) >= 0.1, c("Gene", "FDR", "dPI")]

  # tutte le isoforme testate, dal Robj
  e <- new.env()
  load(paste0(d, "DMSO_STM_10X2_Atleast25reads.Robj"), envir = e)
  pdf <- e$processedDF
  need <- c("Gene", "IsoID", "pi1", "pi2", "delta")
  if (!all(need %in% colnames(pdf)))
    stop(line, ": colonne mancanti. Presenti: ", paste(colnames(pdf), collapse = ", "))
  pdf <- pdf[, need]

  # solo geni significativi e isoforme con |delta| >= soglia
  x <- pdf[pdf$Gene %in% sig$Gene & !is.na(pdf$delta) & abs(pdf$delta) >= thr, ]
  x <- merge(x, k, by.x = c("Gene", "IsoID"), by.y = c("Gene", "iso.id"), all.x = TRUE)
  x <- merge(x, sig, by = "Gene")
  x$cell_line <- line
  # direzione: delta = pi1 - pi2, pi1 = DMSO -> delta > 0 = isoforma che cala con STM
  x$direzione <- ifelse(x$delta > 0, "giu_con_STM", "su_con_STM")
  x[, c("cell_line", "Gene", "Isoform", "IsoID", "pi1", "pi2", "delta", "direzione", "FDR", "dPI")]
}

switch_all <- lapply(lines, switch_line)

