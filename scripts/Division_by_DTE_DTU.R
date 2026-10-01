setwd("/projects/CGS_shared/vfama/BRIGHT_PROJECT/")
library("dplyr")
library("ggExtra")
library("ggrastr")
library("ggplot2")
library("pheatmap")
library("DESeq2")
library("ggrepel")
library("tidyr")
library("rtracklayer")
library("gridExtra")
library("grid")
library("svglite")
library("irr")
library("tidyverse")
library("txdbmaker")   
library("GenomicFeatures")
library("RColorBrewer")
library("purrr")
library("ComplexUpset")
library("patchwork")
library("clusterProfiler")
library("org.Hs.eg.db")

cell_line <- c("MCF7","BT483","T47D","SUM159","MDAMB231","BT549")
palette <- c(MCF7="#CC6677",BT483="#882255",T47D="#AA4499",SUM159="#117733",MDAMB231="#999933",BT549="#44AA99")

#Read the tables
names_tab <- list.files(path=".",pattern="sites_per_transcript_and_features.csv",recursive=TRUE)
tab <- lapply(names_tab,function(i){
		read.csv(i,header=TRUE)
	})
names(tab) <- gsub(".*/([A-Za-z0-9]+)_.*\\.csv$", "\\1", names_tab)

get_transcripts <- function(df, status) {
  unique(df$chrom[df$regulation == status])
}
#Extract up and down regulated transcripts
up_tr   <- lapply(tab[cell_line], get_transcripts, status = "UP")
down_tr <- lapply(tab[cell_line], get_transcripts, status = "DOWN")

Extract_commons_transcripts <- function(sets, luminal, basal, min_n = 5) {
  all_tr <- unique(unlist(sets))
  # matrice presenza/assenza gene x cell line
  m <- sapply(sets, function(s) all_tr%in% s)
  rownames(m) <- all_tr

  n_tot <- rowSums(m)
  n_lum <- rowSums(m[, luminal, drop = FALSE])
  n_bas <- rowSums(m[, basal,   drop = FALSE])

  list(
    # comuni ad almeno 5/6 cell lines
    almeno_5_su_6 = all_tr[n_tot >= min_n],

    # comuni a TUTTE le basal e in NESSUNA luminal (specifici basal)
    basal_only    = all_tr[n_bas >= (length(basal)-1) & n_lum == 0],

    # comuni a TUTTE le luminal e in NESSUNA basal (specifici basal)
    luminal_only  = all_tr[n_lum >= (length(luminal)-1) & n_bas == 0],

    matrice = m
  )
}

luminal <- cell_line[1:3]
basal <- cell_line[4:6]

res_up   <- Extract_commons_transcripts(up_tr,   luminal, basal)
res_down <- Extract_commons_transcripts(down_tr, luminal, basal)

#Summary
recap <- data.frame(
  categoria = c("UP", "DOWN"),
  n_almeno_5_su_6 = c(length(res_up$almeno_5_su_6), length(res_down$almeno_5_su_6)),
  basal_only      = c(length(res_up$basal_only),    length(res_down$basal_only)),
  luminal_only    = c(length(res_up$luminal_only),  length(res_down$luminal_only))
)
print(recap)

#How many DTE are also DGE? (quanti trascritti nella mia classe appartengono a un gene disregolato)
conta_gene_disregolati <- function(res, tab, cell_lines, min_5_6 = 5, min_2_3 = 2) {

  classi <- list(
    almeno_5_su_6 = list(ids = res$almeno_5_su_6, min_n = min_5_6),
    basal_only    = list(ids = res$basal_only,    min_n = min_2_3),
    luminal_only  = list(ids = res$luminal_only,  min_n = min_2_3)
  )

  risultati <- lapply(names(classi), function(k) {
    ids   <- classi[[k]]$ids
    min_n <- classi[[k]]$min_n

    long <- do.call(rbind, lapply(cell_lines, function(cl) {
      presenti <- ids[ids %in% rownames(res$matrice)[res$matrice[, cl]]]
      if (length(presenti) == 0) return(NULL)
      df <- tab[[cl]]
      sub <- df[df$chrom %in% presenti, c("chrom", "gene_id", "status_genes")]
      if (nrow(sub) == 0) return(NULL)
      sub$cell_line <- cl
      sub
    }))

    n_disreg <- tapply(long$status_genes %in% c("UP", "DOWN"), long$chrom, sum)
    transcript_disregolati <- names(n_disreg)[n_disreg >= min_n]

    gene_map <- unique(long[, c("chrom", "gene_id")])
    gene_disregolati <- unique(gene_map$gene_id[gene_map$chrom %in% transcript_disregolati])
    gene_totali       <- unique(gene_map$gene_id[gene_map$chrom %in% ids])

    list(
      summary = data.frame(
        classe = k, soglia_min = min_n,
        n_transcript_totali = length(ids),
        n_transcript_con_gene_disreg = length(transcript_disregolati),
        n_gene_totali = length(gene_totali),
        n_gene_disregolati = length(gene_disregolati),
        perc_gene_disregolati = round(100 * length(gene_disregolati) / length(gene_totali), 1)
      ),
      gene_disregolati = gene_disregolati,
      gene_totali       = gene_totali   # utile come "universe" per la GO
    )
  })
  names(risultati) <- names(classi)
  risultati
}

cell_lines <- dimnames(res_up$matrice)[[2]]

risultati_down <- conta_gene_disregolati(res_down, tab, cell_lines)
risultati_up   <- conta_gene_disregolati(res_up,   tab, cell_lines)

do.call(rbind, lapply(risultati_down, `[[`, "summary"))   # riepilogo numerico
library(clusterProfiler)
library(org.Hs.eg.db)

## universe comune a tutte le analisi: tutti i geni testati, in qualsiasi cell line
universe_tutti <- unique(unlist(lapply(tab, function(df) df$gene_id)))
cat(length(universe_tutti), "geni nell'universe\n")

do_GO <- function(gene_ids, universe_ids = universe_tutti, ont = "BP") {
  gene_ids_clean     <- sub("\\..*", "", gene_ids)
  universe_ids_clean <- sub("\\..*", "", universe_ids)

  enrichGO(
    gene          = gene_ids_clean,
    universe      = universe_ids_clean,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENSEMBL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
}

go_down <- lapply(risultati_down, function(r) {
  if (length(r$gene_disregolati) < 3) return(NULL)
  do_GO(r$gene_disregolati)   # usa universe_tutti di default
})

go_up <- lapply(risultati_up, function(r) {
  if (length(r$gene_disregolati) < 3) return(NULL)
  do_GO(r$gene_disregolati)
})

# controllo rapido
sapply(go_down, function(x) if (is.null(x)) NA else nrow(as.data.frame(x)))
sapply(go_up,   function(x) if (is.null(x)) NA else nrow(as.data.frame(x)))

dotplot(go_down$basal_only, showCategory = 15) +
  ggtitle("GO (BP) - DOWN, basal_only")

dotplot(go_up$almeno_5_su_6, showCategory = 15) +
  ggtitle("GO (BP) - UP, almeno 5/6")

go_all <- list(DOWN = go_down, UP = go_up)

for (direzione in names(go_all)) {
  for (classe in names(go_all[[direzione]])) {

    go_res <- go_all[[direzione]][[classe]]
    if (is.null(go_res) || nrow(as.data.frame(go_res)) == 0) {
      message("Nessun termine significativo per: ", direzione, " - ", classe)
      next
    }

    p <- dotplot(go_res, showCategory = 15) +
      ggtitle(paste0("GO (BP) - ", direzione, ", ", classe))

    file_out <- file.path(paste0("GO_dotplot_", direzione, "_", classe, ".pdf"))
    ggsave(file_out, p, width = 8, height = 7)
    cat("[done] salvato:", file_out, "\n")
  }
}

########################## DTU and DTE #########################################################
DRIMseq_names <- list.files(path=".",pattern="prioritized_DIUs_DMSO_vs_STM",recursive=TRUE)[-1]
DRIMseq_results <- lapply(DRIMseq_names,function(i){
		foe <- read.table(i,sep="\t",header=TRUE)
		foe <- foe$txID
	})
names(DRIMseq_results) <- sapply(strsplit(DRIMseq_names, "/"), `[`, 2)

res_DTU <- Extract_commons_transcripts(DRIMseq_results, luminal = luminal, basal = basal, min_n = 5)
# riepilogo numerico
#sapply(res_DTU[c("almeno_5_su_6","basal_only","luminal_only")], length)

incrocia_DTU_DTE <- function(res_dtu, res_up, res_down) {

  classi <- c("almeno_5_su_6", "basal_only", "luminal_only")

  do.call(rbind, lapply(classi, function(k) {

    dtu_k <- res_dtu[[k]]
    dte_up_k   <- res_up[[k]]
    dte_down_k <- res_down[[k]]
    dte_k <- union(dte_up_k, dte_down_k)   # DTE in quella classe, UP o DOWN indifferentemente

    comuni       <- intersect(dtu_k, dte_k)
    comuni_up    <- intersect(dtu_k, dte_up_k)
    comuni_down  <- intersect(dtu_k, dte_down_k)

    data.frame(
      classe              = k,
      n_DTU               = length(dtu_k),
      n_DTE               = length(dte_k),
      n_DTU_e_DTE         = length(comuni),
      n_DTU_e_DTE_UP      = length(comuni_up),
      n_DTU_e_DTE_DOWN    = length(comuni_down)
    )
  }))
}

riepilogo_incrocio <- incrocia_DTU_DTE(res_DTU, res_up, res_down)
print(riepilogo_incrocio)




