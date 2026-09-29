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
up   <- lapply(tab[cell_line], get_transcripts, status = "UP")
down <- lapply(tab[cell_line], get_transcripts, status = "DOWN")

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

res_up   <- Extract_commons_transcripts(up,   luminal, basal)
res_down <- Extract_commons_transcripts(down, luminal, basal)

#Summary
recap <- data.frame(
  categoria = c("UP", "DOWN"),
  n_almeno_5_su_6 = c(length(res_up$almeno_5_su_6), length(res_down$almeno_5_su_6)),
  basal_tutte     = c(length(res_up$basal_tutte),   length(res_down$basal_tutte)),
  basal_only      = c(length(res_up$basal_only),    length(res_down$basal_only)),
  luminal_tutte   = c(length(res_up$luminal_tutte), length(res_down$luminal_tutte)),
  luminal_only    = c(length(res_up$luminal_only),  length(res_down$luminal_only))
)
print(recap)




