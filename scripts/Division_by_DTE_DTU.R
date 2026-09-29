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

res_up   <- Extract_commons_transcripts(up,   luminal, basal)
res_down <- Extract_commons_transcripts(down, luminal, basal)

#Summary
recap <- data.frame(
  categoria = c("UP", "DOWN"),
  n_almeno_5_su_6 = c(length(res_up$almeno_5_su_6), length(res_down$almeno_5_su_6)),
  basal_only      = c(length(res_up$basal_only),    length(res_down$basal_only)),
  luminal_only    = c(length(res_up$luminal_only),  length(res_down$luminal_only))
)
print(recap)

#How many DTE are also DGE? (quanti trascritti nella mia classe appartengono a un gene disregolato)
conta_concordanza_gene <- function(res, tab, cell_lines, min_5_6 = 5, min_2_3 = 2) {

  classi <- list(
    almeno_5_su_6 = list(ids = res$almeno_5_su_6, min_n = min_5_6),
    basal_only    = list(ids = res$basal_only,    min_n = min_2_3),
    luminal_only  = list(ids = res$luminal_only,  min_n = min_2_3)
  )

  do.call(rbind, lapply(names(classi), function(k) {
    ids   <- classi[[k]]$ids
    min_n <- classi[[k]]$min_n

    # per ogni trascritto della classe, prendi lo status_genes riportato
    # nelle cell line in cui il trascritto e' presente (matrice TRUE)
    long <- do.call(rbind, lapply(cell_lines, function(cl) {
      presenti <- ids[ids %in% rownames(res$matrice)[res$matrice[, cl]]]
      if (length(presenti) == 0) return(NULL)

      df <- tab[[cl]]
      sub <- df[df$chrom %in% presenti, c("chrom", "status_genes")]
      if (nrow(sub) == 0) return(NULL)
      sub$cell_line <- cl
      sub
    }))

    # n. di occorrenze in cui il gene e' UP/DOWN, per trascritto
    n_disreg <- tapply(long$status_genes %in% c("UP", "DOWN"), long$chrom, sum)

    concordanti <- names(n_disreg)[n_disreg >= min_n]

    data.frame(
      classe                   = k,
      soglia_min               = min_n,
      n_transcript_totali      = length(ids),
      n_transcript_concordanti = length(concordanti),
      perc_concordanti         = round(100 * length(concordanti) / length(ids), 1)
    )
  }))
}

cell_lines <- dimnames(res_up$matrice)[[2]]

riepilogo_up   <- conta_concordanza_gene(res_up,   tab, cell_lines)
riepilogo_down <- conta_concordanza_gene(res_down, tab, cell_lines)

print(riepilogo_up)
print(riepilogo_down)

########################################################################################

costruisci_df <- function(res, tab) {

  cell_lines <- dimnames(res$matrice)[[2]]   # nomi delle 6 linee, dalla matrice

  classi <- list(
    almeno_5_su_6 = res$almeno_5_su_6,
    basal_only    = res$basal_only,
    luminal_only  = res$luminal_only
  )

  # formato lungo: trascritto x cell line x classe, solo dove il trascritto e' nella classe e presente in quella cell line (matrice TRUE)
  long <- do.call(rbind, lapply(names(classi), function(k) {
    ids <- classi[[k]]
    do.call(rbind, lapply(cell_lines, function(cl) {
      # trascritti di questa classe presenti in questa cell line (matrice TRUE)
      presenti <- ids[ids %in% rownames(res$matrice)[res$matrice[, cl]]]
      if (length(presenti) == 0) return(NULL)

      df <- tab[[cl]]
      sub <- df[df$chrom %in% presenti, c("chrom", "log2FC", "p.val")]
      if (nrow(sub) == 0) return(NULL)
      sub$classe    <- k
      sub$cell_line <- cl
      sub
    }))
  }))

  # formato largo: una riga per trascritto/classe, colonne log2FC.<cl> e p.val.<cl>
  wide <- reshape(long, idvar = c("chrom", "classe"),
                  timevar = "cell_line", direction = "wide")
  rownames(wide) <- NULL

  lfc_cols <- grep("^log2FC\\.", names(wide), value = TRUE)
  p_cols   <- grep("^p\\.val\\.", names(wide), value = TRUE)

  # mediana solo per almeno_5_su_6
  is5 <- wide$classe == "almeno_5_su_6"
  wide$log2FC_mediana <- NA_real_
  wide$p.val_mediana  <- NA_real_
  wide$log2FC_mediana[is5] <- apply(wide[is5, lfc_cols, drop = FALSE], 1, median, na.rm = TRUE)
  wide$p.val_mediana[is5]  <- apply(wide[is5, p_cols,   drop = FALSE], 1, median, na.rm = TRUE)

  wide$n_cell_lines <- rowSums(!is.na(wide[, lfc_cols, drop = FALSE]))

  wide
}

df_down <- costruisci_df(res_down, tab)
df_up   <- costruisci_df(res_up,   tab)  

table(df_down$classe)
head(df_down)

prepara_plot_df <- function(wide, direzione) {

  lfc_cols <- grep("^log2FC\\.", names(wide), value = TRUE)
  p_cols   <- grep("^p\\.val\\.", names(wide), value = TRUE)

  # --- almeno_5_su_6: un punto per trascritto, valore = mediana ---
  is5 <- wide$classe == "almeno_5_su_6"
  agg <- data.frame(
    chrom      = wide$chrom[is5],
    classe     = "almeno_5_su_6",
    log2FC     = apply(wide[is5, lfc_cols, drop = FALSE], 1, median, na.rm = TRUE),
    p.val      = apply(wide[is5, p_cols,   drop = FALSE], 1, median, na.rm = TRUE)
  )

  # --- basal_only / luminal_only: un punto per trascritto x cell line ---
  altre <- wide[!is5, ]
  long <- do.call(rbind, lapply(seq_len(nrow(altre)), function(i) {
    row <- altre[i, ]
    do.call(rbind, lapply(lfc_cols, function(lc) {
      cl <- sub("^log2FC\\.", "", lc)
      pc <- paste0("p.val.", cl)
      if (is.na(row[[lc]])) return(NULL)
      data.frame(chrom = row$chrom, classe = row$classe,
                log2FC = row[[lc]], p.val = row[[pc]])
    }))
  }))

  out <- rbind(agg, long)
  out$direzione <- direzione
  out
}

plot_df <- rbind(
  prepara_plot_df(df_up,   "UP"),
  prepara_plot_df(df_down, "DOWN")
)

plot_df$neglog10p <- -log10(plot_df$p.val)
plot_df$classe    <- factor(plot_df$classe,
                            levels = c("almeno_5_su_6", "basal_only", "luminal_only"))
plot_df$direzione <- factor(plot_df$direzione, levels = c("UP", "DOWN"))
plot_df <- plot_df[order(plot_df$classe), ]   # disegna prima i più numerosi

cols <- c(almeno_5_su_6 = "grey60", basal_only = "#999933", luminal_only = "#CC6677")

p <- ggplot(plot_df, aes(x = log2FC, y = neglog10p, colour = classe)) +
  geom_point(alpha = 0.6, size = 1.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = cols, name = "Classe",
                      labels = c("≥5/6 (mediana)", "basal only", "luminal only")) +
  facet_wrap(~ direzione, ncol = 2) +
  labs(x = "log2FC", y = "-log10(p.val)") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"))


ggsave("DTE_log2FC_vs_pval.pdf", p, width = 10, height = 8)
