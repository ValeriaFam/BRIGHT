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
library("msigdbr")
library("ggVennDiagram")
library("stringr")

cell_line <- c("MCF7","BT483","T47D","SUM159","MDAMB231","BT549")
palette <- c(MCF7="#CC6677",BT483="#882255",T47D="#AA4499",SUM159="#117733",MDAMB231="#999933",BT549="#44AA99")
luminal <- cell_line[1:3]
basal <- cell_line[4:6]
path <-  "/projects/CGS_shared/vfama/BRIGHT_PROJECT/"

#Read the tables
names_tab <- list.files(path=".",pattern="sites_per_transcript_and_features.csv",recursive=TRUE)
tab <- lapply(names_tab,function(i){
		read.csv(i,header=TRUE)
	})
names(tab) <- gsub(".*/([A-Za-z0-9]+)_.*\\.csv$", "\\1", names_tab)

#Genes and transcripts not tested are NAs -> renamed in not tested
recode_na <- function(x, new = "not_tested") {
  x <- as.character(x)
  x[is.na(x)] <- new
  factor(x, levels = c("UP", "DOWN", "NS", new))
}

tab <- lapply(tab, function(df) {
  df$regulation   <- recode_na(df$regulation)
  df$status_genes <- recode_na(df$status_genes)
  df
})

#Extract up and down regulated transcripts
get_transcripts <- function(df, status) {
  unique(df$chrom[df$regulation == status])
}
up_tr   <- lapply(tab[cell_line], get_transcripts, status = "UP")
down_tr <- lapply(tab[cell_line], get_transcripts, status = "DOWN")

#Extract up and down regulated genes
get_genes_disreg <- function(df) {
  unique(df$gene_id[df$status_genes == "UP" | df$status_genes == "DOWN" ])
}
disreg_genes   <- lapply(tab[cell_line], get_genes_disreg)

### Do scatter plots DGE/DTE
plot_scatter <- function(df, title) {
  ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(title = title, x = "-log10(padj) DTE", y = "-log10(padj) DGE") +
    theme_classic(base_size = 13)
}

for (cl in names(tab)) {
  df <- tab[[cl]] %>%
    filter(!is.na(p.val), !is.na(padj_genes)) %>%     # solo trascritti E geni testati
    distinct(chrom, gene_id, p.val, padj_genes) %>%   # evita righe duplicate da join
    mutate(
      x = -log10(pmax(p.val, 1e-300)),
      y = -log10(pmax(padj_genes, 1e-300))
    )

  p <- plot_scatter(df, paste0(cl, " - DGE vs DTE (trascritti con m6A)"))
  ggsave(paste0(path, cl, "_DGE_vs_DTE.pdf"), p, width = 6, height = 5)
}

## Barplot on DTE in DGE and viceversa
all_df <- bind_rows(tab, .id = "cell_line")

# 2. DTE: quanti sono anche DGE? (unità = trascritti) ----------------------
dte_df <- all_df %>%
  filter(regulation %in% c("UP", "DOWN"), status_genes != "not_tested") %>%
  mutate(DGE = ifelse(status_genes %in% c("UP", "DOWN"),
                      "Also DGE", "Not DGE")) %>%
  dplyr::count(cell_line, DGE) %>%
  group_by(cell_line) %>%
  mutate(perc = 100 * n / sum(n))

dte_df <- dte_df %>%
  mutate(DGE = factor(DGE, levels = c("Not DGE", "Also DGE")))

p1 <- ggplot(dte_df, aes(x = cell_line, y = perc, fill = DGE)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", perc)),
            position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = c("Also DGE" = "#984EA3", "Not DGE" = "grey70")) +
  labs(x = NULL, y = "% of DTEs", fill = NULL,
       title = "DTEs that are also DGEs") +
  theme_classic(base_size = 13)

ggsave(paste0(path, "_DTE_with_DGE_perc.pdf"), p1, width = 6, height = 5)

# 3. DGE: quanti sono anche DTE? (unità = geni) ---------------------------
# un gene è "DTE" se almeno un suo trascritto è DTE
dge_df <- all_df %>%
  filter(status_genes %in% c("UP", "DOWN")) %>%               # solo DGE
  group_by(cell_line, gene_id) %>%
  summarise(tested  = any(regulation != "not_tested"),        # almeno un trascritto testato
            has_DTE = any(regulation %in% c("UP", "DOWN")),   # almeno un trascritto DTE
            .groups = "drop") %>%
  filter(tested) %>%                                          # escludo i geni mai testati a livello di trascritto
  mutate(DTE = factor(ifelse(has_DTE, "Also DTE", "Not DTE"),
                      levels = c("Not DTE", "Also DTE"))) %>% # primo livello = in cima
  dplyr::count(cell_line, DTE) %>%
  group_by(cell_line) %>%
  mutate(perc = 100 * n / sum(n))

dge_df <- dge_df %>%
  mutate(DTE = factor(DTE, levels = c("Not DTE", "Also DTE")))

p2 <- ggplot(dge_df, aes(x = cell_line, y = perc, fill = DTE)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", perc)),
            position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = c("Also DTE" = "#984EA3", "Not DTE" = "grey70")) +
  labs(x = NULL, y = "% of DGEs", fill = NULL,
       title = "DGEs that are also DTEs") +
  theme_classic(base_size = 13)
  ggsave(paste0(path, "_DGE_with_DTE_perc.pdf"), p2, width = 6, height = 5)
############### How many disregulated genes for luminal/basal?
###Gene Ontology
GO_on_subtype <- function(subtype)
{
	m_subtype <- sapply(disreg_genes[subtype], function(s) {
  		all_genes_subtype <- unique(unlist(disreg_genes[subtype]))
  		all_genes_subtype %in% s
	})
	rownames(m_subtype) <- unique(unlist(disreg_genes[subtype]))
	n_subtype <- rowSums(m_subtype)
	subtype_disreg_genes <- rownames(m_subtype)[n_subtype >= 2]
	
	#Universo
	all_genes_tab <- lapply(tab[subtype], function(df) unique(df$gene_id))
	all_genes_subtype <- unique(unlist(all_genes_tab))
	n_subtype <- sapply(all_genes_subtype, function(g) sum(sapply(all_genes_tab, function(s) g %in% s)))
	universe_subtype <- all_genes_subtype[n_subtype >= 2]
	#Clean versions
	gene_ids_clean     <- sub("\\..*", "", subtype_disreg_genes)
	universe_ids_clean <- sub("\\..*", "", universe_subtype)
	#Gene Ontology
	go_luminal <- enrichGO(
	  gene = gene_ids_clean, universe = universe_ids_clean,
	  OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "ALL",
	  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
	)
}

GO_lum <- GO_on_subtype(luminal) #no term enriched
GO_bas <- GO_on_subtype(basal)   #no term enriched

###GSEA
rank_line <- function(line) {
  df <- tab[[line]]
  df <- df[!is.na(df$padj_genes) & !is.na(df$log2FoldChange_genes), ]
  df$gene <- sub("\\..*", "", df$gene_id)
  p <- pmax(df$padj_genes, 1e-300)
  df$score <- sign(df$log2FoldChange_genes) * -log10(p)
  df <- df[order(-abs(df$score)), ]
  df <- df[!duplicated(df$gene), ]
  r <- setNames(df$score, df$gene)
  r <- r[is.finite(r)]
  sort(r, decreasing = TRUE)
}

t2g <- function(collection, subcollection = NULL) {
  msigdbr(species = "Homo sapiens", collection = collection,
          subcollection = subcollection)[, c("gs_name", "ensembl_gene")]
}
hallmark <- t2g("H")
reactome <- t2g("C2", "CP:REACTOME")

run_gsea_line <- function(line, term2gene) {
  GSEA(geneList = rank_line(line),
       TERM2GENE = term2gene,
       minGSSize = 10, maxGSSize = 500,
       pvalueCutoff = 1, pAdjustMethod = "BH",
       eps = 0, seed = TRUE, verbose = FALSE)
}

lines_all <- c(luminal, basal)
gsea_H <- setNames(lapply(lines_all, run_gsea_line, term2gene = hallmark), lines_all)
gsea_R <- setNames(lapply(lines_all, run_gsea_line, term2gene = reactome), lines_all)

sig_ids <- function(gsea_list, thr = 0.05) {
  lapply(gsea_list, function(g) {
    r <- as.data.frame(g)
    r$ID[r$p.adjust < thr]
  })
}
sigH <- sig_ids(gsea_H)
sigR <- sig_ids(gsea_R)
#sapply(sigH, length); sapply(sigR, length)

comuni_subtype <- function(sig_list, lines, min_n = length(lines)) {
  tb <- table(unlist(sig_list[lines]))
  names(tb)[tb >= min_n]
}

lumH <- comuni_subtype(sigH, luminal)
basH <- comuni_subtype(sigH, basal)
lumR <- comuni_subtype(sigR, luminal)
basR <- comuni_subtype(sigR, basal)

classify <- function(df, lum_ids, bas_ids, lum_lines, bas_lines, max_other = 0) {
  concordant <- function(ids, lines) {
    df %>%
      dplyr::filter(ID %in% ids, line %in% lines) %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(
        ok  = dplyr::n_distinct(sign(NES)) == 1,
        dir = sign(dplyr::first(NES)),
        .groups = "drop"
      ) %>%
      dplyr::filter(ok) %>%
      dplyr::select(ID, dir)
  }
  # numero di linee significative per termine in un dato sottotipo
  n_sig <- function(lines) {
    df %>%
      dplyr::filter(line %in% lines) %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(n = sum(sig, na.rm = TRUE), .groups = "drop")
  }
  nsig_lum <- n_sig(lum_lines) %>% dplyr::rename(n_lum = n)
  nsig_bas <- n_sig(bas_lines) %>% dplyr::rename(n_bas = n)

  lum <- concordant(lum_ids, lum_lines) %>% dplyr::rename(dir_lum = dir)
  bas <- concordant(bas_ids, bas_lines) %>% dplyr::rename(dir_bas = dir)

  dplyr::full_join(lum, bas, by = "ID") %>%
    dplyr::left_join(nsig_lum, by = "ID") %>%
    dplyr::left_join(nsig_bas, by = "ID") %>%
    dplyr::mutate(
      n_lum = tidyr::replace_na(n_lum, 0),
      n_bas = tidyr::replace_na(n_bas, 0),
      class = dplyr::case_when(
        !is.na(dir_lum) & !is.na(dir_bas) & dir_lum == dir_bas ~ "Shared",
        !is.na(dir_lum) & !is.na(dir_bas)                       ~ "Opposite",
        !is.na(dir_lum) & n_bas <= max_other                    ~ "Luminal-specific",
        !is.na(dir_bas) & n_lum <= max_other                    ~ "Basal-specific",
        TRUE                                                    ~ "Mixed"   # significativo in un sottotipo e parzialmente nell'altro
      )
    )
}

collect_gsea <- function(gsea_list) {
  bind_rows(lapply(names(gsea_list), function(l) {
    as.data.frame(gsea_list[[l]]) %>%
      mutate(line = l)
  })) %>%
    dplyr::select(line, ID, NES, p.adjust) %>%
    mutate(sig   = p.adjust < 0.05,
           group = ifelse(line %in% luminal, "Luminal", "Basal"))
}

longH <- collect_gsea(gsea_H)
longR <- collect_gsea(gsea_R)

clsH <- classify(longH, lumH, basH, luminal, basal, max_other = 1)
clsR <- classify(longR, lumR, basR, luminal, basal, max_other = 1)
table(clsH$class); table(clsR$class)

plot_heat_cls <- function(df, cls, prefix, title, wrap = 38) {
  d <- df %>%
    dplyr::inner_join(cls %>% dplyr::select(ID, class), by = "ID") %>%
    dplyr::mutate(
      ID = sub(prefix, "", ID) %>%
        str_replace_all("_", " ") %>%
        str_to_sentence() %>%
        str_wrap(width = wrap),
      class = factor(class, levels = c("Shared", "Luminal-specific",
                                       "Basal-specific", "Opposite", "Other")),
      group = factor(group, levels = c("Basal", "Luminal"))
    )
  ord <- d %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(m = mean(NES), .groups = "drop") %>%
    dplyr::arrange(m)
  d$ID <- factor(d$ID, levels = ord$ID)

  ggplot(d, aes(line, ID, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(data = dplyr::filter(d, sig), label = "*", size = 5, vjust = 0.75) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "NES") +
    facet_grid(class ~ group, scales = "free", space = "free") +
    labs(x = NULL, y = NULL, title = title, caption = "* padj < 0.05") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 8, lineheight = 0.85),
          panel.grid = element_blank(),
          panel.spacing = unit(0.4, "lines"),
          strip.text.y = element_text(angle = 0, face = "bold"),
          strip.text.x = element_text(face = "bold"))
}

save_heat <- function(p, file_base, n_terms, n_lines, label_width = 4.5, outdir = ".") {
  h <- max(5, 2 + 0.25 * n_terms)
  w <- label_width + 2.5 + 0.6 * n_lines   # etichette + strip/legenda + celle
  ggsave(file.path(outdir, paste0(file_base, ".pdf")), p,
         width = w, height = h, limitsize = FALSE)
  ggsave(file.path(outdir, paste0(file_base, ".png")), p,
         width = w, height = h, dpi = 300, limitsize = FALSE)
}

outdir <- "."   

p_H <- plot_heat_cls(longH, clsH, "^HALLMARK_", "Hallmark")
p_R <- plot_heat_cls(longR, clsR, "^REACTOME_", "Reactome")

save_heat(p_H, "GSEA_heatmap_Hallmark",
          n_terms = dplyr::n_distinct(clsH$ID), n_lines = length(lines_all), outdir = outdir)
save_heat(p_R, "GSEA_heatmap_Reactome",
          n_terms = dplyr::n_distinct(clsR$ID),
          n_lines = length(lines_all))

################## Cross dysregulated transcripts with dysregulated genes

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

venn_lum_strict <- function(res, lines, other, titolo, colore) {
  m <- res$matrice
  n_this  <- rowSums(m[, lines, drop = FALSE])
  n_other <- rowSums(m[, other, drop = FALSE])
  ids <- rownames(m)[n_this >= 1 & n_other == 0]

  mm <- m[ids, lines, drop = FALSE]
  s  <- setNames(lapply(lines, function(l) rownames(mm)[mm[, l]]), lines)

  ggVennDiagram(s, label_alpha = 0) +
    scale_fill_gradient(low = "white", high = colore) +
    labs(title = titolo,
         subtitle = paste0(sum(n_this >= 2 & n_other == 0), " in ≥2 linee (recap) + ",
                           sum(n_this == 1 & n_other == 0), " in una sola linea")) +
    theme(legend.position = "none")
}

plots <- list(
  UP_luminal_specific   = venn_lum_strict(res_up,   luminal, basal,   "UP, specifici luminal", "tomato"),
  UP_basal_specific     = venn_lum_strict(res_up,   basal,   luminal, "UP, specifici basal",   "tomato"),
  DOWN_luminal_specific = venn_lum_strict(res_down, luminal, basal,   "DOWN, specifici luminal", "steelblue"),
  DOWN_basal_specific   = venn_lum_strict(res_down, basal,   luminal, "DOWN, specifici basal",   "steelblue")
) 

for (n in names(plots)) {
  ggsave(file.path(paste0("venn_", n, ".pdf")), plots[[n]], width = 5, height = 5)
}
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

    tr_di_geni_disreg <- ids[ids %in% gene_map$chrom[gene_map$gene_id %in% gene_disregolati]]
	tr_non_disreg     <- setdiff(ids, tr_di_geni_disreg)

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
      gene_totali       = gene_totali,
      tr_di_geni_disreg =tr_di_geni_disreg,
      tr_non_disreg = tr_non_disreg
    )
  })
  names(risultati) <- names(classi)
  risultati
}

cell_lines <- dimnames(res_up$matrice)[[2]]

risultati_down <- conta_gene_disregolati(res_down, tab, cell_lines)
risultati_up   <- conta_gene_disregolati(res_up,   tab, cell_lines)

do.call(rbind, lapply(risultati_up, `[[`, "summary"))   # riepilogo numerico

universe_tutti <- unique(unlist(lapply(tab, function(df) df$gene_id)))
cat(length(universe_tutti), "geni nell'universe\n")

# do_GO <- function(gene_ids, universe_ids = universe_tutti, ont = "BP") {
#   gene_ids_clean     <- sub("\\..*", "", gene_ids)
#   universe_ids_clean <- sub("\\..*", "", universe_ids)

#   enrichGO(
#     gene          = gene_ids_clean,
#     universe      = universe_ids_clean,
#     OrgDb         = org.Hs.eg.db,
#     keyType       = "ENSEMBL",
#     ont           = ont,
#     pAdjustMethod = "BH",
#     pvalueCutoff  = 0.05,
#     qvalueCutoff  = 0.2,
#     readable      = TRUE
#   )
# }

# go_down <- lapply(risultati_down, function(r) {
#   if (length(r$gene_disregolati) < 3) return(NULL)
#   do_GO(r$gene_disregolati)   # usa universe_tutti di default
# })

# go_up <- lapply(risultati_up, function(r) {
#   if (length(r$gene_disregolati) < 3) return(NULL)
#   do_GO(r$gene_disregolati)
# })

# # controllo rapido
# sapply(go_down, function(x) if (is.null(x)) NA else nrow(as.data.frame(x)))
# sapply(go_up,   function(x) if (is.null(x)) NA else nrow(as.data.frame(x)))

# dotplot(go_down$basal_only, showCategory = 15) +
#   ggtitle("GO (BP) - DOWN, basal_only")

# dotplot(go_up$almeno_5_su_6, showCategory = 15) +
#   ggtitle("GO (BP) - UP, almeno 5/6")

# go_all <- list(DOWN = go_down, UP = go_up)

# for (direzione in names(go_all)) {
#   for (classe in names(go_all[[direzione]])) {

#     go_res <- go_all[[direzione]][[classe]]
#     if (is.null(go_res) || nrow(as.data.frame(go_res)) == 0) {
#       message("Nessun termine significativo per: ", direzione, " - ", classe)
#       next
#     }

#     p <- dotplot(go_res, showCategory = 15) +
#       ggtitle(paste0("GO (BP) - ", direzione, ", ", classe))

#     file_out <- file.path(paste0("GO_dotplot_", direzione, "_", classe, ".pdf"))
#     ggsave(file_out, p, width = 8, height = 7)
#     cat("[done] salvato:", file_out, "\n")
#   }
# }

#Compare the difference in stoichiometry between the groups
media_siti <- function(x) sapply(strsplit(x, ","), function(v) mean(as.numeric(v), na.rm = TRUE))

tabella_23_33 <- function(res, tab, cell_lines, lines, classe, direzione,
                          col = "m6A_diff_dmso_storm", verbose = TRUE) {
  m <- res$matrice
  r <- conta_gene_disregolati(res, tab, cell_lines)[[classe]]

  ids <- c(r$tr_di_geni_disreg, r$tr_non_disreg)
  n   <- rowSums(m[ids, lines, drop = FALSE])
  ids <- ids[n >= 2]                                            # 2/3 e 3/3

  # ---- tabella principale (invariata) ----
  long <- do.call(rbind, lapply(lines, function(cl) {
    df <- tab[[cl]]
    df <- df[df$chrom %in% ids, c("chrom", col)]
    df <- df[!duplicated(df$chrom), ]
    data.frame(transcript = df$chrom, line = cl, value = media_siti(df[[col]]))
  }))

  long$dte    <- m[cbind(long$transcript, long$line)]
  long$classe <- ifelse(n[long$transcript] == length(lines), "3/3", "2/3")
  long$gruppo <- ifelse(long$transcript %in% r$tr_di_geni_disreg,
                        "gene disregolato", "gene non disregolato")

  d <- aggregate(value ~ transcript + dte + classe + gruppo, data = long, FUN = mean)
  d$categoria <- ifelse(d$dte, paste0("DTE ", d$classe), "linea non DTE (2/3)")
  d$categoria <- factor(d$categoria,
                        levels = c("DTE 3/3", "DTE 2/3", "linea non DTE (2/3)"))
  d$classe_trascritto <- classe
  d$direzione         <- direzione

  # ---- CONTROLLI sui trascritti DTE in esattamente 2 linee ----
  ids2 <- ids[rowSums(m[ids, lines, drop = FALSE]) == 2]

  # una riga per sito (trascritto, posizione, valore), per linea
  parse_sites <- function(df) {
    s <- df[!duplicated(df$chrom), c("chrom", "m6A_start_positions", col)]
    s <- s[!is.na(s$m6A_start_positions) & !is.na(s[[col]]), ]
    out <- lapply(seq_len(nrow(s)), function(i) {
      pos <- strsplit(s$m6A_start_positions[i], ",")[[1]]
      val <- as.numeric(strsplit(s[[col]][i], ",")[[1]])
      if (length(pos) != length(val)) return(NULL)
      data.frame(transcript = s$chrom[i], pos = as.integer(pos), value = val)
    })
    do.call(rbind, out)
  }
  sites <- lapply(setNames(lines, lines), function(cl)
    parse_sites(tab[[cl]][tab[[cl]]$chrom %in% ids2, ]))

  get_pos <- function(cl, t) sites[[cl]]$pos[sites[[cl]]$transcript == t]
  get_val <- function(cl, t) sites[[cl]]$value[sites[[cl]]$transcript == t]

  # A: sovrapposizione delle posizioni tra le 2 linee DTE
  ov <- do.call(rbind, lapply(ids2, function(t) {
    dte <- lines[m[t, lines]]
    a <- get_pos(dte[1], t); b <- get_pos(dte[2], t)
    data.frame(transcript = t, linea1 = dte[1], linea2 = dte[2],
               n_linea1 = length(a), n_linea2 = length(b),
               n_condivisi = length(intersect(a, b)),
               jaccard = if (length(union(a, b)) == 0) NA
                         else length(intersect(a, b)) / length(union(a, b)))
  }))

  # B: concordanza delle stechiometrie (media sui siti DENTRO la linea)
  cc <- do.call(rbind, lapply(ids2, function(t) {
    dte <- lines[m[t, lines]]
    data.frame(transcript = t, linea1 = dte[1], linea2 = dte[2],
               val1 = mean(get_val(dte[1], t)),
               val2 = mean(get_val(dte[2], t)))
  }))
  cc <- cc[complete.cases(cc), ]

  rho <- if (nrow(cc) > 2) suppressWarnings(
    cor.test(cc$val1, cc$val2, method = "spearman")) else NULL

  controlli <- data.frame(
    classe = classe, direzione = direzione,
    n_trascritti_2su3      = length(ids2),
    jaccard_mediano        = median(ov$jaccard, na.rm = TRUE),
    quota_senza_siti_comuni = mean(ov$n_condivisi == 0),
    n_confrontabili        = nrow(cc),
    spearman_rho           = if (!is.null(rho)) unname(rho$estimate) else NA,
    spearman_p             = if (!is.null(rho)) rho$p.value else NA,
    quota_stessa_direzione = if (nrow(cc)) mean(sign(cc$val1) == sign(cc$val2)) else NA
  )

  if (verbose) print(controlli, row.names = FALSE)

  attr(d, "controlli")       <- controlli
  attr(d, "overlap_siti")    <- ov
  attr(d, "concordanza_val") <- cc
  d
}


# le 4 combinazioni
d_all <- rbind(
  tabella_23_33(res_up,   tab, cell_lines, luminal, "luminal_only", "UP"),
  tabella_23_33(res_up,   tab, cell_lines, basal,   "basal_only",   "UP"),
  tabella_23_33(res_down, tab, cell_lines, luminal, "luminal_only", "DOWN"),
  tabella_23_33(res_down, tab, cell_lines, basal,   "basal_only",   "DOWN")
)

# numerosità per categoria: guardala prima dei p-value
table(d_all$direzione, d_all$classe_trascritto, d_all$gruppo, d_all$categoria)

cols <- c("DTE 3/3"             = "#C0392B",
          "DTE 2/3"             = "#E67E22",
          "linea non DTE" = "#7F8C8D")

tema_m6A <- theme_bw(base_size = 12) +
  theme(
    panel.grid       = element_blank(),                       # niente griglie
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.spacing    = unit(1.2, "lines"),                    # spazio tra i riquadri
    strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.6),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 25, hjust = 1),
    axis.ticks.length = unit(3, "pt"),
    legend.position  = "none",
    plot.title       = element_text(face = "bold", hjust = 0.5)
  )

for (direzione in c("UP", "DOWN")) {
  p <- ggplot(d_all[d_all$direzione == direzione, ],
              aes(categoria, value, fill = categoria)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, linewidth = 0.5) +
    geom_jitter(width = 0.12, size = 1.2, alpha = 0.6, colour = "black") +
    scale_fill_manual(values = cols) +
    scale_x_discrete(labels = function(x) gsub(" \\(", "\n(", x)) +   # va a capo sull'etichetta lunga
    facet_grid(classe_trascritto ~ gruppo) +
    labs(x = NULL, y = "Riduzione media m6A per trascritto",
         title = paste("DTE", direzione)) +
    tema_m6A

  ggsave(paste0("m6A_reduction_DTE_", direzione, ".pdf"),
         p, width = 9, height = 7)
}

res_test <- lapply(split(d_all, list(d_all$direzione, d_all$classe_trascritto, d_all$gruppo)),
  function(x) {
    a <- x[x$categoria == "DTE 2/3", c("transcript","value")]
    b <- x[x$categoria == "linea non DTE (2/3)", c("transcript","value")]
    k <- merge(a, b, by = "transcript")
    y <- x[x$categoria %in% c("DTE 3/3", "DTE 2/3"), ]
    data.frame(n_appaiati = nrow(k),
      p_2su3_vs_non_DTE = if (nrow(k) > 1) wilcox.test(k$value.x, k$value.y, paired = TRUE)$p.value else NA,
      p_33_vs_23 = if (length(unique(y$categoria)) == 2) wilcox.test(value ~ categoria, data = y)$p.value else NA)
  })
do.call(rbind, res_test)

######
lfc_col <- "log2FoldChange_genes"

lfc_per_tr <- function(risultati, res, tab, cell_lines, lfc_col) {
  dplyr::bind_rows(lapply(names(risultati), function(k) {
    r   <- risultati[[k]]
    ids <- c(r$tr_di_geni_disreg, r$tr_non_disreg)

    long <- dplyr::bind_rows(lapply(cell_lines, function(cl) {
      presenti <- ids[ids %in% rownames(res$matrice)[res$matrice[, cl]]]
      s <- tab[[cl]][tab[[cl]]$chrom %in% presenti, c("chrom", "gene_id", lfc_col)]
      if (nrow(s) == 0) return(NULL)
      names(s)[3] <- "lfc"
      s$cell_line <- cl
      s
    }))

    long %>%
      dplyr::group_by(chrom, gene_id) %>%
      dplyr::summarise(mean_lfc     = mean(lfc, na.rm = TRUE),
                       mean_abs_lfc = mean(abs(lfc), na.rm = TRUE),
                       n_lines      = dplyr::n(),
                       .groups = "drop") %>%
      dplyr::mutate(
        classe = k,
        gruppo = ifelse(chrom %in% r$tr_di_geni_disreg,
                        "Dysregulated genes", "Non dysregulated genes"))
  }))
}

df_lfc_down <- lfc_per_tr(risultati_down, res_down, tab, cell_lines, lfc_col)
df_lfc_up <- lfc_per_tr(risultati_up, res_up, tab, cell_lines, lfc_col)

# p-value Wilcoxon per classe, da annotare sul plot
pv_down <- df_lfc_down %>%
  dplyr::group_by(classe) %>%
  dplyr::summarise(
    p = tryCatch(wilcox.test(mean_abs_lfc ~ gruppo)$p.value, error = function(e) NA_real_),
    y = max(mean_abs_lfc, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::mutate(label = paste0("Wilcoxon p = ", signif(p, 2)))

p <- ggplot(df_lfc_down, aes(gruppo, mean_abs_lfc, fill = gruppo)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_text(data = pv_down, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Dysregulated genes" = "#B2182B", "Non dysregulated genes" = "grey70")) +
  facet_wrap(~ classe, scales = "free_y") +
  labs(x = NULL, y = "Mean log2FC per gene") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("boxplot_lfc_down_disreg_vs_non.pdf", p, width = 9, height = 4.5)

####################DTE per cell line
sites_line <- function(cl, col = "m6A_diff_dmso_storm") {
  df <- tab[[cl]]
  df <- df[!duplicated(df$chrom) & !is.na(df$m6A_start_positions) & !is.na(df[[col]]), ]
  out <- lapply(seq_len(nrow(df)), function(i) {
    pos <- strsplit(df$m6A_start_positions[i], ",")[[1]]
    val <- as.numeric(strsplit(df[[col]][i], ",")[[1]])
    if (length(pos) != length(val)) return(NULL)   # salta righe incoerenti
    data.frame(transcript = df$chrom[i], gene_id = df$gene_id[i],
               regulation = df$regulation[i],
               pos = as.integer(pos), diff_storm = val)
  })
  out <- do.call(rbind, out)
  out$cell_line <- cl
  out
}

sites_all <- do.call(rbind, lapply(names(tab), sites_line))
sites_all$regulation <- factor(sites_all$regulation, levels = c("NS", "UP", "DOWN"))

# controlli
table(sites_all$cell_line, sites_all$regulation) 
sites_all$cell_line <- factor(sites_all$cell_line, levels = names(tab))

ordine <- names(tab)
i <- match(c("T47D", "BT549"), ordine)
ordine[i] <- ordine[rev(i)]          # scambia le due posizioni
ordine                               # controlla il nuovo ordine

sites_all$cell_line <- factor(sites_all$cell_line, levels = ordine)

p <- ggplot(sites_all, aes(regulation, diff_storm, fill = cell_line)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.55,
               linewidth = 0.5, colour = "black") +
  scale_fill_manual(values = palette) +
  scale_colour_manual(values = palette) +
  facet_wrap(~ cell_line, nrow = 2) +
  labs(x = "Regulation del trascritto",
       y = "Differenza stechiometria m6A (DMSO - STORM)") +
  theme_bw(base_size = 12) +
  theme(panel.grid       = element_blank(),
        panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        panel.spacing    = unit(1.2, "lines"),
        strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.6),
        strip.text       = element_text(face = "bold"),
        legend.position  = "none")

ggsave("m6A_diff_by_regulation_DTE_per_line.pdf", p, width = 10, height = 7)

####################################### DTU and DTE ###########################################
DRIMseq_names <- list.files(path=".",pattern="prioritized_DIUs_DMSO_vs_STM",recursive=TRUE)[-1]
DRIMseq_results <- lapply(DRIMseq_names,function(i){
		foe <- read.table(i,sep="\t",header=TRUE)
		foe <- foe$txID
	})
names(DRIMseq_results) <- sapply(strsplit(DRIMseq_names, "/"), `[`, 2)

Scisorseq_results <- readRDS("Scisorseq_ENST.rds")
Scisorseq_iso <- lapply(Scisorseq_results,function(i)unique(i$Isoform))

res_DTU_drimseq <- Extract_commons_transcripts(DRIMseq_results, luminal = luminal, basal = basal, min_n = 5)
res_DTU_scisorseq <- Extract_commons_transcripts(Scisorseq_iso, luminal = luminal, basal = basal, min_n = 5)
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

riepilogo_incrocio_drimseq <- incrocia_DTU_DTE(res_DTU_drimseq, res_up, res_down)
riepilogo_incrocio_scisorseq <- incrocia_DTU_DTE(res_DTU_scisorseq, res_up, res_down)



