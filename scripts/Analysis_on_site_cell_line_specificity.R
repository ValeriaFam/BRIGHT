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

cell_line <- c("MCF7","BT483","T47D","SUM159","MDAMB231","BT549")
palette <- c("#CC6677","#882255","#AA4499","#117733","#999933","#44AA99")

#Read the tables
names_tab <- list.files(path=".",pattern="sites_per_transcript_and_features.csv",recursive=TRUE)
tab <- lapply(names_tab,function(i){
		read.csv(i,header=TRUE)
	})
names(tab) <- gsub(".*/([A-Za-z0-9]+)_.*\\.csv$", "\\1", names_tab)

#Split each row in as many rows as the sites
tab_one_site_per_row <- lapply(tab,function(i){
	i %>% separate_rows(m6A_start_positions, m6A_stoich_dmso_mean, m6A_diff_dmso_storm, sep = ",") %>%
	  mutate(
	    m6A_start_positions = as.numeric(m6A_start_positions),
	    m6A_stoich_dmso_mean = as.numeric(m6A_stoich_dmso_mean),
	    m6A_diff_dmso_storm = as.numeric(m6A_diff_dmso_storm)
	  )
 })

#Extract the keys
keys_list <- imap(tab_one_site_per_row, function(df, nome) {
  df %>%
    distinct(chrom, m6A_start_positions) %>%
    mutate(cell_line = nome)
})

all_keys <- bind_rows(keys_list)

#Table wide
consensus_table <- all_keys %>%
  distinct(chrom, m6A_start_positions, cell_line) %>%
  mutate(presente = 1) %>%
  tidyr::pivot_wider(
    names_from = cell_line,
    values_from = presente,
    values_fill = 0
  )
consensus_table <- consensus_table %>% dplyr::select(chrom, m6A_start_positions, BT483, MCF7, T47D, MDAMB231, SUM159, BT549)

#Upset plot - site shared across cell lines
cell_lines <- colnames(consensus_table)[-(1:2)]
p <- upset(
  consensus_table,
  cell_lines,
  sort_sets = FALSE,
  name = "Cell line",
  width_ratio = 0.15,
  min_size = 1, base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 3, color = "black") 
    )
  )
)
ggsave("Upset_plot_m6A_sites_all_cell_lines.pdf", p, width = 30, height = 8)

#Upset plot - division by subtype
luminal_cols <- c("BT483", "MCF7", "T47D")
basal_cols   <- c("SUM159", "MDAMB231", "BT549")

consensus_table <- consensus_table %>%
  mutate(
    n_luminal = rowSums(across(all_of(luminal_cols))),
    n_basal   = rowSums(across(all_of(basal_cols))),
    subtype = case_when(
      n_luminal >= 2 & n_basal == 0            ~ "Luminal",
      n_luminal == 0 & n_basal >= 2            ~ "Basal",
      n_luminal >= 2 & n_basal >= 2            ~ "Pancancer",
      TRUE                                     ~ "Others"
    )
  ) %>%
  dplyr::select(-n_luminal, -n_basal)

consensus_table <- consensus_table %>%
  mutate(
    Luminal   = as.integer(subtype == "Luminal"),
    Basal     = as.integer(subtype == "Basal"),
    Pancancer = as.integer(subtype == "Pancancer"),
    Others    = as.integer(subtype == "Others")
  )

subtype_cols <- c("Luminal", "Basal", "Pancancer", "Others")

p <- upset(
  consensus_table,
  subtype_cols,
  name = "Subtype",
  width_ratio = 0.15,
  min_size = 1
)

ggsave("Upset_plot_m6A_subtype.pdf", p, width = 10, height = 8)



