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

#Compute medians
medians <- sapply(tab, function(df) {
  pos <- df$m6A_start_positions
  n <- ifelse(is.na(pos) | pos == "", 0L,
              lengths(strsplit(pos, ",", fixed = TRUE)))
  median(n)
})
#BT483    BT549     MCF7 MDAMB231   SUM159     T47D 
#    4        3        3        4        4        3 

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

#Load file with coverage
coverages_boost <- read.table("coverage_summary_boost.tsv", sep = "\t", stringsAsFactors = FALSE)
to_remove <- c("T47D_DMSO_1", "T47D_DMSO_3", "T47D_STM_3")
coverages_boost <- coverages_boost[
  !grepl(paste0("/(", paste(to_remove, collapse = "|"), ")/"), coverages_boost$V1),
  , drop = FALSE
]
coverages <- read.table("coverage_summary.tsv", sep = "\t", stringsAsFactors = FALSE) 
coverages <- coverages[
  !grepl("MDAMB231|T47D", coverages$V1) |
    grepl("/T47D_DMSO_1/", coverages$V1, fixed = TRUE),
  , drop = FALSE
]

coverages_all <- rbind(coverages,coverages_boost) %>%  mutate(
    path      = sub("^\\./", "", V1),
    cell_line = sapply(strsplit(path, "/"), `[`, 1),
    sample    = sapply(strsplit(path, "/"), `[`, 3),
    coverage  = V2
  )

coverage_mean <- coverages_all %>%
  group_by(cell_line) %>%
  summarise(mean_coverage = mean(coverage, na.rm = TRUE))

#Upset plot - site shared across cell lines
cell_lines <- colnames(consensus_table)[-(1:2)]

size_df <- consensus_table %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(cell_lines), sum)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "cell_line",
    values_to = "size"
  )

order_levels <- size_df %>%
  dplyr::arrange(dplyr::desc(size)) %>%
  dplyr::pull(cell_line)

coverage_df <- coverage_mean %>%
  dplyr::mutate(
    cell_line = dplyr::if_else(
      cell_line == "T47D_boost",
      "T47D",
      cell_line
    )
  ) %>%
  dplyr::filter(cell_line %in% cell_lines) %>%
  dplyr::group_by(cell_line) %>%
  dplyr::summarise(
    mean_coverage = mean(mean_coverage, na.rm = TRUE),
    .groups = "drop"
  )

coverage_df <- coverage_df %>% #è stata scalata
  dplyr::mutate(
    coverage_scaled =
      mean_coverage /
      max(mean_coverage, na.rm = TRUE) *
      max(size_df$size, na.rm = TRUE)
  )

coverage_lookup <- coverage_df$coverage_scaled
names(coverage_lookup) <- coverage_df$cell_line


final_plot <- ComplexUpset::upset(
  consensus_table,
  intersect = cell_lines,
  sort_sets = "descending",
  name = "Cell line",
  width_ratio = 0.25,
  min_size = 1,
  encode_sets = FALSE,

  set_sizes = (

    # SIZE
    ComplexUpset::upset_set_size(
      geom = ggplot2::geom_bar(
        aes(
          x = group,
          y = after_stat(count)
        ),
        width = 0.32,
        fill = "grey40",
        position = ggplot2::position_nudge(
          x = -0.18
  	)))

    +
    ggplot2::geom_col(
      aes(
        x = group,
        y = coverage_lookup[as.character(group)]
      ),
      width = 0.32,
      fill = "grey70",
      position = ggplot2::position_nudge(
        x = 0.18
      ))  +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  ),
  base_annotations = list(
    "Intersection size" =
      ComplexUpset::intersection_size(
        text = list(
          size = 3,
          color = "black"
        ))))

ggsave(
  "Upset_plot_m6A_sites_all_cell_lines.pdf",
  final_plot,
  width = 30,
  height = 8
)


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

#Boxplot for the number of sites/length of the transcript
txdb <- makeTxDbFromGFF("/projects/CGS_shared/vfama/BRIGHT_PROJECT/filtered_corrected_assembly.gtf")
tx_lengths_list <- exonsBy(txdb, by = "tx", use.names = TRUE)
tx_lengths <- data.frame(
  ensembl_transcript_id = names(tx_lengths_list),
  transcript_length     = sum(width(tx_lengths_list))
)

#numeri di trascritti con 1,2,3,4,5 o + siti
lapply(names(tab), function(k){

  df <- tab[[k]] %>%
    mutate(
      n_sites = str_count(m6A_start_positions, ",") + 1,
      site_cat = case_when(
        n_sites == 1 ~ "1",
        n_sites == 2 ~ "2",
        n_sites == 3 ~ "3",
        n_sites == 4 ~ "4",
        n_sites >= 5 ~ "5+"
      ),
      site_cat = factor(site_cat, levels = c("1", "2", "3", "4", "5+"))
    )

  counts <- df %>% dplyr::count(site_cat)

  p <- ggplot(counts, aes(x = site_cat, y = n)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    geom_text(aes(label = n), vjust = -0.4, size = 4) +
    labs(
      x = "Numero di siti m6A per trascritto",
      y = "Numero di trascritti",
      title = paste0("Numero di trascritti per numero di siti m6A - ", k)
    ) +
    theme_bw(base_size = 13)

  ggsave(paste0(k, "_nTranscriptsBySiteCount.pdf"), p, width = 7, height = 6)
})

#Boxplot lunghezza dei trascritti con 1,2,3,4,5 o + siti
lapply(names(tab), function(k){

  df <- tab[[k]] %>%
    mutate(
      n_sites = str_count(m6A_start_positions, ",") + 1,
      site_cat = case_when(
        n_sites == 1 ~ "1",
        n_sites == 2 ~ "2",
        n_sites == 3 ~ "3",
        n_sites == 4 ~ "4",
        n_sites >= 5 ~ "5+"
      ),
      site_cat = factor(site_cat, levels = c("1", "2", "3", "4", "5+"))
    )

  df_plot <- df %>%
    left_join(tx_lengths, by = c("chrom" = "ensembl_transcript_id")) %>%
    filter(!is.na(transcript_length))

  p <- ggplot(df_plot, aes(x = site_cat, y = transcript_length)) +
    geom_boxplot(
	  outlier.size = 0.8,
	  fill = palette[k],
	  alpha = 0.6,
	  varwidth = TRUE
	)+
    labs(
      x = "Numero di siti m6A per trascritto",
      y = "Lunghezza trascritto",
      title = paste0("Lunghezza dei trascritti in funzione del numero di siti m6A - ", k)
    ) +
    theme_bw(base_size = 13)

  ggsave(paste0(k, "_lengthOfTranscriptsAndSites.pdf"), p, width = 8, height = 8)
})

#Boxplot lunghezza trascritti divisa per quartili 
lapply(names(tab), function(k){

  df <- tab[[k]] %>%
    mutate(
      n_sites = str_count(m6A_start_positions, ",") + 1,
      site_cat = case_when(
        n_sites == 1 ~ "1",
        n_sites == 2 ~ "2",
        n_sites == 3 ~ "3",
        n_sites == 4 ~ "4",
        n_sites >= 5 ~ "5+"
      ),
      site_cat = factor(site_cat, levels = c("1", "2", "3", "4", "5+"))
    )

  df_plot <- df %>%
    left_join(tx_lengths, by = c("chrom" = "ensembl_transcript_id")) %>%
    filter(!is.na(transcript_length)) %>%
    group_by(site_cat) %>%
    mutate(length_quartile = factor(
      ntile(transcript_length, 4),
      levels = 1:4,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )) %>%
    ungroup()

  p <- ggplot(df_plot, aes(x = site_cat, y = transcript_length, fill = length_quartile)) +
    geom_boxplot(
	  outlier.size = 0.6,
	  alpha = 0.8,
	  varwidth = TRUE,
	  position = position_dodge(width = 0.8)
	)
    labs(
      x = "Numero di siti m6A per trascritto",
      y = "Lunghezza trascritto (nt, log10)",
      fill = "Quartile lunghezza",
      title = paste0("Lunghezza dei trascritti per quartile e numero di siti m6A - ", k)
    ) +
    theme_bw(base_size = 13)

  ggsave(paste0(k, "_lengthQuartilesBySiteCount.pdf"), p, width = 10, height = 8)
})

#Boxplot differenze stoichiometrie DMSO divisi per numero di siti
lapply(names(tab), function(k) {
  df_plot <- tab[[k]] %>%
    mutate(
      n_sites = str_count(m6A_start_positions, ",") + 1,
      site_cat = case_when(
        n_sites == 2 ~ "2",
        n_sites == 3 ~ "3",
        n_sites == 4 ~ "4",
        n_sites >= 5 ~ "5+"
      ),
      stoich_pct_change = sapply(strsplit(m6A_stoich_dmso_mean, ","), function(x) {
        x <- as.numeric(x)
        min_val <- min(x, na.rm = TRUE)
        max_val <- max(x, na.rm = TRUE)
        if (is.finite(max_val) && max_val > 0) (max_val - min_val) / max_val * 100 else NA_real_
      }),
      site_cat = factor(site_cat, levels = c("2", "3", "4", "5+"))
    ) %>%
    filter(!is.na(site_cat), !is.na(stoich_pct_change))

  p <- ggplot(df_plot, aes(x = site_cat, y = stoich_pct_change)) +
    geom_boxplot(outlier.size = 0.6, fill = palette[k], alpha = 0.6, varwidth = TRUE) +
    labs(x = "Numero di siti m6A per trascritto", y = "Variazione relativa della stoichiometria DMSO (%)",
         title = paste0("Variazione della stoichiometria DMSO - ", k)) +
    theme_bw(base_size = 13)

  ggsave(paste0(k, "_stoichPercentageChangeBySiteCount.pdf"), p, width = 8, height = 8)
})

#IQR delle stoichiometrie dei trascritti con 5 o + siti
df_iqr <- lapply(names(tab), function(k) {
  tab[[k]] %>%
    mutate(
      n_sites = str_count(m6A_start_positions, ",") + 1,
      iqr_stoich = sapply(strsplit(m6A_stoich_dmso_mean, ","), function(v) {
        v <- as.numeric(v)
        if (sum(!is.na(v)) >= 2) IQR(v, na.rm = TRUE) else NA_real_
      }),
      cell_line = k
    ) %>%
    filter(n_sites >= 5, !is.na(iqr_stoich)) %>%
    dplyr::select(chrom, cell_line, n_sites, iqr_stoich)
}) %>%
  bind_rows()

p <- ggplot(df_iqr, aes(x = cell_line, y = iqr_stoich, fill = cell_line)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.6, varwidth = TRUE) +
  scale_fill_manual(values = palette) +
  labs(x = NULL,
       y = "IQR della stoichiometria DMSO per trascritto (%)",
       title = "Variabilità della stoichiometria tra siti m6A (trascritti con ≥5 siti)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave("stoich_IQR_per_transcript_5plus.pdf", p, width = 8, height = 8)

#Barplot
names_files_threshold <- list.files(path="./",pattern="threshold",recursive=TRUE)
files <- lapply(names_files_threshold,function(i){
            read.table(i,sep="\t")
        })
names(files) <- sapply(strsplit(names_files_threshold, "/", fixed = TRUE), `[`, 2)

cell_lines <- sub("_.*", "", names(files))
splitted <- split(files,cell_lines)

chiavi <- lapply(splitted, function(reps) {
  lapply(reps, function(df) distinct(df, V1, V2, V3))
})

conteggi <- lapply(chiavi, function(reps) {
  data.frame(
    n_replicati = length(reps),
    n_siti_union = nrow(reduce(reps, full_join, by = c("V1", "V2", "V3")))
  )
}) %>%
  bind_rows(.id = "cell_line")

coverages_dmso <- coverages_all[grep("DMSO",coverages_all$sample),]
coverage_dmso_sum <- coverages_dmso %>%
  group_by(cell_line) %>%
  summarise(
    total_coverage = sum(coverage, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

HCS_names <-list.files(path="./",pattern="High_condifence_sites.csv",recursive=TRUE)
HCS <- lapply(HCS_names,function(i){
		foe <- read.csv(i,header=TRUE)
		nrow(foe)

	})
names(HCS) <- sapply(strsplit(HCS_names, "/", fixed = TRUE), `[`, 1)


df_barplot <- data.frame(cbind("HCS"=unlist(HCS),"Union_sites"=conteggi$n_siti_union,"Coverage"=coverage_dmso_sum$total_coverage))
df_long <- df_barplot %>%
  rownames_to_column("cell_line") %>%
  mutate(cell_line = fct_reorder(cell_line, Coverage, .desc = TRUE)) %>%  # <-- ordinamento qui
  pivot_longer(-cell_line, names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("HCS", "Union_sites", "Coverage"))
  )
scale_factor <- max(df_barplot$Coverage) / max(df_barplot$HCS, df_barplot$Union_sites)

df_long <- df_long %>%
  mutate(value_plot = if_else(metric == "Coverage", value / scale_factor, value))

cols <- c("HCS" = "#1b9e77",
          "Union_sites" = "#d95f02",
          "Coverage" = "#7570b3")

p <- ggplot(df_long, aes(x = cell_line, y = value_plot, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_y_continuous(
    name = "Numero di siti (HCS, Union_sites)",
    labels = label_comma(),
    expand = expansion(mult = c(0, 0.05)),
    sec.axis = sec_axis(~ . * scale_factor,
                         name = "Coverage (reads)",
                         labels = label_comma())
  ) +
  scale_fill_manual(values = cols, name = NULL) +
  labs(x = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

ggsave("barplot_HCS_union_coverage_dualaxis.pdf", p, width = 9, height = 6)

################################################################################################################

#Which percetange of sites undergo a reduction of at least one fold?
HCS_log2FC <- lapply(HCS,function(i){
	dens = i
	dens$mean_DMSO <- rowMeans(dens[, grep("percent_modified_DMSO", colnames(dens))], na.rm = TRUE)
	dens$mean_STORM <- rowMeans(dens[, grep("percent_modified_STORM", colnames(dens))], na.rm = TRUE)
	dens$log2FC <- log2((dens$mean_STORM + 1) / (dens$mean_DMSO + 1))
	round(100*nrow(dens[dens$log2FC<(-1),])/nrow(i),2)
})
#BT483 87.17  #BT549 98.36  #MCF7 96.13  #MDAMB231 96.37  #SUM159 75.56  #T47D 95.81

count_sites <- function(x) {
  # gestisce eventuali NA o stringhe vuote
  ifelse(is.na(x) | x == "", 0, lengths(strsplit(x, ",")))
}

max_per_df <- sapply(tab, function(df) {
  n_siti <- count_sites(df$m6A_start_positions)
  max(n_siti, na.rm = TRUE)
 })

calc_fences <- function(x) {
  q <- quantile(x, probs = c(0.25,0.5, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  lower_fence <- q[1] - 1.5 * iqr
  upper_fence <- q[3] + 1.5 * iqr
  c(lower_fence = lower_fence, upper_fence = upper_fence,q1=q[1],q2=q[2],q3=q[3])
}

fences_per_df <- lapply(tab, function(df) {
  n_siti <- count_sites(df$m6A_start_positions)
  calc_fences(n_siti)
})

df_plot <- do.call(rbind, lapply(names(tab), function(nome) {
  n_siti <- count_sites(tab[[nome]]$m6A_start_positions)
  data.frame(
    cell_line = nome,
    n_siti = n_siti
  )
}))

stats_df <- df_plot %>%
  group_by(cell_line) %>%
  summarise(
    q1 = quantile(n_siti, 0.25, na.rm = TRUE),
    q2 = median(n_siti, na.rm = TRUE),
    q3 = quantile(n_siti, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower_fence = pmax(q1 - 1.5 * iqr, 0),   # non ha senso sotto 0
    upper_fence = q3 + 1.5 * iqr
  )

p <- ggplot(stats_df, aes(x = cell_line)) +
  geom_crossbar(aes(y = q2, ymin = q1, ymax = q3, fill = cell_line), width = 0.5, alpha = 0.6) +
  geom_errorbar(aes(ymin = lower_fence, ymax = upper_fence), width = 0.2) +
  geom_text(aes(y = q1, label = paste0("Q1=", q1)), hjust = 1.3, size = 3) +
  geom_text(aes(y = q2, label = paste0("Q2=", q2)), hjust = 1.3, size = 3, fontface = "bold") +
  geom_text(aes(y = q3, label = paste0("Q3=", q3)), hjust = 1.3, size = 3) +
  geom_text(aes(y = lower_fence, label = paste0("LF=", round(lower_fence,1))), hjust = -0.3, size = 3, color = "blue") +
  geom_text(aes(y = upper_fence, label = paste0("UF=", round(upper_fence,1))), hjust = -0.3, size = 3, color = "blue") +
  theme_minimal() +
  labs(title = "Statistics on number of sites", x = "Cell line", y = "n_siti") +
  theme(legend.position = "none")
ggsave("Boxplot_statistics_on_number_sites.pdf", p, width = 9, height = 6)


#Boxplot lunghezza dei trascritti con (1-2),(3-4),(5-10), (20 o +) siti
df_all <- lapply(names(tab), function(k){
  tab[[k]] %>%
    mutate(
      n_sites = str_count(m6A_start_positions, ",") + 1,
      site_cat = case_when(
        n_sites >= 1  & n_sites <= 2  ~ "1-2",
        n_sites >= 3  & n_sites <= 4  ~ "3-4",
        n_sites >= 5  & n_sites <= 10 ~ "5-10",
        n_sites >= 11 & n_sites <= 20 ~ "11-20",
        n_sites >= 21 ~ "21+"
      ),
      site_cat = factor(site_cat, levels = c("1-2", "3-4", "5-10", "11-20", "21+")),
      cell_line = k
    ) %>%
    left_join(tx_lengths, by = c("chrom" = "ensembl_transcript_id")) %>%
    filter(!is.na(transcript_length))
}) %>% bind_rows()

df_all$cell_line <- factor(df_all$cell_line, levels = names(tab))

p <- ggplot(df_all, aes(x = site_cat, y = log10(transcript_length), fill = cell_line)) +
  geom_boxplot(
    outlier.size = 0.6,
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(values = palette) +
  labs(
    x = "Numero di siti m6A per trascritto",
    y = "Log10(Lunghezza trascritto)",
    fill = "Cell line",
    title = "Lunghezza dei trascritti in funzione del numero di siti m6A, per cell line"
  ) +
  theme_bw(base_size = 13)

ggsave("all_cellLines_lengthOfTranscriptsDividedByRangesOfNumberSites.pdf", p, width = 12, height = 8)

#######################################################################################################
#Stratify the reduction in stoichiometry per number of m6A sites
HCS_names <-list.files(path="./",pattern="High_condifence_sites.csv",recursive=TRUE)
HCS <- lapply(HCS_names,function(i){
		foe <- read.csv(i,header=TRUE)
	})
names(HCS) <- sapply(strsplit(HCS_names, "/", fixed = TRUE), `[`, 1)

m6a_tables <- lapply(names(HCS), function(i){

	print(i)
  m6a_table <- HCS[[i]]
  m6a_table$chrom <- gsub("\\([+-]\\)$", "", m6a_table$chrom)

  m6a_table <- m6a_table %>%
    rowwise() %>%
    mutate(
      stoich_dmso_site  = sum(c_across(matches("^percent_modified_DMSO(_rep\\d+)?$")),  na.rm = TRUE),
      stoich_storm_site = sum(c_across(matches("^percent_modified_STORM(_rep\\d+)?$")), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    group_by(chrom) %>%
    arrange(start_position1, .by_group = TRUE) %>%
    summarise(
      m6A_start_positions = paste(start_position1, collapse = ","),
      total_stoich_dmso   = round(sum(stoich_dmso_site,  na.rm = TRUE), 2),
      total_stoich_storm  = round(sum(stoich_storm_site, na.rm = TRUE), 2),
      diff_dmso_storm     = round(total_stoich_dmso - total_stoich_storm, 2),
      .groups = "drop"
    )
})
names(m6a_tables) <- names(HCS)

delta_stoich_transcripts <- lapply(m6a_tables,function(i){
                i <- data.frame(i %>%
                left_join(tx_lengths, by = c("chrom" = "ensembl_transcript_id")) %>%
                filter(!is.na(transcript_length)))

})
names(delta_stoich_transcripts) <- names(m6a_tables)

df_all_reductions <- bind_rows(delta_stoich_transcripts, .id = "cell_line")
df_all_reductions$cell_line <- factor(df_all_reductions$cell_line, levels = names(m6a_tables))
palette <- palette[levels(df_all_reductions$cell_line)]

# Density plot
p <- ggplot(df_all_reductions, aes(x = diff_dmso_storm, color = cell_line)) +
  geom_density(linewidth = 0.9) +
  scale_color_manual(values = palette) +
  coord_cartesian(xlim = c(0, 5000)) +
  labs(
    x = "Differenza stoichiometria (DMSO - STORM)",
    y = "Densità",
    color = "Cell line",
    title = "Distribuzione della differenza di stoichiometria per trascritto"
  ) +
  theme_bw(base_size = 13)

ggsave("density_diff_dmso_storm_allCellLines_zoom.pdf", p, width = 9, height = 6)

#Boxplot of differences in stoichiometry divided per site number
df_all_delta <- lapply(names(delta_stoich_transcripts), function(k){
  delta_stoich_transcripts[[k]] %>%
    mutate(
      n_sites = str_count(m6A_start_positions, ",") + 1,
      site_cat = case_when(
        n_sites >= 1  & n_sites <= 2  ~ "1-2",
        n_sites >= 3  & n_sites <= 4  ~ "3-4",
        n_sites >= 5  & n_sites <= 10 ~ "5-10",
        n_sites >= 11 & n_sites <= 20 ~ "11-20",
        n_sites >= 21 ~ "21+"
      ),
      site_cat = factor(site_cat, levels = c("1-2", "3-4", "5-10", "11-20", "21+")),
      cell_line = k
    )
}) %>% bind_rows()

nuovo_ordine <- c("MCF7", "T47D", "BT483", "MDAMB231", "BT549", "SUM159")
df_all_delta$cell_line <- factor(df_all_delta$cell_line, levels = nuovo_ordine)
palette <- palette[nuovo_ordine]

p <- ggplot(df_all_delta, aes(x = site_cat, y = log2(diff_dmso_storm/n_sites), fill = cell_line)) +
  geom_boxplot(
    outlier.size = 0.6,
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(values = palette) +
  labs(
    x = "Numero di siti m6A per trascritto",
    y = "Differenza stoichiometria (DMSO - STORM)",
    fill = "Cell line",
    title = "Differenza di stoichiometria in funzione del numero di siti m6A, per cell line"
  ) +
  theme_bw(base_size = 13)

ggsave("boxplot_diff_dmso_storm_bySiteCat_log2_norm_per_site.pdf", p, width = 12, height = 8)

#Scatter reduction stochiometry vs DTE
tab_with_diff_stoich <- lapply(names(tab),function(i){
		left_join(tab[[i]],delta_stoich_transcripts[[i]],by=c("chrom","m6A_start_positions"))
	})
names(tab_with_diff_stoich) <- names(tab)

lapply(names(tab_with_diff_stoich), function(k){
  df <- tab_with_diff_stoich[[k]]
  
  p <- ggplot(df, aes(x = diff_dmso_storm, y = log2FC)) +
    geom_point(alpha = 0.4, size = 1) +
    facet_wrap(~ regulation) +
    scale_color_manual(values = colori_regulation) +
    labs(
      x = "Differenza stoichiometria (DMSO - STORM)",
      y = "log2FC(DTE)",
      title = paste0("log2FC vs differenza di stoichiometria, per classe di regolazione - ", k)
    ) +
    theme_bw(base_size = 13)
  
  ggsave(paste0(k, "_scatter_log2FC_vs_diffStoich_byRegulation.pdf"), p, width = 10, height = 4)
})


colori_regulation <- c("NS" = "grey70", "UP" = "red", "DOWN" = "skyblue")

lapply(names(tab_with_diff_stoich), function(k){
  df_plot <- tab_with_diff_stoich[[k]] %>%
  mutate(n_sites = str_count(m6A_start_positions, ",") + 1) %>%
    filter(!is.na(diff_dmso_storm)) %>%
    arrange(factor(regulation, levels = c("NS", "UP", "DOWN")))  # NS sotto, up/down sopra
  
  p <- ggplot(df_plot, aes(x = diff_dmso_storm/n_sites, y = log2FC, color = regulation)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = colori_regulation) +
    labs(
      x = "Differenza stoichiometria (DMSO - STORM)",
      y = "log2FC",
      color = "Regulation",
      title = paste0("log2FC vs differenza di stoichiometria - ", k)
    ) +
    theme_bw(base_size = 13)
  
  ggsave(paste0(k, "_scatter_log2FC_vs_diffStoich_byRegulation_colored_norm.pdf"), p, width = 8, height = 6)
})

#Gene ontology on up and down-regulated genes
GO <- lapply(names(tab_with_diff_stoich), function(k){
  
  df <- tab_with_diff_stoich[[k]]
  
  universe_genes <- df$gene_id %>% unique() %>% gsub("\\.\\d+$", "", .)

  up_genes <- df %>% 
    filter(regulation == "UP") %>% 
    pull(gene_id) %>% 
    unique() %>% 
    gsub("\\.\\d+$", "", .)   # rimuove la versione Ensembl (.8, .16, ecc.)
  
  down_genes <- df %>% 
    filter(regulation == "DOWN") %>% 
    pull(gene_id) %>% 
    unique() %>% 
    gsub("\\.\\d+$", "", .)
  
  ego_up <- enrichGO(
    gene = up_genes,
    universe = universe_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  ego_down <- enrichGO(
    gene = down_genes,
    universe = universe_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  list(up = ego_up, down = ego_down)
})
names(GO) <- names(tab_with_diff_stoich)


foe <- lapply(names(GO), function(k){
  
  cat("\n=== ", k, " - UP ===\n")
  if(!is.null(GO[[k]]$up) && nrow(GO[[k]]$up@result) > 0){
    res_up <- GO[[k]]$up@result %>% arrange(desc(p.adjust))
    print(res_up$Description)
  } else {
    cat("Nessun termine arricchito\n")
  }
  
  cat("\n=== ", k, " - DOWN ===\n")
  if(!is.null(GO[[k]]$down) && nrow(GO[[k]]$down@result) > 0){
    res_down <- GO[[k]]$down@result %>% arrange(desc(p.adjust))
    print(res_down$Description)
  } else {
    cat("Nessun termine arricchito\n")
  }
})

lapply(names(GO), function(k){
  
  if(!is.null(GO[[k]]$up) && nrow(GO[[k]]$up@result) > 0){
    p_up <- enrichplot::dotplot(GO[[k]]$up, showCategory = 15) + 
      ggtitle(paste0(k, " - UP regulated"))
    ggsave(paste0(k, "_GO_UP_dotplot.pdf"), p_up, width = 8, height = 8)
  }
  
  if(!is.null(GO[[k]]$down) && nrow(GO[[k]]$down@result) > 0){
    p_down <- enrichplot::dotplot(GO[[k]]$down, showCategory = 15) + 
      ggtitle(paste0(k, " - DOWN regulated"))
    ggsave(paste0(k, "_GO_DOWN_dotplot.pdf"), p_down, width = 8, height = 8)
  }
})


