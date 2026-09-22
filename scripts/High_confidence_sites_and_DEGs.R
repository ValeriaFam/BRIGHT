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
library("parallel")

cell_line <- c("MCF7","BT483","T47D","SUM159","MDAMB231","BT549")
palette <- c("#CC6677","#882255","#AA4499","#117733","#999933","#44AA99")
n_cores <- 1
gtf <- import("/projects/CGS_shared/vfama/BRIGHT_PROJECT/filtered_corrected_assembly.gtf")
txdb <- makeTxDbFromGFF("/projects/CGS_shared/vfama/BRIGHT_PROJECT/filtered_corrected_assembly.gtf")
gtf_df <- as.data.frame(gtf) %>%
	  filter(type == "transcript") %>%
	  dplyr::select(transcript_id, gene_id,gene_type,transcript_type) %>%
	  distinct()

mclapply(seq_along(cell_line),function(k){

    path <-  paste0("/projects/CGS_shared/vfama/BRIGHT_PROJECT/",cell_line[k],"/")
    names_files <- list.files(path=path,pattern="Rep.*\\.csv$")
    files <- lapply(names_files,function(i){read.csv(paste0(path,i),header=TRUE)})
    names(files) <- sapply(strsplit(names_files, "\\.csv"), `[`, 1)

    ######################### DENSITY & KENDALL ################################

	df_wt <- lapply(files,function(i){
	    i[,c("chrom", "start_position1", "end_position1","Nvalid_cov_DMSO","percent_modified_DMSO")]
	  })

	#Prendere i siti comuni
	keys <- c("chrom", "start_position1", "end_position1")
	
	common_sites_btw_reps <- Reduce(function(x, y) {
	  i <- which(sapply(df_wt, identical, y))
	  inner_join(x, y, by = keys, suffix = c("", paste0("_rep", i)))
	}, df_wt)
	
	stoichiometry_wt <- common_sites_btw_reps[,grep("percent",colnames(common_sites_btw_reps))] 
	coverage_wt <- common_sites_btw_reps[, grep("Nvalid_cov", colnames(common_sites_btw_reps))]
	
	#UNIVERSO WT: tutte i siti >=20 in almeno un rep, #quanti siti 3 su 4
	n_sites_DMSO <- sum(rowSums(stoichiometry_wt >= 20, na.rm = TRUE) >= (length(df_wt)-1)) 
	row_filter <- rowSums(stoichiometry_wt >= 20, na.rm = TRUE) >= (length(df_wt)-1)
	names_stoich_wt <- common_sites_btw_reps[row_filter, c(1:3)]
	stoichiometry_wt <- stoichiometry_wt[row_filter, ]
	coverage_wt <- coverage_wt[row_filter, ] 

	#Metto ad NA gli zero 
	stoichiometry_wt[stoichiometry_wt==0] <- NA
	
	k_test_dmso <- kendall(stoichiometry_wt, correct = TRUE)
	
	res_kendall_dmso <- data.frame(
	  subjects = k_test_dmso$subjects,
	  value    = k_test_dmso$value,
	  p.value  = k_test_dmso$p.value
	)

	# Salva su csv
	write.csv(res_kendall_dmso, paste0(path,"kendall_results_dmso.csv"), row.names = FALSE)

	df_storm <- lapply(files,function(i){
	    i[,c("chrom", "start_position1", "end_position1","Nvalid_cov_STORM","percent_modified_STORM")]
	  })
	
	common_sites_btw_reps_STORM <- Reduce(function(x, y) {
	  i <- which(sapply(df_storm, identical, y))
	  inner_join(x, y, by = keys, suffix = c("", paste0("_rep", i)))
	}, df_storm)
	
	stoichiometry_storm <- common_sites_btw_reps_STORM[,grep("percent",colnames(common_sites_btw_reps_STORM))] 
	stoichiometry_storm <- stoichiometry_storm[rowSums(stoichiometry_storm > 0,na.rm=TRUE) >= (length(df_wt)-1),]
	
	k_test_storm <- kendall(stoichiometry_storm, correct = TRUE)
	res_kendall_storm <- data.frame(
	  subjects = k_test_storm$subjects,
	  value    = k_test_storm$value,
	  p.value  = k_test_storm$p.value
	)

	# Salva su csv
	write.csv(res_kendall_storm, paste0(path,"kendall_results_storm.csv"), row.names = FALSE)
	
	wt <- cbind(names_stoich_wt, coverage_wt, stoichiometry_wt)
	common_sites_btw_reps_all <- inner_join(common_sites_btw_reps_STORM, names_stoich_wt, by = keys)
	common_sites_btw_reps_all <- inner_join(wt,common_sites_btw_reps_all,by=keys)

	write.csv(common_sites_btw_reps_all,paste0(path,"High_condifence_sites.csv"), row.names = FALSE)

	k_test_combined <- kendall(common_sites_btw_reps_all[,grep("percent",colnames(common_sites_btw_reps_all))], correct = TRUE)
	res_kendall_combined <- data.frame(
	  subjects = k_test_combined$subjects,
	  value    = k_test_combined$value,
	  p.value  = k_test_combined$p.value
	)
	# Salva su csv
	write.csv(res_kendall_combined, paste0(path,"kendall_results_combined.csv"), row.names = FALSE)

	########## DENSITY DMSO & STORM ##############
	df_long_all <- common_sites_btw_reps_all %>%
  			rename_with(~ paste0(.x, "_rep1"), matches("^(percent_modified|Nvalid_cov)_(DMSO|STORM)$"))
  	df_long_all <- df_long_all %>%
  		pivot_longer(
  		  cols = -c(chrom, start_position1, end_position1),
  		  names_to = c(".value", "condition", "rep"),
  		  names_pattern = "(percent_modified|Nvalid_cov)_(DMSO|STORM)_rep(\\d+)"
    )
  	df_long_all$campione <- paste0(df_long_all$condition, "_rep", df_long_all$rep)

	df_long_all$campione <- factor(df_long_all$campione, levels = c(
	  "DMSO_rep1", "DMSO_rep2", "DMSO_rep3", "DMSO_rep4",
	  "STORM_rep1", "STORM_rep2", "STORM_rep3", "STORM_rep4"
	))

	p <- ggplot(df_long_all, aes(x = percent_modified, color = campione, linetype = condition)) +
  		geom_density(alpha = 0, fill = NA, linewidth = 0.8) +
  		scale_linetype_manual(values = c("DMSO" = "solid", "STORM" = "dashed")) +
  		labs(x = "Percent Modified", y = "Density", color = "Sample", linetype = "Treatment") +
  		theme_classic()
	ggsave(paste0(path,cell_line[k],"_density_DMSO_and_STORM_on_high_confidence_sites.pdf"), p, width = 8, height = 8)

	# #Inserire plot con sia wt che storm >=20

	# ########################## BOXPLOT ###############################
	dens <- common_sites_btw_reps_all

	dens$mean_DMSO <- rowMeans(dens[, grep("percent_modified_DMSO", colnames(dens))], na.rm = TRUE)
	dens$mean_STORM <- rowMeans(dens[, grep("percent_modified_STORM", colnames(dens))], na.rm = TRUE)
	dens$log2FC <- log2((dens$mean_STORM + 1) / (dens$mean_DMSO + 1))
	
	dens$bin <- cut(dens$mean_DMSO, 
	              breaks = seq(0, 100, by = 10),
	              labels = paste0(seq(0, 90, by = 10), "-", seq(10, 100, by = 10)),
	              include.lowest = TRUE)
	
	p <- ggplot(dens, aes(x = bin, y = as.numeric(log2FC))) +
	  geom_boxplot(fill = palette[k], alpha = 0.7, outlier.size = 0.5) +
	  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
	  theme_classic() +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	  xlab("DMSO stoichiometry bin (%)") +
	  ylab("log2FC (STORM/DMSO)")
	
	ggsave(filename=paste0(path,cell_line[k],"_Log2FC.pdf"),p)

	# ######################## METAGENE #################################

	# # Get UTR and CDS lengths per transcript
	# # utr5  <- fiveUTRsByTranscript(txdb,  use.names = TRUE)
	# # cds   <- cdsBy(txdb, by = "tx",     use.names = TRUE)
	# # utr3  <- threeUTRsByTranscript(txdb, use.names = TRUE)
	# # utr5_len <- sum(width(utr5))
	# # cds_len  <- sum(width(cds))
	# # utr3_len <- sum(width(utr3))
	# # tx_lengths <- tibble(
	# #   transcript_id = names(utr5_len),
	# #   utr5_len = as.numeric(utr5_len)
	# # ) %>%
	# #   full_join(tibble(transcript_id = names(cds_len),  cds_len  = as.numeric(cds_len)),  by = "transcript_id") %>%
	# #   full_join(tibble(transcript_id = names(utr3_len), utr3_len = as.numeric(utr3_len)), by = "transcript_id") %>%
	# #   replace_na(list(utr5_len = 0, cds_len = 0, utr3_len = 0))
	# # # --- 2. Assuming your df has columns: transcript_id, position_on_transcript ---
	# # # Merge with region lengths
	# # df_sites <- common_sites_btw_reps_all %>%
	# #   left_join(tx_lengths,  by = c("chrom" = "transcript_id")) %>%
	# #   mutate(
	# #     tx_len = utr5_len + cds_len + utr3_len,
	# #     # Classify each site into region
	# #     region = case_when(
	# #       start_position1 <= utr5_len                          ~ "5'UTR",
	# #       start_position1 <= utr5_len + cds_len                ~ "CDS",
	# #       start_position1 <= tx_len                            ~ "3'UTR",
	# #       TRUE ~ NA_character_
	# #     ),
	# #     # Normalize position within each region to [0, 1]
	# #     norm_position = case_when(
	# #       region == "5'UTR" ~ start_position1 / utr5_len,
	# #       region == "CDS"   ~ (start_position1 - utr5_len) / cds_len,
	# #       region == "3'UTR" ~ (start_position1 - utr5_len - cds_len) / utr3_len
	# #     ),
	# #     # Map to metagene scale: 5'UTR=[0,1], CDS=[1,2], 3'UTR=[2,3]
	# #     meta_position = case_when(
	# #       region == "5'UTR" ~ norm_position,
	# #       region == "CDS"   ~ 1 + norm_position,
	# #       region == "3'UTR" ~ 2 + norm_position
	# #     )
	# #   ) %>%
	# #   filter(!is.na(meta_position))

	# # #  # --- 3. Compute density per treatment ---
	# # # Pivot to get a "is this site in WT / STORM" column
	# # df_wt    <- df_sites %>% filter(rowMeans(!is.na(dplyr::select(., contains("percent_modified_DMSO"))))    > 0)
	# # df_storm <- df_sites %>% filter(rowMeans(!is.na(dplyr::select(., contains("percent_modified_STORM")))) > 0)
	
	# # df_density <- bind_rows(
	# #   df_sites %>% 
	# #     dplyr::select(meta_position, matches("percent_modified.*DMSO")) %>%
	# #     pivot_longer(-meta_position, names_to = "sample", values_to = "pct") %>%
	# #     filter(!is.na(pct), !is.na(meta_position)) %>%
	# #     mutate(treatment = "DMSO"),
	  
	# #   df_sites %>% 
	# #     dplyr::select(meta_position, matches("percent_modified.*STORM")) %>%
	# #     pivot_longer(-meta_position, names_to = "sample", values_to = "pct") %>%
	# #     filter(!is.na(pct), !is.na(meta_position)) %>%
	# #     mutate(treatment = "STORM")
	# # )
	
	# # # --- 4. Metagene plot ---
	# # p2 <- ggplot(df_density, aes(x = meta_position, color = treatment, linetype = treatment, group = sample)) +
	# #   geom_density(adjust = 0.5, linewidth = 1) +
	# #   scale_linetype_manual(values = c("DMSO" = "solid", "STORM" = "dashed")) +
	# #   scale_color_manual(values = c("DMSO" = "blue", "STORM" = "red")) +
	# #   scale_x_continuous(
	# #     breaks = c(0.5, 1.5, 2.5),
	# #     labels = c("5'UTR", "CDS", "3'UTR"),
	# #     limits = c(0, 3)
	# #   ) +
	# #   geom_vline(xintercept = c(1, 2), linetype = "dotted", color = "grey40") +
	# #   labs(x = "", y = "m6A site density",
	# #        color = "Treatment", linetype = "Treatment") +
	# #   theme_classic() +
	# #   theme(axis.text.x = element_text(size = 12, face = "bold"))
	
	# # ggsave(paste0(path,cell_line,"_metagene_m6A_density.pdf"), plot=p2,width = 8, height = 5)


	# ######################### TABLE HCS AND VOLCANO PLOTS ###########################
	
	#Load transcript files
	tr_tpm_dmso <- read.table(paste0("/projects/CGS_shared/vfama/BRIGHT_PROJECT/IsoQuant_BRIGHT/",cell_line[k],"_DMSO/",cell_line[k],"_DMSO.transcript_grouped_tpm.tsv"),header=TRUE)
	tr_tpm_stm <- read.table(paste0("/projects/CGS_shared/vfama/BRIGHT_PROJECT/IsoQuant_BRIGHT/",cell_line[k],"_STM/",cell_line[k],"_STM.transcript_grouped_tpm.tsv"),header=TRUE)
	tr_counts_dmso <- read.table(paste0("/projects/CGS_shared/vfama/BRIGHT_PROJECT/IsoQuant_BRIGHT/",cell_line[k],"_DMSO/",cell_line[k],"_DMSO.transcript_grouped_counts.tsv"),header=TRUE)
	tr_counts_stm <- read.table(paste0("/projects/CGS_shared/vfama/BRIGHT_PROJECT/IsoQuant_BRIGHT/",cell_line[k],"_STM/",cell_line[k],"_STM.transcript_grouped_counts.tsv"),header=TRUE)

	#Load gene files
	gene_tpm_dmso <- read.table(paste0("/projects/CGS_shared/vfama/BRIGHT_PROJECT/IsoQuant_BRIGHT/",cell_line[k],"_DMSO/",cell_line[k],"_DMSO.gene_grouped_tpm.tsv"),header=TRUE)
	gene_tpm_stm <- read.table(paste0("/projects/CGS_shared/vfama/BRIGHT_PROJECT/IsoQuant_BRIGHT/",cell_line[k],"_STM/",cell_line[k],"_STM.gene_grouped_tpm.tsv"),header=TRUE)
	gene_counts_dmso <- read.table(paste0("/projects/CGS_shared/vfama/BRIGHT_PROJECT/IsoQuant_BRIGHT/",cell_line[k],"_DMSO/",cell_line[k],"_DMSO.gene_grouped_counts.tsv"),header=TRUE)
	gene_counts_stm <- read.table(paste0("/projects/CGS_shared/vfama/BRIGHT_PROJECT/IsoQuant_BRIGHT/",cell_line[k],"_STM/",cell_line[k],"_STM.gene_grouped_counts.tsv"),header=TRUE)

	#Run DESeq2 on transcript raw counts
	tr_counts <- inner_join(tr_counts_dmso,tr_counts_stm,by="gene_id")
	rownames(tr_counts) <- tr_counts$gene_id
	tr_counts <- tr_counts[,-1]

	coldata <- data.frame(
	  row.names = colnames(tr_counts),
	  condition = factor(c(rep("DMSO", length(names_files)), rep("STM", length(names_files))),
	                     levels = c("DMSO", "STM"))  # DMSO = reference
	)
	dds <- DESeqDataSetFromMatrix(
	  countData = tr_counts,
	  colData   = coldata,
	  design    = ~ condition
	)
	dds <- DESeq(dds)
	res <- results(dds,
	               contrast = c("condition", "STM", "DMSO"),
	               alpha = 0.05)
	log2FC <- res$log2FoldChange
	names(log2FC) <- names(dds)
	pvals <- res$padj
	names(pvals) <- names(dds)

	#Volcano plot on transcripts
	res_df <- as.data.frame(res)
	res_df$gene <- rownames(res_df)
	#res_df <- na.omit(res_df) # Rimuovi NA (geni filtrati da DESeq2)
	
	res_df$status <- "NS" #Gene classification
	res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange >  1] <- "UP"
	res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "DOWN"
	res_df$status <- factor(res_df$status, levels = c("UP", "DOWN", "NS"))
	
	top_genes <- res_df[res_df$status != "NS", ] #Geni da etichettare (top 10 per padj)
	top_genes <- top_genes[order(top_genes$padj), ][1:10, ]
	
	volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
	  geom_point(alpha = 0.6, size = 1.5) +
	  scale_color_manual(values = c("UP" = "#E41A1C", "DOWN" = "#377EB8", "NS" = "grey70")) +
	  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.4) +
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.4) +
	  geom_text_repel(
	    data = top_genes,
	    aes(label = gene),
	    size = 3,
	    color = "black",
	    max.overlaps = 20
	  ) +
	  labs(
	    title = paste0(cell_line[k]," - DTE STM vs DMSO"),
	    x = "log2 Fold Change",
	    y = "-log10(padj)",
	    color = "Status"
	  ) +
	  theme_classic(base_size = 13)
	ggsave(filename=paste0(path,cell_line[k],"_Volcano_plot_on_transcript.pdf"),volcano)


	### Build the table
	#common_sites_btw_reps_all <- read.csv(paste0(path,"High_condifence_sites.csv"),header=TRUE)
	m6a_table <- common_sites_btw_reps_all
	m6a_table$chrom <- gsub("\\([+-]\\)$", "", m6a_table$chrom)
	m6a_table <- m6a_table %>%
	  rowwise() %>%
	  mutate(
	    stoich_dmso_mean = mean(c_across(matches("^percent_modified_DMSO(_rep\\d+)?$")), na.rm = TRUE),
	    stoich_storm_mean = mean(c_across(matches("^percent_modified_STORM(_rep\\d+)?$")), na.rm = TRUE),
	    diff_dmso_storm = stoich_dmso_mean - stoich_storm_mean
	  ) %>%
	  ungroup() %>%
	  group_by(chrom) %>%
	  summarise(
	    m6A_start_positions = paste(start_position1, collapse = ","),
	    m6A_stoich_dmso_mean = paste(round(stoich_dmso_mean, 2), collapse = ","),
	    m6A_diff_dmso_storm = paste(round(diff_dmso_storm, 2), collapse = ","),
	    .groups = "drop"
	  )
	colnames(tr_tpm_dmso)[-1] <- paste0("TMP_",colnames(tr_tpm_dmso)[-1])
	transcript_features <- data.frame("id"=names(log2FC),"log2FC"=log2FC,"p-val"=pvals,"regulation"=res_df$status)
	m6a_table_with_features <- m6a_table %>% inner_join(tr_tpm_dmso,by=c("chrom"="gene_id"))  %>% inner_join(transcript_features,by=c("chrom"="id"))

	
	# Join con tabella
	m6a_table_with_features_complete <- left_join(m6a_table_with_features, gtf_df,by = c("chrom" = "transcript_id"))

	### DESeq2 and volcano on gene counts
	gene_counts <- inner_join(gene_counts_dmso,gene_counts_stm,by="gene_id")
	rownames(gene_counts) <- gene_counts$gene_id
	gene_counts <- gene_counts[,-1]
	
	coldata_genes <- data.frame(
	  row.names = colnames(gene_counts),
	  condition = factor(c(rep("DMSO", length(names_files)), rep("STM", length(names_files))),
	                     levels = c("DMSO", "STM"))  # DMSO = reference
	)
	dds_genes <- DESeqDataSetFromMatrix(
	  countData = gene_counts,
	  colData   = coldata_genes,
	  design    = ~ condition
	)
	dds_genes <- DESeq(dds_genes)
	res_genes <- results(dds_genes,
	               contrast = c("condition", "STM", "DMSO"),
	               alpha = 0.05)
	log2FC_genes <- res_genes$log2FoldChange
	names(log2FC_genes) <- names(dds_genes)
	pvals_genes <- res_genes$padj
	names(pvals_genes) <- names(dds_genes)

	#Volcano plot
	res_df_genes <- as.data.frame(res_genes)
	res_df_genes$gene <- rownames(res_df_genes)
	#res_df <- na.omit(res_df) # Rimuovi NA (geni filtrati da DESeq2)
	
	res_df_genes$status <- "NS" #Gene classification
	res_df_genes$status[res_df_genes$padj < 0.05 & res_df_genes$log2FoldChange >  1] <- "UP"
	res_df_genes$status[res_df_genes$padj < 0.05 & res_df_genes$log2FoldChange < -1] <- "DOWN"
	res_df_genes$status <- factor(res_df_genes$status, levels = c("UP", "DOWN", "NS"))
	
	top_genes <- res_df_genes[res_df_genes$status != "NS", ] #Geni da etichettare (top 10 per padj)
	top_genes <- top_genes[order(top_genes$padj), ][1:10, ]
	
	volcano_genes <- ggplot(res_df_genes, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
	  geom_point(alpha = 0.6, size = 1.5) +
	  scale_color_manual(values = c("UP" = "#E41A1C", "DOWN" = "#377EB8", "NS" = "grey70")) +
	  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.4) +
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.4) +
	  geom_text_repel(
	    data = top_genes,
	    aes(label = gene),
	    size = 3,
	    color = "black",
	    max.overlaps = 20
	  ) +
	  labs(
	    title = paste0(cell_line[k]," - DEG STM vs DMSO"),
	    x = "log2 Fold Change",
	    y = "-log10(padj)",
	    color = "Status"
	  ) +
	  theme_classic(base_size = 13)
	ggsave(paste0(path,cell_line[k],"_Volcano_plot_on_genes.pdf"),volcano_genes)

	#Faccio il join con la tabella degli high confidence
	colnames(res_df_genes) <- paste0(colnames(res_df_genes),"_genes")
	tab_m6a_complete <- left_join(m6a_table_with_features_complete,
                     res_df_genes %>% dplyr::select(gene_genes, log2FoldChange_genes, padj_genes, status_genes),
                     by = c("gene_id" = "gene_genes"))

	write.csv(tab_m6a_complete,paste0(path,cell_line[k],"_sites_per_transcript_and_features.csv"),row.names=TRUE)

},mc.cores=n_cores)


















