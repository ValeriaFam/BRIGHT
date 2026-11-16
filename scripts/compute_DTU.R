# ================================================================================
# Script: compute_DTU.R
# Purpose: Compute DTU between cell lines
# Author: Camilla Ugolini
# ================================================================================

# -----------------------------------------------------
# Load argument parsing libraries
# -----------------------------------------------------
library(argparse)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 6) {
  stop("Six arguments must be provided")
}

# Assign arguments to variables
BASEDIR <- args[1]
quantification_dir <- args[2]
assembly_annotation_file <- args[3]
outputdir <- args[4]
ORF_annotation_file <- args[5]
SCRIPTDIR <- paste0(args[6],"/plots")
dir.create(paste0(outputdir,"/plots_and_tables/"), showWarnings = FALSE,recursive = TRUE)

# Load general configuration and functions
source(paste0(BASEDIR, "/BRIGHT_Dataset/general/config_R.R"))
source(paste0(BASEDIR, "/BRIGHT_Dataset/general/functions_R.R"))
source(paste0(BASEDIR, "/BRIGHT_Dataset/general/external_datasets.R"))

# ---------------
# Parameters
# ---------------

threads <- 8
FDR <- 0.05
cell_line_order <- c("SUM159_DMSO", "MDAMB231_DMSO","BT549_DMSO",
                     "MCF7_DMSO","T47D_DMSO","BT483_DMSO")
basal_cell_lines <- c("SUM159", "MDAMB231","BT549")
luminal_cell_lines <- c("MCF7","T47D","BT483")
subtype_order <- c("basal","luminal")

# ---------------
# Functions
# ---------------

# Source external functions
source(paste0(SCRIPTDIR, "/DGE_DTE_DTU_functions.R"))

# Function to compute DTU across all samples
DTU_processing <- function(print_tables_and_plots){
  
  res_list <- list()                                                              # empty list to store results
  kept_list <- list()
  
  quantification_list <- list.files(quantification_dir,pattern = "*_DMSO.transcript_grouped_counts.tsv",
                                    full.names = T,recursive = T)
  
  
  # Load each file once, not 15 times
  quantification_data <- lapply(quantification_list, function(x){
    read.table(x,header = T)
  })
  names(quantification_data) <- sub("\\..*", "", basename(quantification_list))
  s <- names(quantification_data)
  
  pairs <- c(lapply(s, function(x) c(x, x)),
             combn(s, 2, simplify = FALSE))
  length(pairs)   # 21
  
  
  # Perform differential analysis
  for (p in pairs) {
    cat(p[1], "vs", p[2], "\n")
    name1 <- p[1]
    name2 <- p[2]
    df1 <- quantification_data[[name1]]
    df2 <- quantification_data[[name2]]
    df1_sample_number <- length(names(df1 %>% dplyr::select(-gene_id)))
    df2_sample_number <- length(names(df2 %>% dplyr::select(-gene_id)))
    cat(df1_sample_number, "replicates vs", df2_sample_number,"replicates", "\n")
    comp_name <- paste0(name1, "_vs_", name2)        
    
    ###############################################
    ###### Prepare raw counts table ###############
    ###############################################
    tr_counts <- full_join(df1,df2,by="gene_id") %>%
      dplyr::rename("transcript_id"="gene_id") %>%
      left_join(transcript_gene_id,by="transcript_id") %>% 
      relocate(gene_id) %>%
      dplyr::rename("feature_id"="transcript_id")
    tr_counts[is.na(tr_counts)] <- 0
    
    cols <- setdiff(names(tr_counts), c("gene_id", "feature_id"))
    metadata <- data.frame(
      sample_id = cols,
      group = c(rep("Cell_line_1", df1_sample_number), rep("Cell_line_2", df2_sample_number))  
    )     
    
    d <- dmDSdata(counts = tr_counts, samples = metadata)                       # creating DRIMSeq data structure
    # histogram of the number of transcripts per gene
    p1 <- plotData(d)
    if (print_tables_and_plots == TRUE){
      ggsave(filename = paste0(outputdir,"/plots_and_tables","/number_of_transcripts_per_gene_",name1,"_vs_",name2,".pdf"), p1, 
             width = 10, height = 6)}
    
    
    ###############################################
    ###### Filtering ##############################
    ###############################################
    ## Check what is the minimal number of replicates per condition
    # this is useful for setting min_samps_feature_expr because then 
    # it allows us to detect features/transcripts that are expressed 
    # in one condition but not in the other one
    table(DRIMSeq::samples(d)$group)
    
    min_sample_gene <- sum(table(DRIMSeq::samples(d)$group))                    # total number of samples for which the minimum number of reads is required
    min_sample_tx <- min(table(DRIMSeq::samples(d)$group))                      # the smallest between the number of samples
    
    d <- dmFilter(d, 
                  min_samps_gene_expr = min_sample_gene, 
                  min_samps_feature_expr = min_sample_tx,
                  min_gene_expr = 10, 
                  min_feature_expr = 10)
    
    ################################################
    ##### Presicion estimation #####################
    ################################################
    ## dispersion = 1/(1 + precision)
    
    ## to parallelize 
    BPPARAM = BiocParallel::MulticoreParam(workers = threads)
    ## define design
    design_full <- model.matrix(~ group, data = DRIMSeq::samples(d))
    design_full
    ## To make the analysis reproducible (the dmPrecision takes just the 10% of the genes to make the estimation)
    set.seed(123)
    ## Calculate precision
    d <- dmPrecision(d, design = design_full, BPPARAM = BPPARAM)
    
    p2 <- plotPrecision(d)                                                      # QC control on dispersion 
    if (print_tables_and_plots == TRUE){
      ggsave(filename = paste0(outputdir,"/plots_and_tables","/precision_plot_",name1,"_vs_",name2,".pdf"), p2, 
             width = 10, height = 6)}
    
    ###################################
    ##### Tx Proporion estimation #####
    ###################################
    d <- dmFit(d, design = design_full, verbose = 1)
    
    cat("Print regression coefficients\n")
    print(head(proportions(d)))
    ## Get the DM regression coefficients (gene-level)
    print(head(coefficients(d)))
    ## Get the BB regression coefficients (feature-level)
    print(head(coefficients(d), level = "feature"))
    
    #################################################
    ### Testing for differential transcript usage ###
    #################################################
    design_null <- model.matrix(~ 1, data = DRIMSeq::samples(d))
    d <- dmTest(d, design = design_null)
    
    ## histogram of gene-level and transcript-level p-values
    p3 <- plotPValues(d)
    p4 <- plotPValues(d, level = "feature")
    if (print_tables_and_plots == TRUE){
      ggsave(filename = paste0(outputdir,"/plots_and_tables","/gene_level_pvalue_",name1,"_vs_",name2,".pdf"), p3, 
             width = 10, height = 6)}
    if (print_tables_and_plots == TRUE){
      ggsave(filename = paste0(outputdir,"/plots_and_tables","/transcript_level_pvalue_",name1,"_vs_",name2,".pdf"), p4, 
             width = 10, height = 6)}
    
    ## Store results
    res <- results(d)
    res.t <- results(d, level = "feature")
    
    ######################################################################
    ########## stageR to obtain padj on gene-level and Tx-level #########
    ######################################################################
    
    no.na <- function(x) ifelse(is.na(x), 1, x)
    res$pvalue <- no.na(res$pvalue)
    res.t$pvalue <- no.na(res.t$pvalue)
    
    ## Assign gene-level pvalues to the screening stage
    pScreen <- res$pvalue
    names(pScreen) <- res$gene_id
    ## Assign transcript-level pvalues to the confirmation stage
    pConfirmation <- matrix(res.t$pvalue, ncol = 1)
    rownames(pConfirmation) <- res.t$feature_id
    
    ## Create the gene-transcript mapping
    tx2gene <- res.t[, c("feature_id", "gene_id")]
    
    ## Run stageR
    stageRObj <- stageRTx(
      pScreen = pScreen,
      pConfirmation = pConfirmation,
      pScreenAdjusted = FALSE,
      tx2gene = tx2gene
    )
    
    stageRObj <- stageWiseAdjustment(
      stageRObj,
      method = "dtu",
      alpha = FDR
    )
    
    padj <- getAdjustedPValues(
      stageRObj,
      order = TRUE,
      onlySignificantGenes = FALSE
    )
    
    head(padj)
    ################################################################
    ## Derive DTU effect sizes from fitted proportions by DRIMSeq ##
    ################################################################
    
    prop <- DRIMSeq::proportions(d)
    sample_info <- DRIMSeq::samples(d)
    
    sample_info$sample_id <- as.character(sample_info$sample_id)
    
    # Replace these with your exact group names
    group_1 <- "Cell_line_1"
    group_2 <- "Cell_line_2"
    
    ids_1 <- sample_info$sample_id[sample_info$group == group_1]
    
    ids_2 <- sample_info$sample_id[sample_info$group == group_2]
    
    stopifnot(all(c(ids_1, ids_2) %in% colnames(prop)))
    
    # even without rowMeans, for each group we have the same value for one transcript and same for the other group
    prop_1 <- rowMeans(prop[, ids_1, drop = FALSE]) 
    prop_2 <- rowMeans(prop[, ids_2, drop = FALSE])
    
    effect_table <- data.frame(
      geneID = as.character(prop$gene_id),
      txID = as.character(prop$feature_id),
      
      group_1 = group_1,
      group_2 = group_2,
      
      proportion_group_1 = prop_1,
      proportion_group_2 = prop_2,
      
      # Signed: group 2 minus group 1
      delta_usage = prop_2 - prop_1,
      
      # Magnitude, regardless of direction
      abs_delta_usage = abs(prop_2 - prop_1)
    )
    
    head(effect_table)
    
    ########################################################################
    ### Combine effects table with DRIMSeq statistics and stageR results ###
    ########################################################################
    padj_df <- as.data.frame(padj)
    padj_df$geneID <- as.character(padj_df$geneID)
    padj_df$txID <- as.character(padj_df$txID)
    
    # Original transcript-level DRIMSeq statistics
    tx_res <- results(d, level = "feature")
    
    tx_stats <- data.frame(
      geneID = as.character(tx_res$gene_id),
      txID = as.character(tx_res$feature_id),
      tx_lr = tx_res$lr,
      tx_pvalue_raw = tx_res$pvalue
    )
    
    # Original gene-level DRIMSeq statistics
    gene_res <- results(d)
    
    gene_stats <- data.frame(
      geneID = as.character(gene_res$gene_id),
      gene_lr = gene_res$lr,
      gene_pvalue_raw = gene_res$pvalue
    )
    
    # Merge everything
    final_table <- merge(
      padj_df,
      effect_table,
      by = c("geneID", "txID"),
      all.x = TRUE,
      sort = FALSE
    )
    
    final_table <- merge(
      final_table,
      tx_stats,
      by = c("geneID", "txID"),
      all.x = TRUE,
      sort = FALSE
    )
    
    final_table <- merge(
      final_table,
      gene_stats,
      by = "geneID",
      all.x = TRUE,
      sort = FALSE
    )
    
    # stageR significance
    alpha <- getAlpha(stageRObj)
    
    final_table$gene_significant <- !is.na(final_table$gene) & final_table$gene <= alpha
    final_table$tx_significant <- !is.na(final_table$transcript) & final_table$transcript <= alpha
    final_table <- final_table %>%
      mutate(group_1 = name1) %>%
      mutate(group_2 = name2)
    
    ## prioritize results
    priority <- subset(
      final_table,
      tx_significant &
        abs(delta_usage) >= 0.10
    )
    
    priority_table <- priority[
      order(
        -abs(priority$delta_usage),
        priority$transcript
      ),
    ]
    priority_table <- priority_table %>%
      mutate(group_1 = name1) %>%
      mutate(group_2 = name2)
    
    length(unique(final_table[final_table$tx_significant==TRUE,"geneID"]))
    
    ## write tables
    if (print_tables_and_plots == TRUE){
      write.table(priority_table,paste0(outputdir,"/plots_and_tables/prioritized_DIUs_DRIMSeq_stageR_",name1,"_vs_",name2,".tsv"), 
                  sep = "\t", row.names = F, quote = F)
      write.table(final_table,paste0(outputdir,"/plots_and_tables/merged_DIUs_DRIMSeq_stageR_",name1,"_vs_",name2,".tsv"), 
                  sep = "\t", row.names = F, quote = F)
    }
    
    res_list[[comp_name]] <- final_table                        # store the results
    kept_list[[comp_name]] <- priority_table
  }
  return(list(results = kept_list, kept = res_list))
} 

# Function to compute DGE and DTE across grouped samples
DTU_all_processing <- function(print_tables_and_plots){
  
  quantification_list <- list.files(quantification_dir,pattern = "*_DMSO.transcript_grouped_counts.tsv",
                                    full.names = T,recursive = T)
  
  # Load each file once, not 15 times
  quantification_data <- lapply(quantification_list, function(x){
    read.table(x,header = T)
  })
  merged <- reduce(quantification_data, full_join, by = "gene_id")
  merged[is.na(merged)] <- 0
  merged <- merged %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "count") %>%
    mutate(cell_line = sub("_\\d+$", "", sample)) %>%
    group_by(gene_id, cell_line) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    pivot_wider(names_from = cell_line, values_from = count)
  basal <- merged %>% 
    dplyr::select(gene_id, starts_with(basal_cell_lines))
  luminal <- merged %>% 
    dplyr::select(gene_id, starts_with(luminal_cell_lines))
  
  # Perform differential analysis
  cat("basal vs luminal\n")
  name1 <- "basal"
  name2 <- "luminal"
  basal_sample_number <- length(names(basal %>% dplyr::select(-gene_id)))
  luminal_sample_number <- length(names(luminal %>% dplyr::select(-gene_id)))
  cat(basal_sample_number, "replicates vs", luminal_sample_number,"replicates", "\n")
  
  ###############################################
  ###### Prepare raw counts table ###############
  ###############################################
  tr_counts <- full_join(basal,luminal,by="gene_id") %>%
    dplyr::rename("transcript_id"="gene_id") %>%
    left_join(transcript_gene_id,by="transcript_id") %>% 
    relocate(gene_id) %>%
    dplyr::rename("feature_id"="transcript_id")
  tr_counts[is.na(tr_counts)] <- 0
  tr_counts <- as.data.frame(tr_counts)
  
  cols <- setdiff(names(tr_counts), c("gene_id", "feature_id"))
  metadata <- data.frame(
    sample_id = cols,
    group = c(rep("Basal", basal_sample_number), rep("Luminal", luminal_sample_number))  
  )     
  
  d <- dmDSdata(counts = tr_counts, samples = metadata)                       # creating DRIMSeq data structure
  
  # histogram of the number of transcripts per gene
  p1 <- plotData(d)
  if (print_tables_and_plots == TRUE){
    ggsave(filename = paste0(outputdir,"/plots_and_tables","/number_of_transcripts_per_gene_",name1,"_vs_",name2,".pdf"), p1, 
           width = 10, height = 6)}
  
  
  ###############################################
  ###### Filtering ##############################
  ###############################################
  ## Check what is the minimal number of replicates per condition
  # this is useful for setting min_samps_feature_expr because then 
  # it allows us to detect features/transcripts that are expressed 
  # in one condition but not in the other one
  table(DRIMSeq::samples(d)$group)
  
  min_sample_gene <- sum(table(DRIMSeq::samples(d)$group))                    # total number of samples for which the minimum number of reads is required
  min_sample_tx <- min(table(DRIMSeq::samples(d)$group))                      # the smallest between the number of samples
  
  d <- dmFilter(d, 
                min_samps_gene_expr = min_sample_gene, 
                min_samps_feature_expr = min_sample_tx,
                min_gene_expr = 10, 
                min_feature_expr = 10)
  
  ################################################
  ##### Presicion estimation #####################
  ################################################
  ## dispersion = 1/(1 + precision)
  
  ## to parallelize 
  BPPARAM = BiocParallel::MulticoreParam(workers = threads)
  ## define design
  design_full <- model.matrix(~ group, data = DRIMSeq::samples(d))
  design_full
  ## To make the analysis reproducible (the dmPrecision takes just the 10% of the genes to make the estimation)
  set.seed(123)
  ## Calculate precision
  d <- dmPrecision(d, design = design_full, BPPARAM = BPPARAM)
  
  p2 <- plotPrecision(d)                                                      # QC control on dispersion 
  if (print_tables_and_plots == TRUE){
    ggsave(filename = paste0(outputdir,"/plots_and_tables","/precision_plot_",name1,"_vs_",name2,".pdf"), p2, 
           width = 10, height = 6)}
  
  ###################################
  ##### Tx Proporion estimation #####
  ###################################
  d <- dmFit(d, design = design_full, verbose = 1)
  
  cat("Print regression coefficients\n")
  print(head(proportions(d)))
  ## Get the DM regression coefficients (gene-level)
  print(head(coefficients(d)))
  ## Get the BB regression coefficients (feature-level)
  print(head(coefficients(d), level = "feature"))
  
  #################################################
  ### Testing for differential transcript usage ###
  #################################################
  design_null <- model.matrix(~ 1, data = DRIMSeq::samples(d))
  d <- dmTest(d, design = design_null)
  
  ## histogram of gene-level and transcript-level p-values
  p3 <- plotPValues(d)
  p4 <- plotPValues(d, level = "feature")
  if (print_tables_and_plots == TRUE){
    ggsave(filename = paste0(outputdir,"/plots_and_tables","/gene_level_pvalue_",name1,"_vs_",name2,".pdf"), p3, 
           width = 10, height = 6)}
  if (print_tables_and_plots == TRUE){
    ggsave(filename = paste0(outputdir,"/plots_and_tables","/transcript_level_pvalue_",name1,"_vs_",name2,".pdf"), p4, 
           width = 10, height = 6)}
  
  ## Store results
  res <- results(d)
  res.t <- results(d, level = "feature")
  
  
  ######################################################################
  ########## stageR to obtain padj on gene-level and Tx-level #########
  ######################################################################
  
  no.na <- function(x) ifelse(is.na(x), 1, x)
  res$pvalue <- no.na(res$pvalue)
  res.t$pvalue <- no.na(res.t$pvalue)
  
  ## Assign gene-level pvalues to the screening stage
  pScreen <- res$pvalue
  names(pScreen) <- res$gene_id
  ## Assign transcript-level pvalues to the confirmation stage
  pConfirmation <- matrix(res.t$pvalue, ncol = 1)
  rownames(pConfirmation) <- res.t$feature_id
  
  ## Create the gene-transcript mapping
  tx2gene <- res.t[, c("feature_id", "gene_id")]
  
  ## Run stageR
  stageRObj <- stageRTx(
    pScreen = pScreen,
    pConfirmation = pConfirmation,
    pScreenAdjusted = FALSE,
    tx2gene = tx2gene
  )
  
  stageRObj <- stageWiseAdjustment(
    stageRObj,
    method = "dtu",
    alpha = FDR
  )
  
  padj <- getAdjustedPValues(
    stageRObj,
    order = TRUE,
    onlySignificantGenes = FALSE
  )
  
  head(padj)
  
  
  ################################################################
  ## Derive DTU effect sizes from fitted proportions by DRIMSeq ##
  ################################################################
  
  prop <- DRIMSeq::proportions(d)
  sample_info <- DRIMSeq::samples(d)
  
  sample_info$sample_id <- as.character(sample_info$sample_id)
  
  # Replace these with your exact group names
  group_1 <- "Basal"
  group_2 <- "Luminal"
  
  ids_1 <- sample_info$sample_id[sample_info$group == group_1]
  ids_2 <- sample_info$sample_id[sample_info$group == group_2]
  stopifnot(all(c(ids_1, ids_2) %in% colnames(prop)))
  
  # even without rowMeans, for each group we have the same value for one transcript and same for the other group
  prop_1 <- rowMeans(prop[, ids_1, drop = FALSE]) 
  prop_2 <- rowMeans(prop[, ids_2, drop = FALSE])
  
  effect_table <- data.frame(
    geneID = as.character(prop$gene_id),
    txID = as.character(prop$feature_id),
    
    group_1 = group_1,
    group_2 = group_2,
    
    proportion_group_1 = prop_1,
    proportion_group_2 = prop_2,
    
    # Signed: group 2 minus group 1
    delta_usage = prop_2 - prop_1,
    
    # Magnitude, regardless of direction
    abs_delta_usage = abs(prop_2 - prop_1)
  )
  
  head(effect_table)
  
  ########################################################################
  ### Combine effects table with DRIMSeq statistics and stageR results ###
  ########################################################################
  padj_df <- as.data.frame(padj)
  padj_df$geneID <- as.character(padj_df$geneID)
  padj_df$txID <- as.character(padj_df$txID)
  
  # Original transcript-level DRIMSeq statistics
  tx_res <- results(d, level = "feature")
  
  tx_stats <- data.frame(
    geneID = as.character(tx_res$gene_id),
    txID = as.character(tx_res$feature_id),
    tx_lr = tx_res$lr,
    tx_pvalue_raw = tx_res$pvalue
  )
  
  # Original gene-level DRIMSeq statistics
  gene_res <- results(d)
  
  gene_stats <- data.frame(
    geneID = as.character(gene_res$gene_id),
    gene_lr = gene_res$lr,
    gene_pvalue_raw = gene_res$pvalue
  )
  
  # Merge everything
  final_table <- merge(
    padj_df,
    effect_table,
    by = c("geneID", "txID"),
    all.x = TRUE,
    sort = FALSE
  )
  
  final_table <- merge(
    final_table,
    tx_stats,
    by = c("geneID", "txID"),
    all.x = TRUE,
    sort = FALSE
  )
  
  final_table <- merge(
    final_table,
    gene_stats,
    by = "geneID",
    all.x = TRUE,
    sort = FALSE
  )
  
  # stageR significance
  alpha <- getAlpha(stageRObj)
  
  final_table$gene_significant <- !is.na(final_table$gene) & final_table$gene <= alpha
  final_table$tx_significant <- !is.na(final_table$transcript) & final_table$transcript <= alpha
  final_table <- final_table %>%
    mutate(group_1 = name1) %>%
    mutate(group_2 = name2)
  
  ## prioritize results
  priority <- subset(
    final_table,
    tx_significant &
      abs(delta_usage) >= 0.10
  )
  
  priority_table <- priority[
    order(
      -abs(priority$delta_usage),
      priority$transcript
    ),
  ]
  priority_table <- priority_table %>%
    mutate(group_1 = name1) %>%
    mutate(group_2 = name2)
  
  length(unique(final_table[final_table$tx_significant==TRUE,"geneID"]))
  
  ## write tables
  if (print_tables_and_plots == TRUE){
    write.table(priority_table,paste0(outputdir,"/plots_and_tables/prioritized_DIUs_DRIMSeq_stageR_",name1,"_vs_",name2,".tsv"), 
                sep = "\t", row.names = F, quote = F)
    write.table(final_table,paste0(outputdir,"/plots_and_tables/merged_DIUs_DRIMSeq_stageR_",name1,"_vs_",name2,".tsv"), 
                sep = "\t", row.names = F, quote = F)
  }
  
  return(list(results = priority_table, kept = final_table))
} 

# -------------------------------
# Annotation datasets processing
# -------------------------------

# Load assembly annotation file
assembly_annotation <- rtracklayer::import(assembly_annotation_file)
transcripts_assembly <- assembly_annotation[assembly_annotation$type == "transcript"]
transcripts_assembly_df <- as.data.frame(transcripts_assembly) %>%
  dplyr::select(c("seqnames","start","end","width","strand","gene_id",
                  "gene_type","gene_name","level","tag","transcript_id",
                  "transcript_type","transcript_name")) 
transcript_gene_id <- transcripts_assembly_df %>%
  dplyr::select(c("gene_id","transcript_id")) 

genes_assembly <- assembly_annotation[assembly_annotation$type == "gene"]
genes_assembly_df <- as.data.frame(genes_assembly) %>%
  dplyr::select(c("gene_id","gene_type")) %>%
  mutate(gene_type = if_else(startsWith(gene_id, "ENSG"), gene_type, "Novel"))

# Transcripts with only one isoform per gene - sanity check
single_tx_genes <- transcript_gene_id %>%
  dplyr::count(gene_id, name = "n_tx") %>%
  dplyr::filter(n_tx == 1) %>%
  dplyr::pull(gene_id)

# Load ORF annotation
ORF_annotation <- read_tsv(ORF_annotation_file)

# Assess subtype per cell line
cline_group_map <- c(
  setNames(rep("Luminal", length(luminal_cell_lines)),
           paste0(luminal_cell_lines, "_DMSO")),
  setNames(rep("Basal", length(basal_cell_lines)),
           paste0(basal_cell_lines, "_DMSO"))
)

# ------------------------------------------------
# Compute DTU for all comparisons
# ------------------------------------------------

#----------- Results ------------------------------
# Compute DGE and DTE removing low expressed items
out_full <- DTU_processing(TRUE)
res_list_full  <- out_full$results     
kept_list_full <- out_full$kept 

#----------- Heatmap (all) ------------------------------
# Compute heatmap of DTU
DTU_heatmap <- heatmap_DTU_DRIMSeq("all_cell_lines","DTU","ENSG",
                                   res_list_full,kept_list_full,
                                   cell_line_order,
                                   biotype_col     = "biotype",
                                   biotype_value   = NULL,
                                   biotype_palette = GENCODE_biotypes_palette)

#----------- Heatmap (biotypes) ------------------------------
# Compute heatmap of DGE and DTE for protein coding and lncRNAs 
kept_list_full <- lapply(kept_list_full, function(df) {
  df <- df %>%
    left_join(ORF_annotation, by = c("txID" = "transcript_id")) %>%
    dplyr::select(-c(prediction_type,transcript_type_updated)) %>%
    left_join(genes_assembly_df,by = c("geneID" = "gene_id")) 
  df <- assign_gene_biotype(df,"geneID")
})


res_list_full <- lapply(res_list_full, function(df) {
  df <- df %>%
    left_join(ORF_annotation, by = c("txID" = "transcript_id")) %>%
    dplyr::select(-c(prediction_type,transcript_type_updated)) %>%
    left_join(genes_assembly_df,by = c("geneID" = "gene_id")) 
  df <- assign_gene_biotype(df,"geneID")
})

DTU_heatmap_pc <- heatmap_DTU_DRIMSeq("all_cell_lines","DTU","ENSG",
                              res_list_full,kept_list_full,cell_line_order,
                              biotype_col     = "gene_biotype",
                              biotype_value   = "Protein coding",
                              biotype_palette = GENCODE_biotypes_palette)
DTU_heatmap_lncRNA <- heatmap_DTU_DRIMSeq("all_cell_lines","DTU","ENSG",
                                  res_list_full,kept_list_full,cell_line_order,
                                  biotype_col     = "gene_biotype",
                                  biotype_value   = "lncRNA",
                                  biotype_palette = GENCODE_biotypes_palette)

# ------------------------------------------
# Compute DTU for basal vs luminal
# ------------------------------------------

#----------- Results ------------------------------
out_cell_type <- DTU_all_processing(TRUE)
res_list_cell_type  <- out_cell_type$results     
kept_list_cell_type <- out_cell_type$kept   

