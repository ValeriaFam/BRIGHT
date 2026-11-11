# ==============================================================
# Script: methods_isoform_switch.R
# Purpose: Script to find isoform switches based on delta for 
# both scisorseqR and DRIMSeq
# Author: Fatemeh Kordevani
# ==============================================================

# -----------------------------------------------------
# Load necessary libraries
# -----------------------------------------------------
BASEDIR <- "/projects/CGS_shared/fkordevani/BC_atlas"
source(paste0(BASEDIR, "/BRIGHT_Dataset/general/config_R.R"))
source(paste0(BASEDIR, "/BRIGHT_Dataset/general/functions_R.R"))
cell_line <- c("MCF7","BT483","T47D","SUM159","MDAMB231","BT549")
# -------------------------------------
# Import and process files (scISOrseqR)
# -------------------------------------
## Load dataframe output that contains all isoforms of a given gene, not only the top two
for (cell in cell_line){
  pairwise_dir <- paste0("/projects/CGS_shared/vfama/BRIGHT_PROJECT/DIU/",cell,"TreeTraversal_Iso/DMSO_STM_10/")
  all_iso_files <- list.files(pairwise_dir, recursive = TRUE, pattern = ".Robj$")
  all_iso_files_lst <- list()
  for (all_iso_file in all_iso_files) {
    load(file.path(pairwise_dir, all_iso_file))
    all_iso_files_lst[[paste0("all_iso_",all_iso_file)]] <- processedDF
    processedDF <- NULL
  }
  ## keep the name simple and put list items (DFs) in the R environment
  names(all_iso_files_lst) <- sub("/.*", "", names(all_iso_files_lst))
  list2env(all_iso_files_lst, envir = .GlobalEnv)

  ## load results file and perform GO analysis
  results_files <- list.files(pairwise_dir, recursive = TRUE, pattern = "results.csv$")
  signif_res_files_lst <- list()
  for (file in results_files) {
    name <- sub("_results\\.csv$", "", basename(file))
    result_file <- read.table(file.path(pairwise_dir, file), sep = "\t", header = TRUE)
    signif <- result_file[result_file$FDR <= 0.05 & abs(result_file$dPI) >= 0.1, ]
    signif_res_files_lst[[paste0("sigres_",name)]] <- signif
  }
  list2env(signif_res_files_lst, envir = .GlobalEnv)

  ### Load file that matches isoform numbers of each gene to its TranscriptID
  ranking_isoforms <- read.table("/projects/CGS_shared/vfama/BRIGHT_PROJECT/DIU/IsoQuantOutput/isoform_ranks.tsv"), header = T)
  
  ### load gtf that was used for quantfication
  gtf_file <- "/projects/CGS_shared/vfama/BRIGHT_PROJECT/filtered_corrected_assembly.gtf"
  colnames_GTF <- c("seqname", "source", "feature","start","end","score","strand","frame","attributes")
  gtf <- fread(gtf_file, sep = "\t", header = FALSE)
  colnames(gtf) <- colnames_GTF
  head(gtf)
  gtf_exons <- gtf[gtf$feature == "exon",]
  
  # extract desired columns
  gtf_exons$transcript_id <- str_extract(gtf_exons$attributes, 'transcript_id "[^"]+"')
  gtf_exons$transcript_id <- str_replace_all(gtf_exons$transcript_id, 'transcript_id "|"','')
  
  ##############################################################
  ########## merging different datasets to obtain transcriptID, 
  ########## transcript biotype, isoform number and deltaPI
  ##############################################################

  ## we want to replace iso.id with transcript id in the signif res df
  iso_sigres_basalvsLum <- sigres_Basal_Luminal_25 %>%
    left_join(
      ranking_isoforms %>%
        dplyr::select(Gene, iso.id, Isoform) %>%
        rename(
          transcriptID_ix1 = Isoform,
        ),
      by = c("Gene", "maxDeltaPI_ix1" = "iso.id")
    ) %>%
    left_join(
      ranking_isoforms %>%
        dplyr::select(Gene, iso.id, Isoform) %>%
        rename(
          transcriptID_ix2 = Isoform,
        ),
      by = c("Gene", "maxDeltaPI_ix2" = "iso.id")
    )
  
  head(iso_sigres_basalvsLum)

### all isofroms with biotype and transcriptID and protein ID
# we want to find the first basal-favored and luminal-favored in each gene, then see whether they produce different proteins 
all_iso_Basal_Luminal_10 <- left_join(all_iso_Basal_Luminal_10, ranking_isoforms, 
                                           by = c("Gene", "IsoID" = "iso.id"))
all_iso_Basal_Luminal_10$ENSEMBL <- sub("\\..*$","",all_iso_Basal_Luminal_10$Gene)
head(all_iso_Basal_Luminal_10)
# only in significant genes not all genes, significant genes can be found in sigres_Basal_Luminal_25
sig_genes_basalvsLum <- unique(sigres_Basal_Luminal_25$Gene)
all_iso_Basal_Luminal_10_signif <- all_iso_Basal_Luminal_10[all_iso_Basal_Luminal_10$Gene %in% sig_genes_basalvsLum,]
## make cols homgeneous among scisorseqr and drimseq
scisors_cols_to_keep <- c("Gene", "Basal", "Luminal", "pi1", "pi2", "delta", "Isoform")
all_iso_Basal_Luminal_10_signif <- all_iso_Basal_Luminal_10_signif[,scisors_cols_to_keep]
colnames(all_iso_Basal_Luminal_10_signif) <- c("Gene", "basal_reads", "luminal_reads", "basal_prop", "luminal_prop", "delta", "Isoform")
# find top basal-favored transcripts
basal_Tx <- all_iso_Basal_Luminal_10_signif %>%
  filter(
    delta > 0,
    !is.na(Isoform)
  ) %>%
  group_by(Gene) %>%
  slice_max(
    order_by = delta,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  dplyr::select(
    Gene,
    basal_Isoform = Isoform,
    basal_reads,
    basal_prop,
    basal_delta = delta)

# find top luminal-favored transcripts
luminal_Tx <- all_iso_Basal_Luminal_10_signif %>%
  filter(
    delta < 0,
    !is.na(Isoform)
  ) %>%
  group_by(Gene) %>%
  slice_min(
    order_by = delta,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  dplyr::select(
    Gene,
    luminal_Isoform = Isoform,
    luminal_reads,
    luminal_prop,
    luminal_delta = delta
  )

Tx_switches_scrisors <- inner_join(
  basal_Tx,
  luminal_Tx,
  by = "Gene"
)


Tx_switches_scrisors <- Tx_switches_scrisors %>%
  filter(
    basal_reads >= 5,
    luminal_reads >= 5
  )
Tx_switches_scrisors$basal_reads <- NULL
Tx_switches_scrisors$luminal_reads <- NULL
Tx_switches_scrisors$source <- "scISOrSeq"
# -------------------------------------
# Import and process files (DRIMSeq)
# -------------------------------------

drimseq_priority <- read.table(paste0(BASEDIR,"/analysis/all_cell_lines/dRNA/DIFFERENTIAL/DIU_DRIMSeq/prioritized_DIUs_DRIMSeq_stageR_Basal_vs_Luminal.tsv"), header=T)
drims_cols_to_keep <- c("geneID", "txID", "proportion_group_1", "proportion_group_2", "delta_usage")
drimseq_priority <- drimseq_priority[,drims_cols_to_keep]
### make colnames similar across scisorseqr and drimseq tables
colnames(drimseq_priority) <- c("Gene", "Isoform", "luminal_prop", "basal_prop", "delta")

# find top basal-favored transcripts
basal_Tx_drimseq <- drimseq_priority %>%
  filter(
    delta > 0,
    !is.na(Isoform)
  ) %>%
  group_by(Gene) %>%
  slice_max(
    order_by = delta,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  dplyr::select(
    Gene,
    basal_Isoform = Isoform,
    basal_prop,
    basal_delta = delta)

# find top luminal-favored transcripts
luminal_Tx_drimseq <- drimseq_priority %>%
  filter(
    delta < 0,
    !is.na(Isoform)
  ) %>%
  group_by(Gene) %>%
  slice_min(
    order_by = delta,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  dplyr::select(
    Gene,
    luminal_Isoform = Isoform,
    luminal_prop,
    luminal_delta = delta
  )

Tx_switches_drimseq <- inner_join(
  basal_Tx_drimseq,
  luminal_Tx_drimseq,
  by = "Gene"
)

Tx_switches_drimseq$source <- "DRIMSeq"


# ----------------------------------------------
# Merge switches from DRIMseq and scisorseqr
# ----------------------------------------------

colnames(Tx_switches_scrisors)
colnames(Tx_switches_drimseq)

combined_DRIMSeq_scISOrSeq <- rbind(Tx_switches_scrisors, Tx_switches_drimseq)

write.table(combined_DRIMSeq_scISOrSeq,paste0(BASEDIR,"/analysis/all_cell_lines/dRNA/DIFFERENTIAL/DIU_Combined/signif_isoform_switch_DRIMSeq_scISOrSeq_basal_vs_luminal.tsv"), sep = "\t", row.names = F)


}

 