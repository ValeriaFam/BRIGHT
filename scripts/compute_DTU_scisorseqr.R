# ================================================================================
# Script: compute_DTU_scisorseqr.R
# Purpose: Compute DTU with scisorseqR
# Author: Camilla Ugolini
# ================================================================================

# ---------------------
# Load libraries
# ---------------------

library("data.table")
library("rtracklayer")
library("scisorseqr")
library("ggnewscale")


# Assign arguments to variables
BASEDIR <- "/projects/cgs/cugolini/BC_atlas"
quantification_dir <- "/projects/cgs/cugolini/BC_atlas/analysis/all_cell_lines/dRNA/quantification_IsoQuant/all_cell_lines/IsoQuant_BRIGHT"
assembly_annotation_file <- "/projects/cgs/cugolini/BC_atlas/analysis/all_cell_lines/CFCseq/transcript_gene_intersection/reorder_gtf/filtered_corrected_assembly.sorted.gtf"
outputdir <- "/projects/cgs/cugolini/BC_atlas/analysis/all_cell_lines/dRNA/quantification_IsoQuant/all_cell_lines/DTE_DGE_DTU/DTU/scisorseqr/basal_vs_luminal"
out_prefix <- "IsoQuantOutput"
ORF_annotation_file <- "/projects/cgs/cugolini/BC_atlas/analysis/all_cell_lines/CFCseq/transcript_gene_intersection/ORF_annotator/assembly_gtf_CDS_coding_prediction.tsv"
SCRIPTDIR <- (paste0(BASEDIR,"/BRIGHT_Dataset/CFCseq/assembly/assemblies/assembly_quantification/dRNA_quantification/plots"))
dir.create(paste0(outputdir,"/",out_prefix), showWarnings = FALSE,recursive = TRUE)


# Load general configuration and functions
source(paste0(BASEDIR, "/BRIGHT_Dataset/general/config_R.R"))
source(paste0(BASEDIR, "/BRIGHT_Dataset/general/functions_R.R"))
source(paste0(BASEDIR, "/BRIGHT_Dataset/general/external_datasets.R"))

# Source external functions
source(paste0(SCRIPTDIR, "/DGE_DTE_DTU_functions.R"))


# ---------------
# Parameters
# ---------------

cell_line_order <- c("SUM159_DMSO", "MDAMB231_DMSO","BT549_DMSO",
                     "MCF7_DMSO","T47D_DMSO","BT483_DMSO")
basal_cell_lines <- c("SUM159", "MDAMB231","BT549")
luminal_cell_lines <- c("MCF7","T47D","BT483")
subtype_order <- c("basal","luminal")

# ---------------
# Run ScisorseqR
# ---------------

#----------- Read the GTF and build a transcript -> gene lookup ------------------------------

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

#----------- Aggregate all samples ------------------------------

quantification_list <- list.files(quantification_dir,pattern = "*_DMSO.transcript_grouped_counts.tsv",
                                    full.names = T,recursive = T)

quantification_data <- lapply(quantification_list, function(x){
  read.table(x,header = T)
})

merged <- quantification_data %>%
  reduce(full_join, by = "gene_id") %>%
  mutate(across(-gene_id, ~ replace(.x, is.na(.x), 0))) %>%
  dplyr::rename( feature_id=gene_id)

long_df <- merged %>%
  pivot_longer(
    cols = -feature_id,
    names_to = "group_id",
    values_to = "count"
  ) 

sample_order <- colnames(merged[,-1])
line <- sub("_[0-9]+$", "", sample_order)
rep  <- as.integer(sub(".*_", "", sample_order))

samples_sorted <- sample_order[order(match(line, cell_line_order), rep)]


#----------- Sum counts per transcript, attach gene id, filter ------------------------------

tot_df <- aggregate(count ~ feature_id, data = long_df, FUN = sum)
names(tot_df) <- c("tx", "n")

tot_df <- tot_df %>%
  left_join(transcript_gene_id, by = c("tx" = "transcript_id")) %>%
  dplyr::rename(gene=gene_id)


#----------- Rank isoforms within each gene (0 = most abundant), then filter further ---------

tot_df <- tot_df[order(tot_df$gene, -tot_df$n, tot_df$tx), ]         # sort first ...
tot_df$iso.id <- ave(tot_df$n, tot_df$gene, FUN = seq_along) - 1  # ... then number within each gene
tot_df$pos <- seq_len(nrow(tot_df))              # remembers the final gene/rank order for file 2
tot_df$key <- paste0(tot_df$gene, "::", tot_df$iso.id)


#----------- Write file 1: Gene  Isoform  n  iso.id ---------

file1 <- data.frame(Gene = tot_df$gene, Isoform = tot_df$tx, n = tot_df$n, iso.id = tot_df$iso.id)
write.table(file1, paste0(outputdir, "/",out_prefix,"/isoform_ranks.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)


#----------- Write file 2: GeneID::isoID  sample  count ---------

grid <- long_df
grid$pos   <- tot_df$pos[match(grid$feature_id, tot_df$tx)]
grid$key   <- tot_df$key[match(grid$feature_id, tot_df$tx)]
grid$lorder <- match(grid$group_id, sample_order)
grid <- grid[order(grid$pos, grid$lorder), ]
file2 <- data.frame(iso = grid$key, sample = grid$group_id, count = grid$count)
write.table(file2, paste0(outputdir, "/",out_prefix,"/NumIsoPerCluster"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("[done]", nrow(tot_df), "isoforms across", length(unique(tot_df$gene)), "genes and", length(samples_sorted), "samples\n")


#----------- write the config file for DiffSplicingAnalysis ---------

line <- sub("_DMSO_[0-9]+$", "", samples_sorted)
## group names should not contain underscores
config <- data.frame(
  c1 = "Basal",
  g1 = paste(samples_sorted[line %in% basal_cell_lines],   collapse = ","),
  c2 = "Luminal",
  g2 = paste(samples_sorted[line %in% luminal_cell_lines], collapse = ",")
)

write.table(config, paste0(outputdir, "/",out_prefix,"/config_Basal_vs_Luminal.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#----------- DiffSplicingAnalysis ---------

setwd(outputdir)
DiffSplicingAnalysis(paste0(outputdir, "/",out_prefix,"/config_Basal_vs_Luminal.tsv"),
                     numIsoforms = 10, minNumReads = 25,
                     typeOfTest = 'Iso', numThreads = 4)
writeLines(c("Basal", "Luminal"), "groups.txt")
triHeatmap("TreeTraversal_Iso/", "groups.txt",
           typeOfTest = "Iso", outName = "Heatmap_Basal_vs_Luminal")

#----------- Heatmap ---------

res_full <- read.table(paste0(outputdir, "/TreeTraversal_Iso/Basal_Luminal_10/Basal_Luminal_25_results.csv"),header = T)
res_sign <- read.table(paste0(outputdir, "/TreeTraversal_Iso/Basal_Luminal_10/sigGenes_Basal_Luminal25"),header = T)

res_sign <- res_sign %>%
    left_join(genes_assembly_df,by = c("x" = "gene_id")) %>%
  dplyr::rename(gene=x)
res_sign <- assign_gene_biotype(res_sign)

res_full <- res_full %>%
  left_join(genes_assembly_df,by = c("Gene" = "gene_id")) %>%
  dplyr::rename(gene=Gene)
res_full <- assign_gene_biotype(res_full)

heatmap_DTU(res_full,res_sign, "basal","luminal","DTU",
            annotation_pattern = "ENSG",
            biotype_col     = "biotype",
            biotype_value   = NULL,   # NULL = all biotypes
            biotype_palette = GENCODE_biotypes_palette)

heatmap_DTU(res_full,res_sign, "basal","luminal","DTU",
            annotation_pattern = "ENSG",
            biotype_col     = "gene_biotype",
            biotype_value   = "lncRNA",   # NULL = all biotypes
            biotype_palette = GENCODE_biotypes_palette)

heatmap_DTU(res_full,res_sign, "basal","luminal","DTU",
            annotation_pattern = "ENSG",
            biotype_col     = "gene_biotype",
            biotype_value   = "Protein coding",   # NULL = all biotypes
            biotype_palette = GENCODE_biotypes_palette)

