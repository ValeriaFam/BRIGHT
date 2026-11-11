# ================================================================================
# Script: external_datasets.R
# Purpose: Script containing relevant datasets for the manual analysis
# Author: Camilla Ugolini
# ================================================================================

# ------------------------------------------------------------------------------------------------------------
# Dataset of Breast Cancer relevant genes - somatic mutations (drivers) + germline mutations (by F.Nicassio)
# ------------------------------------------------------------------------------------------------------------

# ── Germline / predisposition lists (from before) ────────────────────────────
breast_cancer  <- c("BRCA1","BRCA2","PALB2","ATM","CHEK2","TP53","PTEN",
                    "BARD1","RAD51C","RAD51D","CDH1","STK11")

hrr_hrd        <- c("BRCA1","BRCA2","PALB2","ATM","CHEK2","BARD1",
                    "BRIP1","RAD51C","RAD51D")

mmr            <- c("MLH1","MSH2","MSH6","PMS2","EPCAM")

pan_cancer     <- c("APC","BAP1","CDKN2A","DICER1","FH","FLCN","MUTYH",
                    "NF1","POLD1","POLE","PTCH1","RB1","RET","SDHA",
                    "SDHAF2","SDHB","SDHC","SDHD","SMAD3","SMARCA4",
                    "SMARCB1","SUFU","TMEM127","TSC2","VHL")

# ── Somatic / breast cancer driver lists (new) ────────────────────────────────
core_recurrent       <- c("TP53","PIK3CA","GATA3")

luminal_enriched     <- c("PIK3CA","MAP3K1","GATA3","TP53","CDH1","MAP2K4",
                          "AKT1","RUNX1","CBFB","FOXA1","KMT2C","CDKN1B")

basal_tnbc           <- c("TP53","PIK3CA","BRCA1","BRCA2","ATM","RB1")

her2_enriched        <- c("ERBB2","TP53","PIK3CA","PIK3R1","EGFR","ERBB3")

pi3k_akt             <- c("PIK3CA","PIK3R1","PTEN","AKT1","AKT3","INPP4B",
                          "EGFR","ERBB2","IGF1R","KRAS","BRAF","FGFR1",
                          "FGFR2","KIT","MET","PDGFRA")

tp53_dna_damage      <- c("TP53","ATM","CHEK2","MDM2","MDM4","BRCA1","BRCA2")

rb1_cell_cycle       <- c("RB1","CCND1","CDK4","CDK6","CCNE1","CDKN1B",
                          "CDKN2A","CDKN2C")

basal_cnv            <- c("BRCA1","BRCA2","RB1","CCNE1","MYC","AKT3","TP53")

extended_smg         <- c("TBX3","RUNX1","CBFB","AFF2","PIK3R1","PTPN22",
                          "PTPRD","NF1","SF3B1","CCND3","CTCF","FOXA1")

# ── Union of all genes ────────────────────────────────────────────────────────
all_genes <- sort(unique(c(
  breast_cancer, hrr_hrd, mmr, pan_cancer,
  core_recurrent, luminal_enriched, basal_tnbc, her2_enriched,
  pi3k_akt, tp53_dna_damage, rb1_cell_cycle, basal_cnv, extended_smg
)))

germline <- c("breast_cancer_susc", "hrr_hrd", "mmr", "pan_cancer_syndromic")
somatic  <- c("core_recurrent", "luminal_enriched", "basal_tnbc", "her2_enriched",
              "pi3k_akt_pathway", "tp53_dna_damage", "rb1_cell_cycle",
              "basal_cnv", "extended_smg")

# ── Build dataframe ───────────────────────────────────────────────────────────
bc_gene_df <- data.frame(
  gene_name               = all_genes,
  # germline / predisposition
  breast_cancer_susc      = all_genes %in% breast_cancer,
  hrr_hrd                 = all_genes %in% hrr_hrd,
  mmr                     = all_genes %in% mmr,
  pan_cancer_syndromic    = all_genes %in% pan_cancer,
  # somatic / driver
  core_recurrent          = all_genes %in% core_recurrent,
  luminal_enriched        = all_genes %in% luminal_enriched,
  basal_tnbc              = all_genes %in% basal_tnbc,
  her2_enriched           = all_genes %in% her2_enriched,
  pi3k_akt_pathway        = all_genes %in% pi3k_akt,
  tp53_dna_damage         = all_genes %in% tp53_dna_damage,
  rb1_cell_cycle          = all_genes %in% rb1_cell_cycle,
  basal_cnv               = all_genes %in% basal_cnv,
  extended_smg            = all_genes %in% extended_smg,
  stringsAsFactors        = FALSE
) %>%
  dplyr::mutate(
    germline_susceptibility = dplyr::if_any(dplyr::all_of(germline)),
    somatic_driver          = dplyr::if_any(dplyr::all_of(somatic))
  )

