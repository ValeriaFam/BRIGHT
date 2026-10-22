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
library("edgeR")
library("data.table")
library("Biostrings")
library("BSgenome")
library("AnnotationDbi")
library("org.Hs.eg.db")

cell_line <- c("MCF7","BT483","T47D","SUM159","MDAMB231","BT549")
palette <- c("#CC6677","#882255","#AA4499","#117733","#999933","#44AA99")
path <-  "/projects/CGS_shared/vfama/BRIGHT_PROJECT/"
names_files <- list.files(path=path,pattern="Rep.*\\.csv$",recursive=TRUE)

##############################################################################################
##############################             TABLE B         ###################################
##############################################################################################

#Create table with all sites from all cell lines with coverage greater than 20
key_cols <- c("chrom", "start_position1", "end_position1")
val_cols <- c("Nvalid_cov_DMSO", "percent_modified_DMSO",
              "Nvalid_cov_STORM", "percent_modified_STORM")

keys <- rbindlist(lapply(names_files, fread, select = key_cols))
keys_2 <- unique(keys)
setkeyv(keys_2, key_cols)

for (f in names_files){
	print(f)
  tag <- sub("\\.csv$", "", basename(f))              #es. "BT483_Rep1"
  dt  <- fread(f, select = c(key_cols, val_cols))
  setkeyv(dt, key_cols)
  new_names <- paste0(val_cols, "_", tag)             #es. "Nvalid_cov_DMSO_BT483_Rep1"
  keys_2[dt, (new_names) := mget(paste0("i.", val_cols)), on = key_cols]
  rm(dt); gc()
  message("fatto: ", tag)
}

keys_2[, genomic_strand := fifelse(grepl("\\([+-]\\)$", chrom),
                                   sub(".*\\(([+-])\\)$", "\\1", chrom),
                                   NA_character_)]
keys_2[, chrom := sub("\\([+-]\\)$", "", chrom)]

#Attach class, class A classification and cell line specific classification
classA_sites <- read.table("ClassA_sites_divided_by_specificity.csv",sep=",",header=TRUE)
classA_dt <- as.data.table(classA_sites)[, .(chrom,
                                            start = as.integer(m6A_start_positions),
                                            subtype)]

stopifnot(!anyDuplicated(classA_dt, by = c("chrom", "start")))

keys_2[classA_dt, Classification := i.subtype, on = c("chrom", "start_position1" = "start")]
keys_2[, class := fifelse(is.na(Classification), "B", "A")]

#Divide cell line specific sites in methylation driven or transcript driven
cov_min    <- 20L
stoich_max <- 20

groups <- list(
  Basal   = c("MDAMB231", "SUM159","BT549"),   
  Luminal = c("MCF7", "T47D","BT483")           
)

idx <- keys_2[class == "A" &
              (grepl("_specific$", Classification) |
               Classification %in% names(groups)), which = TRUE]
sp  <- keys_2[idx]
sp[, own := sub("_specific$", "", Classification)]

members_of <- function(own) if (own %in% names(groups)) groups[[own]] else own

classify_site <- function(i) {
  oth <- setdiff(cell_line, members_of(sp$own[i]))

  cov_all <- c(); pct_all <- c()
  for (cl in oth) {
    cov_cols <- grep(paste0("^Nvalid_cov_DMSO_", cl, "_Rep"), names(sp), value = TRUE)
    pct_cols <- sub("^Nvalid_cov_DMSO", "percent_modified_DMSO", cov_cols)
    cov_all  <- c(cov_all, unlist(sp[i, ..cov_cols]))
    pct_all  <- c(pct_all, unlist(sp[i, ..pct_cols]))
  }

  all_nocov <- all(is.na(cov_all) | cov_all < cov_min)
  all_low   <- all(!is.na(cov_all) & cov_all >= cov_min &
                   !is.na(pct_all) & pct_all < stoich_max)

  if (all_low) "methylation_driven"
  else if (all_nocov) "transcript_driven"
  else "mixed"
}

sp[, specific_type := vapply(seq_len(.N), classify_site, character(1))]

keys_2[, specific_type := NA_character_]
keys_2[idx, specific_type := sp$specific_type]

table(keys_2$Classification[idx], keys_2$specific_type[idx])

#Attach gene_id and gene_type
gtf <- import("/projects/CGS_shared/vfama/BRIGHT_PROJECT/filtered_corrected_assembly.gtf")
gtf_df <- as.data.frame(gtf) %>%
	  filter(type == "transcript") %>%
	  dplyr::select(transcript_id, gene_id,gene_type,transcript_type) %>%
	  distinct()
#in place conversion
setDT(keys_2) 
setDT(gtf_df)

keys_2[gtf_df, `:=`(gene_id = i.gene_id, gene_type = i.gene_type),
       on = c("chrom" = "transcript_id")]

fwrite(keys_2, file.path(path, "all_sites_joined.csv"))

#Find DRACH motifs
assembly_gr_ex <- subset(gtf, type == "exon")

genome_fa <- readDNAStringSet("GRCh38.primary_assembly.genome.fa")
names(genome_fa) <- sub("\\s.*", "", names(genome_fa))

chromosomes <- unique(as.character(seqnames(assembly_gr_ex)))
genome_fa_subset <- genome_fa[chromosomes]

ex <- assembly_gr_ex[order(mcols(assembly_gr_ex)$transcript_id,
                           ifelse(as.character(strand(assembly_gr_ex)) == "-",
                                  -start(assembly_gr_ex),
                                   start(assembly_gr_ex)))] #Exons in order of transcription: increasing for +, decreasing for -
assembly_gr_list <- split(ex, mcols(ex)$transcript_id)

txs_seqs <- extractTranscriptSeqs(genome_fa_subset, assembly_gr_list)

sites <- unique(keys_2[, .(chrom, pos = as.integer(start_position1) + 1L)]) # start_position1 is 0-based (bedMethyl), subseq is 1-based
sites <- sites[chrom %in% names(txs_seqs)]
sites[, tx_len := width(txs_seqs[chrom])]

ok <- sites$pos >= 3 & sites$pos + 2 <= sites$tx_len
sites[, kmer := NA_character_]
sites[ok, kmer := as.character(subseq(txs_seqs[chrom], start = pos - 2, end = pos + 2))]
sites[, DRACH := grepl("^[AGT][AG]AC[ACT]$", kmer)]
sites[is.na(kmer), DRACH := NA]

# Check
# str_dt <- unique(data.table(chrom = mcols(ex)$transcript_id,
#                             gtf_strand = as.character(strand(ex))))
# sites[str_dt, gtf_strand := i.gtf_strand, on = "chrom"]
# sites[!is.na(kmer), .(fracA = mean(substr(kmer, 3, 3) == "A"), n = .N),
#       by = gtf_strand]          # atteso ~0.99 per entrambi i filamenti

sites[, start_position1 := pos - 1L]
keys_2[sites, `:=`(kmer = i.kmer, DRACH = i.DRACH),
       on = c("chrom", "start_position1")]
keys_2[, motif := fifelse(DRACH %in% TRUE, kmer, NA_character_)] #attach also the motif

#Compute distance from the splice junction
ex_dt <- data.table(tx = as.character(mcols(ex)$transcript_id), w = width(ex))
ex_dt[, J := cumsum(w), by = tx]          # J = last exon base (1-based)

# chk <- ex_dt[, .(L = max(J)), by = tx] # check: last J of every transcript must be equal to the sequence length
# chk[, fa := as.integer(width(txs_seqs[tx]))]
# mean(chk$L == chk$fa)                      # atteso 1

junc <- ex_dt[, .SD[-.N, .(J)], by = tx] # junction = end of every exon except the last one
setkey(junc, tx, J)
junc[, Jx := J]
# --- Dist for every site ----------------------------------------------
# pos is 1-based; junction a monte: J <= pos - 1, a valle: J >= pos
sites[, J_up := junc[.(sites$chrom, sites$pos - 1L), on = .(tx, J), roll = Inf,
                     rollends = c(FALSE, FALSE)]$J]
sites[, J_down := junc[.(sites$chrom, sites$pos), on = .(tx, J), roll = -Inf,
                       rollends = c(FALSE, FALSE)]$J]
sites[, dist_up   := pos - J_up - 1L]      # nt dalla junction a monte
sites[, dist_down := J_down - pos]         # nt dalla junction a valle
sites[, dist_junc := pmin(dist_up, dist_down, na.rm = TRUE)]	#min distance
sites[is.na(J_up) & is.na(J_down), dist_junc := NA_integer_]   #monoexonic transcripts

keys_2[sites, `:=`(dist_junc = i.dist_junc, dist_up = i.dist_up, dist_down = i.dist_down),
       on = c("chrom", "start_position1")]

fwrite(keys_2, file.path(path, "all_sites_joined.csv"))

#Insert genomic positions
txdb <- makeTxDbFromGFF("/projects/CGS_shared/vfama/BRIGHT_PROJECT/filtered_corrected_assembly.gtf")
exons_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)

# Tieni solo i siti di trascritti presenti nel GTF
ok <- keys_2$chrom %in% names(exons_by_tx)
k  <- keys_2[ok, ]

# Posizione 1-based sul trascritto (vedi nota sotto)
sites_tx <- GRanges(
  seqnames = k$chrom,
  ranges   = IRanges(start = k$end_position1, width = 1)
)

gen <- mapFromTranscripts(sites_tx, exons_by_tx)

# xid = indice della riga di partenza
coords <- data.frame(
  idx        = gen$xHits,
  genome_chr = as.character(seqnames(gen)),
  genome_pos = start(gen),
  strand     = as.character(strand(gen))
)

k$idx <- seq_len(nrow(k))
k_gen <- dplyr::left_join(k, coords, by = "idx")

k_gen <- k_gen[,-ncol(k_gen)]
k_gen <- k_gen %>% dplyr::select(-idx)

fwrite(k_gen, file.path(path, "all_sites_joined.csv"))
write_xlsx(k_gen, path = "TableA.xlsx")
##############################################################################################
##############################             TABLES A        ###################################
##############################################################################################

#Definition of m6A score
# val_cols   <- c("Nvalid_cov_DMSO", "Nmod_DMSO", "Nvalid_cov_STORM", "Nmod_STORM")
# cols2      <- c("chrom", "start_position1", val_cols, "percent_modified_DMSO")
# stoich_min <- 20
# n_req_of   <- function(cl) if (cl == "T47D") 2L else 3L   # 3/4 rep, 2/3 per T47D

# cl_of_file  <- sub("_.*$", "", basename(names_files))
# rep_of_file <- sub("\\.csv$", "", basename(names_files))

# res <- lapply(unique(cl_of_file), function(cl) {
#   idx <- which(cl_of_file == cl)

#   #all reps of the same line
#   d <- rbindlist(lapply(idx, function(i) {
#     x <- fread(names_files[i], select = cols2)
#     x[, chrom := sub("\\([+-]\\)$", "", chrom)]
#     x[, replicate := rep_of_file[i]]
#     x
#   }))

#   #stoichiometry DMSO >= 20
#   d[, hi := !is.na(percent_modified_DMSO) & percent_modified_DMSO >= stoich_min]
#   n_files <- length(idx)
#   keep <- d[, .(n_rep = .N, n_hi = sum(hi)), by = .(chrom, start_position1)][
#           n_rep == n_files & n_hi >= n_req_of(cl), .(chrom, start_position1)]
#   d <- d[keep, on = .(chrom, start_position1), nomatch = 0L]

#   s <- d[, lapply(.SD, sum, na.rm = TRUE),
#          by = .(chrom, start_position1), .SDcols = val_cols]

#   a <- s[, c(lapply(.SD, sum),
#       list(n_sites       = .N,
#            n_sites_DMSO  = sum(Nvalid_cov_DMSO  > 0),
#            n_sites_STORM = sum(Nvalid_cov_STORM > 0))),
#   by = chrom, .SDcols = val_cols]
#   a[, cell_line := cl]

#   message("fatto: ", cl, " | siti tenuti: ", nrow(keep))
#   a
# })

# tx_scores <- rbindlist(res)
# tx_scores[, `:=`(
#   percent_DMSO  = fifelse(Nvalid_cov_DMSO  > 0, 100 * Nmod_DMSO  / Nvalid_cov_DMSO,  NA_real_),
#   percent_STORM = fifelse(Nvalid_cov_STORM > 0, 100 * Nmod_STORM / Nvalid_cov_STORM, NA_real_)
# )]

# tx_scores[, score := percent_DMSO - percent_STORM]
# tx_scores[, site_cat := fcase(
#   n_sites >= 1  & n_sites <= 2,  "1-2",
#   n_sites >= 3  & n_sites <= 4,  "3-4",
#   n_sites >= 5  & n_sites <= 10, "5-10",
#   n_sites >= 11 & n_sites <= 20, "11-20",
#   n_sites >= 21,                 "21+"
# )]
# tx_scores[, site_cat := factor(site_cat,
#                                levels = c("1-2", "3-4", "5-10", "11-20", "21+"))]

# palette <- palette[levels(tx_scores$cell_line)]
# nuovo_ordine <- c("MCF7", "T47D", "BT483", "MDAMB231", "BT549", "SUM159")
# tx_scores$cell_line <- factor(tx_scores$cell_line, levels = nuovo_ordine)
# palette <- palette[nuovo_ordine]

# fwrite(tx_scores, file.path(path, "transcript_m6A_scores.csv"))

# p <- ggplot(tx_scores[!is.na(score) & !is.na(site_cat)],
#             aes(x = site_cat, y = score, fill = cell_line)) +
#   geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
#   geom_boxplot(outlier.size = 0.6, alpha = 0.7,
#                position = position_dodge(width = 0.8)) +
#   scale_fill_manual(values = palette) +
#   labs(x = "Number of m6A per transcript",
#        y = "Score",
#        fill = "Cell line",
#        title = "Score as a function of number of site per transcript") +
#   theme_bw(base_size = 13)
# ggsave("New_deltascore_by_number_of_sites_per_transcript.pdf", p, width = 12, height = 8)

# p <- ggplot(tx_scores[!is.na(score) & !is.na(site_cat)],
#             aes(x = site_cat, y = percent_DMSO, fill = cell_line)) +
#   geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
#   geom_boxplot(outlier.size = 0.6, alpha = 0.7,
#                position = position_dodge(width = 0.8)) +
#   scale_fill_manual(values = palette) +
#   labs(x = "Number of m6A per transcript",
#        y = "Score",
#        fill = "Cell line",
#        title = "Score as a function of number of site per transcript") +
#   theme_bw(base_size = 13)
# ggsave("New_DMSOscore_by_number_of_sites_per_transcript.pdf", p, width = 12, height = 8)


# txdb <- makeTxDbFromGFF("/projects/CGS_shared/vfama/BRIGHT_PROJECT/filtered_corrected_assembly.gtf")
# tx_lengths_list <- exonsBy(txdb, by = "tx", use.names = TRUE)
# tx_lengths <- data.frame(
#   ensembl_transcript_id = names(tx_lengths_list),
#   transcript_length     = sum(width(tx_lengths_list))
# )

# tx_lengths <- as.data.table(tx_lengths)[, .(ensembl_transcript_id, transcript_length)]
# stopifnot(!anyDuplicated(tx_lengths, by = "ensembl_transcript_id"))
# tx_scores[tx_lengths, transcript_length := i.transcript_length,
#           on = c("chrom" = "ensembl_transcript_id")]


# p <- ggplot(tx_scores[!is.na(score) & !is.na(site_cat)],
#             aes(x = site_cat, y = score/transcript_length, fill = cell_line)) +
#   geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
#   geom_boxplot(outlier.size = 0.6, alpha = 0.7,
#                position = position_dodge(width = 0.8)) +
#   scale_fill_manual(values = palette) +
#   labs(x = "Number of m6A per transcript",
#        y = "Score",
#        fill = "Cell line",
#        title = "Score as a function of number of site per transcript normalized over the length of the transcript") +
#   theme_bw(base_size = 13)
# ggsave("New_deltascore_by_number_of_sites_per_transcript_normalized.pdf", p, width = 12, height = 8)

# p <- ggplot(tx_scores[!is.na(score) & !is.na(site_cat)],
#             aes(x = site_cat, y = percent_DMSO/transcript_length, fill = cell_line)) +
#   geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
#   geom_boxplot(outlier.size = 0.6, alpha = 0.7,
#                position = position_dodge(width = 0.8)) +
#   scale_fill_manual(values = palette) +
#   labs(x = "Number of m6A per transcript",
#        y = "Score",
#        fill = "Cell line",
#        title = "Score as a function of number of site per transcript") +
#   theme_bw(base_size = 13)
# ggsave("New_DMSOscore_by_number_of_sites_per_transcript_normalized.pdf", p, width = 12, height = 8)


#Making the table
tx_scores <- read.table("transcript_m6A_scores.csv",sep=",",header=TRUE)

HC_sites_names <- list.files(path=path,pattern="_sites_per_transcript_and_features.csv",recursive=TRUE)
HC_sites <- lapply(HC_sites_names,function(i){
              read.table(paste0(path,i),sep=",",header=TRUE)
       })
names(HC_sites) <- sapply(strsplit(HC_sites_names, "/"), `[`, 1)

merged_list <- imap(HC_sites, function(df, cl) {

  # 1. Select columns from HC_sites
  base <- df %>%
    dplyr::select(chrom, gene_id, starts_with("TMP_"),m6A_start_positions,m6A_stoich_dmso_mean,m6A_diff_dmso_storm)   # prende tutti i TMP_*_DMSO_*

  # 2. Takes rows of that cell line from the second df
  m6a <- tx_scores %>%
    filter(cell_line == cl) %>%
    transmute(
      chrom,
      m6a_load   = n_sites,
      m6a_score  = percent_DMSO,
      m6a_score_kb = percent_DMSO / (transcript_length / 1000),
      delta_score  = score,
      tx_length = transcript_length
    )

  # 3. Merge by chrom
  base %>%
    inner_join(m6a, by = "chrom") %>%   # left_join se vuoi tenere anche i trascritti senza m6A
    mutate(cell_line = cl)
})

merged_list <- map(merged_list, function(df) {
  df %>%
    mutate(gene_id_clean = sub("\\..*$", "", gene_id),
           gene_name = mapIds(org.Hs.eg.db,
                              keys = gene_id_clean,
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")) %>%
    dplyr::select(-gene_id_clean)
})
saveRDS(merged_list, file = "TableB.rds")

write_xlsx(merged_list, path = "TableB.xlsx")

#import external_dataset.R to see BC relevant genes
load("scripts/external_datasets.R")

