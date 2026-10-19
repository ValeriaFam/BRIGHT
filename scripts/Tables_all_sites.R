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

cell_line <- c("MCF7","BT483","T47D","SUM159","MDAMB231","BT549")
palette <- c("#CC6677","#882255","#AA4499","#117733","#999933","#44AA99")
path <-  "/projects/CGS_shared/vfama/BRIGHT_PROJECT/"
names_files <- list.files(path=path,pattern="Rep.*\\.csv$",recursive=TRUE)

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
n_low_req  <- function(cl) if (cl == "T47D") 1L else 2L

idx <- keys_2[class == "A" & grepl("_specific$", Classification), which = TRUE] #subset on class A sites
sp  <- keys_2[idx]
sp[, own := sub("_specific$", "", Classification)]

nocov_m <- low_m <- matrix(FALSE, nrow(sp), length(cell_line),
                           dimnames = list(NULL, cell_line))
for (cl in cell_line) {
  cov_cols <- grep(paste0("^Nvalid_cov_DMSO_", cl, "_Rep"), names(sp), value = TRUE)
  pct_cols <- sub("^Nvalid_cov_DMSO", "percent_modified_DMSO", cov_cols)
  cov <- as.matrix(sp[, ..cov_cols]); pct <- as.matrix(sp[, ..pct_cols])

  cov_ok <- !is.na(cov) & cov >= cov_min
  low    <- cov_ok & !is.na(pct) & pct < stoich_max

  nocov_m[, cl] <- rowSums(cov_ok) == 0
  low_m[, cl]   <- rowSums(low) >= n_low_req(cl)
}

own_m <- outer(sp$own, cell_line, "==")
n_other <- length(cell_line) - 1L

all_low   <- rowSums(low_m   & !own_m) == n_other
all_nocov <- rowSums(nocov_m & !own_m) == n_other

sp[, specific_type := fifelse(all_low, "methylation_driven",
                      fifelse(all_nocov, "transcript_driven", "mixed"))]
keys_2[, specific_type := NA_character_]
keys_2[idx, specific_type := sp$specific_type]

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








