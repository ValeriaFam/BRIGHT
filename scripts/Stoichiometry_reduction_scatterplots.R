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
library("writexl")
library("parallel")

cell_line <- c("MCF7","BT483","T47D","SUM159","MDAMB231","BT549")
palette <- c("#CC6677","#882255","#AA4499","#117733","#999933","#44AA99")
n_cores <- 6

mclapply(seq_along(cell_line),function(k){

    path <-  paste0("/projects/CGS_shared/vfama/BRIGHT_PROJECT/",cell_line[k],"/")
    names_files <- list.files(path=path,pattern="filtCov.m6A.bed",recursive=TRUE)
    files <- lapply(names_files,function(i){
            read.table(paste0(path,i),sep="\t")
        })
    names(files) <- unname(sapply(names_files,function(i){
        sub("/.*", "", i)
        }))
    
    mod_sites <- lapply(files[grep("DMSO",names(files))],function(i){
            colnames(i) <- c("chrom","start_position1","end_position1","modified_base" ,"score","strand","start_position2","end_position2","color","Nvalid_cov","percent_modified","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall")
            i[i$Nvalid_cov>=20,]
        })
    names(mod_sites) <- names(files[grep("DMSO",names(files))])
    
    mod_sites_storm <- lapply(files[grep("STM",names(files))],function(i){
            colnames(i) <- c("chrom","start_position1","end_position1","modified_base" ,"score","strand","start_position2","end_position2","color","Nvalid_cov","percent_modified","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall")
            i[i$Nvalid_cov>=20,]
        })
    names(mod_sites_storm) <- names(files[grep("STM",names(files))])
    
    
    #DF for plotting the stoichiometry curve -  all sites common and non
    df_all <- lapply(seq_along(mod_sites),function(i){
        tbl <- full_join(
                mod_sites[[i]],
                mod_sites_storm[[i]],
                by = c("chrom", "start_position1", "end_position1"),
                suffix = c("_DMSO", "_STORM"))
        data.table::fwrite(tbl, paste0(path,cell_line[k],"_Rep", i, ".csv"))
        tbl
      })
    
    df_all <- lapply(df_all,function(i){
        i[,c("chrom","start_position1","end_position1","Nvalid_cov_DMSO","Nvalid_cov_STORM","percent_modified_DMSO","percent_modified_STORM")]
    })
    names(df_all) <- sapply(1:length(mod_sites),function(i)paste0("Rep",i))
    
    
    ###################### PLOT DMSO vs STORM ####################################
    
    lapply(names(df_all), function(i) {
      # Prepara i dati
      plot_data <- df_all[[i]] %>%
        mutate(
          x_plot = ifelse(is.na(percent_modified_DMSO), 0, percent_modified_DMSO),
          y_plot = ifelse(is.na(percent_modified_STORM), 0, percent_modified_STORM),
          status = case_when(
            is.na(percent_modified_DMSO) ~ "DMSO_NA",
            is.na(percent_modified_STORM) ~ "STORM_NA",
            TRUE ~ "common"
          )
        )
      
      # Jitter per NA
      set.seed(123)
      plot_data <- plot_data %>%
        mutate(
          x_plot = ifelse(status == "DMSO_NA", jitter(x_plot, amount = 0.15), x_plot),
          y_plot = ifelse(status == "STORM_NA", jitter(y_plot, amount = 0.15), y_plot)
        )
      
      max_val <- max(c(plot_data$x_plot, plot_data$y_plot), na.rm = TRUE)
      
      # Separa i dati per categoria e conta
      data_common <- plot_data %>% filter(status == "common")
      data_DMSO_NA <- plot_data %>% filter(status == "DMSO_NA")
      data_STORM_NA <- plot_data %>% filter(status == "STORM_NA")
      
      n_common <- nrow(data_common)
      n_DMSO_NA <- nrow(data_DMSO_NA)
      n_STORM_NA <- nrow(data_STORM_NA)
      
      # Crea etichette con conteggi
      legend_labels <- c(
        sprintf("Common (n=%s)", format(n_common, big.mark=",")),
        sprintf("DMSO NA (n=%s)", format(n_DMSO_NA, big.mark=",")),
        sprintf("STORM NA (n=%s)", format(n_STORM_NA, big.mark=","))
      )
    
      print(legend_labels)
      
      # Plot con geom_point invece di scattermore
      p <- ggplot(plot_data, aes(x = x_plot, y = y_plot)) +
        rasterize(geom_point(data = data_common, 
                   color = palette[k], 
                   size = 0.5),dpi=300) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                    linewidth = 1, color = "black") +
        coord_fixed(ratio = 1, xlim = c(0, max_val), ylim = c(0, max_val)) +
        scale_x_continuous(breaks = sort(unique(c(seq(0, max_val, by = 5), 5, 10, 15, 20)))) +
      scale_y_continuous(breaks = sort(unique(c(seq(0, max_val, by = 5), 5, 10, 15, 20)))) +
        labs(x = "DMSO", y = "STORM", title = "Stoichiometry DMSO vs STORM") +
        theme_classic() +
        # annotate("point", x = max_val*0.08, y = max_val*c(0.95, 0.91, 0.87),
        #          color = c("blue", "red", "black"), size = 3) +
        annotate("text", x = max_val*0.11, y = max_val*c(0.95, 0.91, 0.87),
                 label = legend_labels,
                 hjust = 0, size = 3.5)
      
      # Aggiungi violin plots marginali
      p_with_marginals <- ggMarginal(p, 
                                      type = "violin",
                                      fill = palette[k],
                                      alpha = 0.5,
                                      size = 4)
      
      ggsave(paste0(path,"DMSOvsSTORMviolin_rast_", i, ".pdf"), 
             plot = p_with_marginals, 
             width = 8, height = 8)
    
    
      ggsave(paste0(path,"DMSOvsSTORMviolin_rast_", i, ".svg"), 
             plot = p_with_marginals, 
             width = 8, height = 8)
    
    })
}, mc.cores = n_cores)






