### This script is used to generate gmt signatures, which will then be submitted to the CLUE platform to compute connectivity profiles. Each gmt signature contains n most upregulated landmark genes, and n most downregulated landmark genes. The selection of n is arbitrary, and the default setting used by CMap for internal queries is n = 50

library(dplyr)
library(stringr)
library(data.table)

gmt_dir <- './3_gmt_files/'

###########################################################################
### Defining functions to generate the top regulated and down regulated genes
###########################################################################
Getting_top_upreg <- function(in_fil_df, number_of_top_genes) {
  Upreg_fil_df <- in_fil_df[in_fil_df$Zval > 0, ]
  print(paste(c('Total number of Upreg gene', length(row.names(Upreg_fil_df))), collapse = ' '))
  Upreg_fil_df_sorted <- Upreg_fil_df[order(Upreg_fil_df$Zval, decreasing = TRUE), ]
  Top_upreg_genes <- (Upreg_fil_df_sorted$genes)[1:number_of_top_genes] 
  return(c(Top_upreg_genes))
}
Getting_top_downreg <- function(in_fil_df, number_of_top_genes) {
  Downreg_fil_df <- in_fil_df[in_fil_df$Zval < 0, ]
  print(paste(c('Total number of Downreg gene', length(row.names(Downreg_fil_df))), collapse = ' '))
  Downreg_fil_df_sorted <- Downreg_fil_df[order(Downreg_fil_df$Zval, decreasing = FALSE), ]
  Top_downreg_genes <- (Downreg_fil_df_sorted$genes)[1:number_of_top_genes]
  return(c(Top_downreg_genes))
}  

###########################################################################
### Generating query signatures
###########################################################################
generating_gmt <- function(in_lm_file, in_top_gene_count, in_subdir, outfile){
  ### Import dataframe
  in_lm_df <- fread(in_lm_file, data.table = FALSE)
  
  ### Generating a gmt signature for each of the drugs
  uplist <- list()
  downlist <- list()
  for (in_sig_name in (colnames(in_lm_df)[2:length(colnames(in_lm_df))]) ) {
    
    ### Formatting the sig names
    curr_cp <- str_split(in_sig_name, '_')[[1]][1]
    curr_cell <- str_split(in_sig_name, '_')[[1]][3]
    curr_sig_name <- paste(c(str_to_title(curr_cp), '_', curr_cell), collapse = '')
    print(head(curr_sig_name))
    
    ### Getting the gene and sig columns
    curr_df <- in_lm_df[, c('V1', in_sig_name)]
    colnames(curr_df) <- c('genes', 'Zval')
    
    curr_query_up<-Getting_top_upreg(curr_df, in_top_gene_count)
    curr_query_down<-Getting_top_downreg(curr_df, in_top_gene_count)
    
    ### Saving the gmt file for each signature
    curr_up_gmt_list <- list('head' = curr_sig_name, 'desc' = in_sig_name, 'len' = in_top_gene_count, 'entry' = curr_query_up)
    curr_up_gmt_list_out <- list(curr_sig_name = curr_up_gmt_list)
    up_gmt_file <- paste(c(gmt_dir, in_subdir, curr_sig_name, '_', in_top_gene_count, '_lm_genes_24h_10uM_up_query.gmt'), collapse = "")
    write_gmt(curr_up_gmt_list_out, up_gmt_file)

    curr_down_gmt_list <- list('head' = curr_sig_name, 'desc' = in_sig_name, 'len' = in_top_gene_count, 'entry' = curr_query_down)
    curr_down_gmt_list_out <- list(curr_sig_name = curr_down_gmt_list)
    down_gmt_file <- paste(c(gmt_dir, in_subdir, curr_sig_name, '_', in_top_gene_count, '_lm_genes_24h_10uM_down_query.gmt'), collapse = "")
    write_gmt(curr_down_gmt_list_out, down_gmt_file)
    
    ### Appending the individual gmt to a master list 
    up_gmt <- parse_gmt(up_gmt_file)
    down_gmt <- parse_gmt(down_gmt_file)
        
    uplist <- c(uplist, curr_sig_name = up_gmt)
    downlist <- c(downlist, curr_sig_name = down_gmt)
  }
  
  ### Collating the individual gmt lists into a gmt file
  all_up_gmt_file <- paste(c(gmt_dir, in_subdir, outfile, '_up.gmt'), collapse = '')
  write_gmt(uplist, all_up_gmt_file)
  print(length(uplist))
  
  all_down_gmt_file <- paste(c(gmt_dir, in_subdir, outfile, '_down.gmt'), collapse = '')
  write_gmt(downlist, all_down_gmt_file)
  print(length(downlist))
}


generating_gmt('./2_z_score_files/Statin_HA1E_10uM_24h_lm_genes.txt', 50, 'Statin_HA1E_10uM_24h_50/', 'Statin_HA1E_10uM_24h_50')
generating_gmt('./2_z_score_files/Statin_HA1E_10uM_24h_lm_genes.txt', 100, 'Statin_HA1E_10uM_24h_100/', 'Statin_HA1E_10uM_24h_100')
generating_gmt('./2_z_score_files/Statin_HA1E_10uM_24h_lm_genes.txt', 150, 'Statin_HA1E_10uM_24h_150/', 'Statin_HA1E_10uM_24h_150')
generating_gmt('./2_z_score_files/Statin_NPC_10uM_24h_lm_genes.txt', 50, 'Statin_NPC_10uM_24h_50/', 'Statin_NPC_10uM_24h_50')

###########################################################################
### Submitted to CLUE query, compared against v1.0 without the fastgutc tool
###########################################################################
### Gene expression (L1000)
### Touchstone
### Batch query
### 1.0
