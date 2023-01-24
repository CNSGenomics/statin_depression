### This script is used to perform pathway enrichment analysis (gProfiler2) using differentially expressed genes from statin signatures. The z-score threshold used to define differentially expressed genes is arbitrary. In our analysis we performed analysis using z-score thresholds of 1, 1.5 and 2.

library(gprofiler2)
library(dplyr)
library(data.table)
library(stringr)

statin_gProfiler_dir <- './5_Statin_DEG_gProfiler/'

############################################################################
### Setting the bg genes, which are all the genes profiled by LINCS
############################################################################
### Getting background genes
lincs_level5_gene_file <- './Data/LINCS/GSE92742_Broad_LINCS_gene_info.txt'
lincs_level5_gene_df <- read.csv(lincs_level5_gene_file, sep = '\t', header = TRUE)
bg_gene_list <- lincs_level5_gene_df$pr_gene_id

### The function for running gProfiler2
Running_gProfiler2 <- function(in_query){
  gost_result <- gost(query = in_query, custom_bg = bg_gene_list, correction_method = 'gSCS', domain_scope = 'custom', organism = 'hsapiens')
  gost_results_df <- gost_result$result
  gost_results_df <- data.frame(lapply(gost_results_df, as.character), stringsAsFactors=FALSE)
  return(gost_results_df)
}

############################################################################
### Running gProfiler2 analysis using Ensembl 104
############################################################################
statin_pathway <- function(statin_z_file, in_sub_dir, in_cell, in_thresh) {
  ### Set the version of the gProfiler data to Ensembl 104 to keep it consistent
  set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e104_eg51_p15")
  
  ### Getting the DEGs and running gProfilers
  summary_line <- paste(c('Statin', 'Cell', 'Threshold', 'Up', "Down"), collapse = '\t')
  summary_list <- c(summary_line)
  
  ### Getting DEGs
  statin_z_df <- fread(statin_z_file)
  
  for (curr_statin in (colnames(statin_z_df)[colnames(statin_z_df) != 'V1'])) {
    if (!curr_statin %like% 'pravastatin') { # Removing pravastatin due to its inconsistency with the other statins
      curr_statin_str <- str_split(curr_statin, '_')[[1]][1]
      print(curr_statin_str)

      ### Setting the z-score column to numeric values
      statin_z_df[[curr_statin]] <- as.numeric(statin_z_df[[curr_statin]])
      
      ### Up-regulated genes
      up_df <- statin_z_df[statin_z_df[[curr_statin]] > in_thresh, ]
      up_genes <- up_df[['V1']]
      count_up <- length(up_genes)
  
      out_gProfiler2_file <- paste(c(statin_gProfiler_dir, in_sub_dir, curr_statin_str, '_', in_cell, '_thresh', in_thresh, '_DEGs_gSCS_gProfiler2_up.txt'), collapse = '')
      write.table(Running_gProfiler2(up_genes), out_gProfiler2_file, sep = '\t', row.names = FALSE)
     
      ### Down-regulated genes
      down_df <- statin_z_df[statin_z_df[[curr_statin]] < ((in_thresh)*(-1)), ]
      down_genes <- down_df[['V1']]
      count_down <- length(down_genes)
      
      out_gProfiler2_file <- paste(c(statin_gProfiler_dir, in_sub_dir, curr_statin_str, '_', in_cell, '_thresh', in_thresh, '_DEGs_gSCS_gProfiler2_down.txt'), collapse = '')
      write.table(Running_gProfiler2(down_genes), out_gProfiler2_file, sep = '\t', row.names = FALSE) 
      
      ### Summarising the gene counts
      summary_line <- paste(c(curr_statin, in_cell, in_thresh, count_up, count_down), collapse = '\t')
      summary_list <- c(summary_list, summary_line)
    }
  }

  summary_file <- paste(c(statin_gProfiler_dir, in_sub_dir, 'Statins_', in_cell, '_thresh', in_thresh, '_DEGs_count_summary.txt'), collapse = '')
  fwrite(as.list(summary_list), summary_file, sep = '\n', quote = FALSE)
}

### HA1E
statin_pathway('./2_z_score_files/Statin_HA1E_10uM_24h_all_genes.txt', 'HA1E_statins/', 'HA1E', 1) 
statin_pathway('./2_z_score_files/Statin_HA1E_10uM_24h_all_genes.txt', 'HA1E_statins/', 'HA1E', 1.5)
statin_pathway('./2_z_score_files/Statin_HA1E_10uM_24h_all_genes.txt', 'HA1E_statins/', 'HA1E', 2) 

### NPC
statin_pathway('./2_z_score_files/Statin_NPC_10uM_24h_all_genes.txt', 'NPC_statins/', 'NPC', 1) 
statin_pathway('./2_z_score_files/Statin_NPC_10uM_24h_all_genes.txt', 'NPC_statins/', 'NPC', 1.5) 
statin_pathway('./2_z_score_files/Statin_NPC_10uM_24h_all_genes.txt', 'NPC_statins/', 'NPC', 2) 
