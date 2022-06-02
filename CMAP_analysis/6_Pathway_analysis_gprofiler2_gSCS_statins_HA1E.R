library(gprofiler2)
library(dplyr)
library(data.table)

DEG_dir <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/4_HA1E_24_DEGs_thresh1/'
gProfiler_dir <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/5_gProfiler2_analysis/'

###########################################################################
### 1. Getting the background genes
###########################################################################
lincs_level5_gene_file <- '/Users/uqjjian3/Data/LINCS/GSE92742_Broad_LINCS_gene_info.txt'

lincs_level5_gene_df <- read.csv(lincs_level5_gene_file, sep = '\t', header = TRUE)
All_gene_list <- lincs_level5_gene_df$pr_gene_id

###########################################################################
### 2. Performing pathway enrichment analysis
###########################################################################

Getting_ind_statin_gene_list <- function(in_cell_line, in_statin, gene_group, condition){
  statin_cell_DEG_file <- paste(c(DEG_dir, gene_group, '/', in_statin, '_', in_cell_line, '_24h_', gene_group, '_', condition, '.txt'), collapse = '')
  statin_cell_DEG_list_Df <- read.csv(statin_cell_DEG_file, sep = '\t', header = FALSE)
  colnames(statin_cell_DEG_list_Df) = c('ID')
  statin_cell_DEG_list <- statin_cell_DEG_list_Df$ID
  return(statin_cell_DEG_list)
}


Running_gprofiler2 <- function(in_cell_line, in_statin, gene_group, condition, in_bg_list, in_gene_space, in_thresh, out_dir){
  print(paste(c(in_cell_line, in_statin, gene_group, condition, in_thresh), collapse = ' '))
  ### Getting the query gene list
  in_query_list <- Getting_ind_statin_gene_list(in_cell_line, in_statin, gene_group, condition)

  Number_of_query_genes <-  length(in_query_list)
  print(paste(c('Number of query genes',Number_of_query_genes), collapse = ' '))
  
  ### Running gProfiler2
  gost_results <- gost(query = in_query_list, custom_bg = in_bg_list, correction_method = 'gSCS', domain_scope = 'custom', organism = 'hsapiens')
  gost_results_sig_results <- gost_results$result
  gost_results_sig_meta <- gost_results$meta
  
  actual_gene_count <- length(gost_results_sig_meta$genes_metadata$query$query_1$mapping)
  ambiguous_genes_list <- paste(c(names(gost_results_sig_meta$genes_metadata$ambiguous)), collapse = '|')
  failed_genes_list <- paste(c(gost_results_sig_meta$genes_metadata$failed), collapse = '|')
  dup_genes_list <- paste(c(names(gost_results_sig_meta$genes_metadata$duplicates)), collapse = '|')
  
  gost_results_sig_results = data.frame(lapply(gost_results_sig_results, as.character), stringsAsFactors=FALSE)
  results_out_file <- paste(c(gProfiler_dir, out_dir, in_statin, '_', in_cell_line, '_', gene_group, '_', condition, '_gProfiler_gSCS.txt'), collapse = '')
  write.table(gost_results_sig_results,results_out_file, sep = '\t', row.names = FALSE)

  if ('GO:BP' %in% gost_results_sig_results$source) {
    gost_results_fil <- filter(gost_results_sig_results, gost_results_sig_results$source == 'GO:BP')
    gost_results_GP_list <- paste(gost_results_fil$term_name, collapse = '|')
    summary_line <- c(in_statin, in_cell_line, condition, in_thresh, Number_of_query_genes,actual_gene_count, ambiguous_genes_list, failed_genes_list, dup_genes_list, length(rownames(gost_results_fil)), gost_results_GP_list)
  } else {
    summary_line <- c(in_statin, in_cell_line, condition, in_thresh, Number_of_query_genes,actual_gene_count, ambiguous_genes_list, failed_genes_list, dup_genes_list, '0')
  }
  return(paste(summary_line,collapse = '\t'))
}

Statins_HA1E_list <- c("atorvastatin", "fluvastatin", "lovastatin", "mevastatin", "pravastatin", "rosuvastatin", "simvastatin")
DE_dir_list <- c('up', 'down')


############################################################################
### Run for all genes
summary_header_line <- paste(c('statin_name', 'cell_line', 'DE_direction', 'z_thresh', 'Number_of_input_query_genes', 'Number_of_actual_query_genes', 'Ambiguous_genes_list', 'Failed_genes_list', 'Dup_genes_list', 'Number_of_sig_GO_BP', 'List_of_sig_GP_BP'), collapse = '\t')
summary_line_list <- c()
summary_line_list <- c(summary_line_list,summary_header_line)
for (query_statin in Statins_HA1E_list){
  Cell = 'HA1E'
  name_to_check = paste(c(query_statin, Cell), collapse = '_')
  for (DE_dir in DE_dir_list) {
    summary_line_1 <- Running_gprofiler2(Cell, query_statin, 'DEG_thresh1_all_genes', DE_dir, All_gene_list, 'all_genes', '1', 'gProfiler2_results_gSCS_files_all_genes/')
    summary_line_list <- c(summary_line_list, summary_line_1)
  }
}

results_out_file <- paste(c(gProfiler_dir, 'Statins_HA1E_all_genes_DEG_gProfiler_gSCS_results_summary.txt'), collapse = '')
fwrite(as.list(summary_line_list), results_out_file, sep = '\n')



