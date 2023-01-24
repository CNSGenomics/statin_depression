### This script is used to perform pathway enrichment analysis for genes perturbed in the same and opposite directions by statins and antidepressants

library(gprofiler2)
library(dplyr)
library(data.table)
library(stringr)

Statin_antidep_gprofiler2_dir <- './8_Statins_antidepressant_DEG_gProfiler2/' 


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
### A function for performing gProfiler2 for DEGs shared between a pair of statin and antidepressant
Getting_common_DEGs_z <- function(in_statin, in_drug, in_merged_df, in_sub_dir, in_cell){
  ### Set the version of the gProfiler data to Ensembl 104 to keep it consistent with the other analysis
  set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e104_eg51_p15")
  
  print(paste(c(in_statin, in_drug), collapse = ' '))
  curr_df <- in_merged_df[,c(in_statin,in_drug)]

  ### genes perturbed in the same direction
  common_up_df <- curr_df[curr_df[[in_statin]]> 1 & curr_df[[in_drug]]> 1, ]
  common_up_genes <- row.names(common_up_df)
  common_down_df <- curr_df[curr_df[[in_statin]] < -1 & curr_df[[in_drug]] < -1, ]
  common_down_genes <- row.names(common_down_df)

  same_direction_genes <- unique(c(common_up_genes, common_down_genes))
  out_gProfiler2_file <- paste(c(Statin_antidep_gprofiler2_dir, in_sub_dir, in_statin, '_', in_drug , '_', in_cell, '_thresh1_DEGs_gSCS_gProfiler2_same_direction.txt'), collapse = '')
  write.table(Running_gProfiler2(same_direction_genes), out_gProfiler2_file, sep = '\t', row.names = FALSE)

  ### genes perturbed in the opposite direction
  statin_up_drug_down_df <- curr_df[curr_df[[in_statin]]> 1 & curr_df[[in_drug]]< -1, ]
  statin_up_drug_down_genes <- row.names(statin_up_drug_down_df)
  statin_down_drug_up_df <- curr_df[curr_df[[in_statin]] < -1 & curr_df[[in_drug]] > 1, ]
  statin_down_drug_up_genes <- row.names(statin_down_drug_up_df)

  opposite_direction_genes <- unique(c(statin_up_drug_down_genes, statin_down_drug_up_genes))
  out_gProfiler2_file <- paste(c(Statin_antidep_gprofiler2_dir, in_sub_dir, in_statin, '_', in_drug , '_', in_cell, '_thresh1_DEGs_gSCS_gProfiler2_opposite_direction.txt'), collapse = '')
  write.table(Running_gProfiler2(opposite_direction_genes), out_gProfiler2_file, sep = '\t', row.names = FALSE)

  ## Summarising the gene counts
  summary_line <- paste(c(in_statin, in_drug, length(same_direction_genes), length(opposite_direction_genes)), collapse = '\t')
  return(summary_line)
}

### A function for performing gProfiler2 for DEGs shared between statins and antidepressants from a cell type
comparing_statin_antidepressant <- function(in_statin_z_file, in_antidepressant_z_file, in_sub_dir, in_cell) {
  ### Importing the statin and antidepressant z-scores dataframe
  in_statin_z_df <- fread(in_statin_z_file, data.table = FALSE)
  in_antidepressant_z_df <- fread(in_antidepressant_z_file, data.table = FALSE)
  
  ### Removing pravastatin from analysis as it displayed much inconsistency with the other statins
  in_statin_z_df <- in_statin_z_df[ , c(str_subset(colnames(in_statin_z_df), 'pravastatin', negate = TRUE))]
  
  ### Formatting the colnames
  for (curr_col in colnames(in_statin_z_df)) {
    new_col <- str_split(curr_col, '_')[[1]][1]
    new_col <- str_to_title(new_col)
    colnames(in_statin_z_df) <- str_replace(colnames(in_statin_z_df), curr_col, new_col)
  } 
  for (curr_col in colnames(in_antidepressant_z_df)) {
    new_col <- str_split(curr_col, '_')[[1]][1]
    new_col <- str_to_title(new_col)
    colnames(in_antidepressant_z_df) <- str_replace(colnames(in_antidepressant_z_df), curr_col, new_col)
  } 
  
  ### Ordering the columns to the alphabetical order
  in_statin_z_df <- in_statin_z_df[, order(colnames(in_statin_z_df))]
  in_antidepressant_z_df <- in_antidepressant_z_df[, order(colnames(in_antidepressant_z_df))]
  
  ### merging the dataframe and setting the gene names to row.names
  merged_df <- merge(in_statin_z_df, in_antidepressant_z_df, by = 'V1')
  row.names(merged_df) <- merged_df$V1
  merged_df$V1 <- NULL
  
  ### Getting the common DEGs and running gProfilers
  summary_line <- paste(c('Statin', 'Antidepressant', "Same_direction", "Opposite_direction"), collapse = '\t')
  summary_list <- c(summary_line)

  combination_df <- expand.grid(colnames(in_statin_z_df[, colnames(in_statin_z_df) != 'V1']), colnames(in_antidepressant_z_df[, colnames(in_antidepressant_z_df) != 'V1']), stringsAsFactors = FALSE)

  i = 1
  while (i <= length(row.names(combination_df))) {
    curr_combo <- as.character(combination_df[i,])
    curr_summary_line <- Getting_common_DEGs_z(curr_combo[1], curr_combo[2], merged_df, in_sub_dir, in_cell)

    summary_list <- c(summary_list, curr_summary_line)
    i = i + 1
  }

  summary_file <- paste(c(Statin_antidep_gprofiler2_dir, in_sub_dir, 'Individual_statins_vs_antidepressants_', in_cell, '_thresh1_DEGs_count_summary.txt'), collapse = '')
  fwrite(as.list(summary_list), summary_file, sep = '\n', quote = FALSE)

}

comparing_statin_antidepressant('./2_z_score_files/Statin_HA1E_10uM_24h_all_genes.txt', './2_z_score_files/Antidepresssants_HA1E_10uM_24h_all_genes.txt', 'HA1E_statin_antidepressant/', 'HA1E')
comparing_statin_antidepressant('./2_z_score_files/Statin_NPC_10uM_24h_all_genes.txt', './2_z_score_files/Antidepressants_NPC_10uM_24h_all_genes.txt', 'NPC_statin_antidepressant/', 'NPC')
