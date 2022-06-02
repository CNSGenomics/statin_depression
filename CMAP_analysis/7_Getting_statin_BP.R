library(ggplot2)
library(dplyr)
library(data.table)
library(ggrepel)
library(cowplot)
library(egg)

########################################################################
### Defining directorys and lists
########################################################################
statin_gProfiler_dir <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/5_gProfiler2_analysis/gProfiler2_results_gSCS_files_all_genes/'

Statins_HA1E_list_no_prava <- c("atorvastatin", "fluvastatin", "lovastatin", "mevastatin", "rosuvastatin", "simvastatin")
gene_dir_vector <- c('up', 'down')

########################################################################
### There are too many terms for the down-regulated group, so i'm gonna group the GO terms into parent terms using the revigo web server
########################################################################
All_merge_Df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(All_merge_Df) <- c('term_id', 'term_size', 'term_name', 'minuslog10P', 'statin', 'direction')

for (in_statin in Statins_HA1E_list_no_prava) {
  for (in_gene_direction in gene_dir_vector){
    print(paste(c(in_statin, in_gene_direction), collapse = '+'))
    curr_statn_gProfiler_file <- paste(c(statin_gProfiler_dir, in_statin, '_HA1E_DEG_thresh1_all_genes_', in_gene_direction, '_gProfiler_gSCS.txt'), collapse ='')
    curr_df <- fread(curr_statn_gProfiler_file)
    curr_df_BP <- curr_df[which(curr_df$source == 'GO:BP'),]
    curr_df_BP$minuslog10P <- -log10(curr_df_BP$p_value)
    curr_df_BP_fil <- curr_df_BP[,c('term_id', 'term_size', 'term_name', 'minuslog10P')]
    curr_df_BP_fil$statin <- in_statin
    curr_df_BP_fil$direction <- in_gene_direction
    
    All_merge_Df <- rbind(All_merge_Df, curr_df_BP_fil)
  }
}

GO_BP_term <- All_merge_Df[,c('term_id')]
unique_GO_BP_term <- GO_BP_term[!duplicated(GO_BP_term)]

unique_GO_BP_term_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/5_gProfiler2_analysis/gProfiler2_results_gSCS_files_all_genes/Statins_HA1E_GO_BP_sig_in_atleast1_up_down.csv'
write.csv(unique_GO_BP_term, unique_GO_BP_term_file, quote = FALSE, row.names = FALSE)
