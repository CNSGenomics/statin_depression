library(gprofiler2)
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(ggrepel)

Common_DEG_dir <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/10_Common_DEGs_statin_antidepressant/'
Ind_dir <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/10_Common_DEGs_statin_antidepressant/Individual_statin_vs_antidepressant/'
Ind_gProfiler2_dir <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/10_Common_DEGs_statin_antidepressant/Individual_statin_vs_antidepressant/Individual_statin_vs_antidepressant_gProfiler2/'
Ind_scatter_dir <-'/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/10_Common_DEGs_statin_antidepressant/Individual_statin_vs_antidepressant/Individual_statin_vs_antidepressant_scatter/'

Statins_HA1E_list_no_prava <- c("atorvastatin", "fluvastatin", "lovastatin", "mevastatin", "rosuvastatin", "simvastatin")
antidepressant_drugs <- c('desipramine', 'trimipramine', 'nortriptyline', 'paroxetine', 'sertraline')
############################################################################
### Getting background genes
############################################################################
lincs_level5_gene_file <- '/Users/uqjjian3/Data/LINCS/GSE92742_Broad_LINCS_gene_info.txt' # Retreived from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
lincs_level5_gene_df <- read.csv(lincs_level5_gene_file, sep = '\t', header = TRUE)
All_gene_list <- lincs_level5_gene_df$pr_gene_id

############################################################################
### The function for running gProfiler2
############################################################################
Running_gProfiler2 <- function(in_query){
  gost_result <- gost(query = in_query, custom_bg = All_gene_list, correction_method = 'gSCS', domain_scope = 'custom', organism = 'hsapiens')
  gost_results_df <- gost_result$result
  gost_results_df = data.frame(lapply(gost_results_df, as.character), stringsAsFactors=FALSE)
  return(gost_results_df)
}

############################################################################
### Z-score plots of the statin and antidepressant DEG genes
############################################################################
### Getting the antidepressant z-scores
Sig_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/7_Analysis_psyc_drugs_z_scores/z_score_files/HA1E_antidepressant_High_Tau_z_scores.txt'
drug_z_df <- read.csv(Sig_file, sep = '\t', header = TRUE, row.names = 1)
old_colnames <- colnames(drug_z_df)
antidepressant_vector <- c()
for (old_name in old_colnames){
  new_name <- strsplit(old_name, '_')[[1]][1]
  antidepressant_vector <- c(antidepressant_vector, new_name)
}
colnames(drug_z_df) <- antidepressant_vector

### Getting the statin z-scores
statin_sig_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/1_GSE92742_sig/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_Statins_HA1E_all_genes.txt'
statin_zscores_df <- read.csv(statin_sig_file, sep = '\t', header = TRUE, row.names = 1)
statin_zscores_df <- subset(statin_zscores_df, select = -c(pravastatin_HA1E) )
old_colnames <- colnames(statin_zscores_df)
statin_vector <- c()
for (old_name in old_colnames){
  new_name <- strsplit(old_name, '_')[[1]][1]
  statin_vector <- c(statin_vector, new_name)
}
colnames(statin_zscores_df) <- statin_vector
### Combining statin and antidepressant df
statin_antidepressant_z_df <- cbind(statin_zscores_df, drug_z_df)

### Getting the common DEGs
summary_line <- paste(c('Statin', 'Antidepressant', 'Both_up', "Both_down", "Same_direction", "Statin_up_drug_down", "Statin_down_drug_up", "Opposite_direction"), collapse = '\t')
summary_list <- c(summary_line)
Getting_common_DEGs_z <- function(in_statin, in_drug){
  print(paste(c(in_statin, in_drug), collapse = ' '))
  curr_df <- statin_antidepressant_z_df[,c(in_statin,in_drug)]
  
  common_up_df <- curr_df[curr_df[[in_statin]]>1 & curr_df[[in_drug]]>1, ]
  common_up_genes <- row.names(common_up_df)

  common_down_df <- curr_df[curr_df[[in_statin]] < -1 & curr_df[[in_drug]] < -1, ]
  common_down_genes <- row.names(common_down_df)

  same_direction_genes <- unique(c(common_up_genes, common_down_genes))
  out_gProfiler2_file <- paste(c(Ind_gProfiler2_dir, in_statin, '_', in_drug, '_thresh1_DEGs_gSCS_gProfiler2_same_direction.txt'), collapse = '')
  write.table(Running_gProfiler2(same_direction_genes), out_gProfiler2_file, sep = '\t', row.names = FALSE)

  statin_up_drug_down_df <- curr_df[curr_df[[in_statin]]>1 & curr_df[[in_drug]]< -1, ]
  statin_up_drug_down_genes <- row.names(statin_up_drug_down_df)

  statin_down_drug_up_df <- curr_df[curr_df[[in_statin]] < -1 & curr_df[[in_drug]]>1, ]
  statin_down_drug_up_genes <- row.names(statin_down_drug_up_df)

  opposite_direction_genes <- unique(c(statin_up_drug_down_genes, statin_down_drug_up_genes))
  out_gProfiler2_file <- paste(c(Ind_gProfiler2_dir, in_statin, '_', in_drug, '_thresh1_DEGs_gSCS_gProfiler2_opposite_direction.txt'), collapse = '')
  write.table(Running_gProfiler2(opposite_direction_genes), out_gProfiler2_file, sep = '\t', row.names = FALSE)
  
  summary_line <- paste(c(in_statin, in_drug, length(common_up_genes), length(common_down_genes), length(same_direction_genes), length(statin_up_drug_down_genes), length(statin_down_drug_up_genes), length(opposite_direction_genes)), collapse = '\t')
  summary_list <<- c(summary_list, summary_line)
}

combination_df <- expand.grid(Statins_HA1E_list_no_prava, antidepressant_drugs, stringsAsFactors = FALSE)

i = 1
while (i <= length(row.names(combination_df))) {
  curr_combo <- as.character(combination_df[i,])
  Getting_common_DEGs_z(curr_combo[1], curr_combo[2])
  i = i + 1
}

summary_file <- paste(c(Ind_dir, 'Individual_statins_vs_antidepressants_DEGs_count_summary.txt'), collapse = '')
fwrite(as.list(summary_list), summary_file, sep = '\n', quote = FALSE)


### Drawing scatter of DEGs
Drawing_scatter_DEGs_z <- function(in_statin, in_drug){
  print(paste(c(in_statin, in_drug), collapse = ' '))
  curr_df <- statin_antidepressant_z_df[,c(in_statin,in_drug)]
  
  common_up_df <- curr_df[curr_df[[in_statin]]>1 & curr_df[[in_drug]]>1, ]
  common_down_df <- curr_df[curr_df[[in_statin]] < -1 & curr_df[[in_drug]] < -1, ]
  statin_up_drug_down_df <- curr_df[curr_df[[in_statin]]>1 & curr_df[[in_drug]]< -1, ]
  statin_down_drug_up_df <- curr_df[curr_df[[in_statin]] < -1 & curr_df[[in_drug]]>1, ]
  
  all_DEGs_DF <- rbind(common_up_df, common_down_df)
  all_DEGs_DF <- rbind(all_DEGs_DF, statin_up_drug_down_df)
  all_DEGs_DF <- rbind(all_DEGs_DF, statin_down_drug_up_df)
  
  all_DEGs_count <- length(row.names(all_DEGs_DF))
  
  out_pdf <- paste(c(Ind_scatter_dir, in_statin, '_', in_drug, '_thresh1_common_DEGs_z_scores.pdf'),collapse = '')
  pdf(out_pdf, width = 4.5, height = 5)
  par(mgp=c(1.2,0.5,0))
  plot(all_DEGs_DF[[in_statin]], all_DEGs_DF[[in_drug]], xlab = in_statin, ylab=in_drug, xlim = c(-10,10), ylim = c(-10,10), col = 'grey')
  abline(h = 0, col = "black")
  abline(v = 0, col = "black")
  abline(h = c(-1, 1), col = "red", lty=2)
  abline(v = c(-1, 1), col = "red", lty=2)
  text(x = -10, y = 10, paste(c(in_statin, '\n', in_drug), collapse = ''), adj = c(0,1))
  text(x = 10, y = 10, paste(c('n = ', all_DEGs_count), collapse = ''), adj = c(1,1))
  # legend('topleft', 'text', lwd = 0, xjust = 0, yjust = 0)
  # legend('topright', 'textq', lwd = 0, xjust = 0, yjust = 0)
  dev.off()
}

i = 1
while (i <= length(row.names(combination_df))) {
  curr_combo <- as.character(combination_df[i,])
  Drawing_scatter_DEGs_z(curr_combo[1], curr_combo[2])
  i = i + 1
}

