library(cmapR)
library(dplyr)
library(stringr)
library(data.table)

sig_dir <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/1_GSE92742_sig/'
all_z_file <- paste(c(sig_dir, 'GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_Statins_HA1E_all_genes.txt'), collapse = '')

###########################################################################
### 1. Preparing dataframes
###########################################################################
Preparing_df <- function(in_file){
  z_scores_df <- read.csv(in_file, sep = '\t', header = TRUE, row.names = 1)
  print(dim(z_scores_df))
  return(z_scores_df)
}
Preparing_df_no_prava <- function(in_file){
  z_scores_df <- read.csv(in_file, sep = '\t', header = TRUE, row.names = 1)
  z_scores_df <- z_scores_df[ , -which(names(z_scores_df) %in% c('pravastatin_HA1E'))]
  print(dim(z_scores_df))
  return(z_scores_df)
}

###########################################################################
### 2. Plotting pairwise z scores
###########################################################################
### A function to plot the lm regression line in the subplots
lm_reg <- function(x,y,...){
  points(x,y)
  lm_model <- lm(y~x)
  abline(lm_model, col = 'red')
  abline (h = 0, col = "grey")
  abline (v = 0, col = "grey")
  beta_Value <- coef(lm_model)[[2]]
  p_beta <- summary(lm_model)$coefficient[2,4]
  R2_adj <- summary(lm_model)$adj.r.squared
  text(x = -5, y = 7, labels = paste(c('R2=',round(R2_adj,2),'\n', 'b=', round(beta_Value,2), '\n', 'p=', format(p_beta, digits=3)), collapse= '') ,col='red')
}
drawing_pairwise_LR <- function(in_cell_line, in_df, gene_group, prava_state){
  ### Remove the cell line names and reordering the columns
  to_remove <- paste(c('_', in_cell_line), collapse = '')
  cell_line_z_df <- in_df
  colnames(cell_line_z_df) <- colnames(cell_line_z_df) %>% str_replace(to_remove, "")
  
  #Plotting the pair-wise z scatter plots
  out_pdf = paste(c(Pairwise_dir, '/Statins_', in_cell_line, '_', gene_group, '_z_scores_pairwise_LR_24h_', prava_state, '.pdf'), collapse = '')
  
  plot_title <- paste(c(in_cell_line, gene_group, '24h', prava_state), collapse = '_')
  pdf(out_pdf, width = 10, height = 10)
  pairs(cell_line_z_df, lower.panel=NULL, panel = 
          lm_reg, xlim= c(-10,10) , ylim= c(-10, 10), main = plot_title, cex.axis=0.8, 
        )
  dev.off()
}  

###########################################################################
### 3. Running commands
###########################################################################
Pairwise_dir <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/3_Z_scores_analysis/Statins_HA1E_pairwise_z/'
Selected_cell_id_TS <- c("HA1E")
for (cell_name in Selected_cell_id_TS) {
  print(cell_name)
  print('With pravastatin')
  print('All genes')
  z_df_all <- Preparing_df(all_z_file)
  drawing_pairwise_LR(cell_name, z_df_all, "all_genes",'')

  print('Without pravastatin')
  print('All genes')
  z_df_all <- Preparing_df_no_prava(all_z_file)
  drawing_pairwise_LR(cell_name, z_df_all, "all_genes", 'no_prava')
}




