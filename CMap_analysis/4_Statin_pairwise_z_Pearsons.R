### This script is used to calculate the Peasons' correlation between statin gene expression signatures

library(stringr)
library(data.table)

Cor_plot_dir <- './4_statin_z_score_Pearsons/'

############################################################
### A function for drawing the scatter plots with the Pearson's correlation coefficients annotated
############################################################
cor_pearson <- function(x,y,...){
  points(x,y)
  abline (h = 0, col = "grey")
  abline (v = 0, col = "grey")
  cor_res <- cor.test(x, y, method = 'pearson')
  r_value <- round(cor_res$estimate, 2)
  p_r <- formatC(cor_res$p.value, digits = 2, format = 'g')# Rounding the two-sided p values to two significant digits
  text(x = -5, y = 7, labels = paste(c('r=', r_value, '\n', 'p=', p_r), collapse= '') ,col='red')
}

############################################################
### Pairwise Pearson's correlation for six statins 
############################################################
drawing_Pearsons_pair <- function(in_z_file, in_cell) {
  ### Importing the dataframe
  in_z_df <- fread(in_z_file, data.table = FALSE)
  
  ### Removing pravastatin from analysis as it displayed much inconsistency with the other statins
  in_z_df <- in_z_df[ , c(str_subset(colnames(in_z_df), 'pravastatin', negate = TRUE))]
  
  ### Formatting the colnames
  for (curr_col in colnames(in_z_df)) {
    new_col <- str_split(curr_col, '_')[[1]][1]
    new_col <- str_to_title(new_col)
    colnames(in_z_df) <- str_replace(colnames(in_z_df), curr_col, new_col)
  } 
  
  ### Ordering the column alphabetically
  in_z_df <- in_z_df[, order(colnames(in_z_df))]
  
  ### Drawing the pairwise pearson's correlation plots and saving the plot
  out_pdf <- paste(c(Cor_plot_dir, 'Statins_', in_cell,'_all_genes_z_scores_pairwise_Pearsons.pdf'), collapse = '')
  pdf(out_pdf, width = 10, height = 10)
  pairs(in_z_df[, colnames(in_z_df) != 'V1'], lower.panel=NULL, panel = cor_pearson, xlim= c(-10,10) , ylim= c(-10, 10), main = in_cell, cex.axis=0.8)
  dev.off()

}


drawing_Pearsons_pair('./2_z_score_files/Statin_HA1E_10uM_24h_all_genes.txt', 'HA1E')
drawing_Pearsons_pair('./2_z_score_files/Statin_NPC_10uM_24h_all_genes.txt', 'NPC')
