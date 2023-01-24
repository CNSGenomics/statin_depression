### This script is used to compute the Pearsons' correlation between statin and antidepressant gene expression signatures

library(stringr)
library(data.table)
library(GSA)
library(gprofiler2)

Cor_dir <- './7_Statin_antidepressant_pairwise_z_Pearsons/'

############################################################
### Pairwise Pearson's correlation for six statins 
############################################################
drawing_Pearsons_pair <- function(in_statin_z_file, in_antidepressant_z_file, in_cell) {
  
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
  
  ### Ordering the rows based on gene IDs
  in_statin_z_df <- in_statin_z_df[order(in_statin_z_df$V1), ]
  in_antidepressant_z_df <- in_antidepressant_z_df[order(in_antidepressant_z_df$V1), ]

  ### Drawing the pairwise pearson's correlation plots and saving the plot
  out_pdf = paste(c(Cor_dir, 'Statins_antidepressants_', in_cell,'_all_genes_z_scores_pairwise_Pearsons.pdf'), collapse = '')
  pdf(out_pdf, width = 9.5, height = 10)
  par(mfrow=c((ncol(in_statin_z_df) - 1), (ncol(in_antidepressant_z_df) - 1)), mar=c(2.5,5,1.2,0.5),mgp=c(1.5,0.5,0))
  i <- 1
  while (i <= (ncol(in_statin_z_df) - 1)) {
    curr_statin <- colnames(in_statin_z_df)[i]
    j <- 1
    while (j <= (ncol(in_antidepressant_z_df) - 1)){
      curr_drug <- colnames(in_antidepressant_z_df)[j]

      print(paste(c(curr_statin, curr_drug), collapse = ' '))

      plot(in_statin_z_df[[curr_statin]], in_antidepressant_z_df[[curr_drug]], xlim = c(-10,10), ylim = c(-10,10), xlab=curr_statin, ylab=curr_drug)
      abline (h = 0, col = "grey")
      abline (v = 0, col = "grey")
      cor_res <- cor.test(in_statin_z_df[[curr_statin]], in_antidepressant_z_df[[curr_drug]], method = 'pearson')
      r_value <- round(cor_res$estimate, 2)
      p_r <- formatC(cor_res$p.value, 2, format = 'g')
      text(x = -5, y = 7, labels = paste(c('r=', r_value, '\n', 'p=', p_r), collapse= '') ,col='red')

      j <- j + 1
    }
    i <- i + 1
  }
  
  title(paste(c('Statin vs antidepressant ', in_cell, ' z-scores'), collapse = ''), line = -1, outer = TRUE)
  dev.off()
}


drawing_Pearsons_pair('./2_z_score_files/Statin_HA1E_10uM_24h_all_genes.txt', './2_z_score_files/Antidepresssants_HA1E_10uM_24h_all_genes.txt', 'HA1E')
drawing_Pearsons_pair('./2_z_score_files/Statin_NPC_10uM_24h_all_genes.txt', './2_z_score_files/Antidepressants_NPC_10uM_24h_all_genes.txt', 'NPC')
