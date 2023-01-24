### This script is used to analyse the pre-compiled connectivity profiles of "HMGCR inhibitor" downloaded from the CLUE platform. The distributions of Tau scores amongst statins, and between statin and non-statin compounds are plotted separately for each cell line, in order to identify a cell line where statins induced the most specific gene expression signatures

library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)

out_fig_dir <- './1_CLUE_conn_HMGCR_inhibitor/'
Touchstone_gct <- './1_CLUE_conn_HMGCR_inhibitor/conn_download/conn_HMGCR_inhibitor.gct'

### importing the dataframe and removing unnecessary columns (e.g. description)
Touchstone_gct_df <- fread(Touchstone_gct)
Touchstone_gct_df <- Touchstone_gct_df[-c(1,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), -c(2,4,5,6,7,8,9)]

### getting a list of statinS
statin_vector <- unique(as.character(Touchstone_gct_df[Touchstone_gct_df$`2429` == 'name',]))
statin_vector <- statin_vector[!statin_vector %in% c('na', 'name')]

### getting a list of cell type
cell_vector <- unique(as.character(Touchstone_gct_df[Touchstone_gct_df$`2429` == 'cell_id',]))
cell_vector <- cell_vector[!cell_vector %in% c('cell_id', 'na', 'summary')]

Touchstone_gct_df <- as.data.frame(Touchstone_gct_df)

drawing_density <- function(in_cell_id, in_statin_vector){
  ### selecting the right statins
  sel_Touchstone_gct_df <- Touchstone_gct_df[, which(Touchstone_gct_df[2, ] %in% c(in_statin_vector, 'name', 'na'))]
  
  # setting the colnames
  colnames(sel_Touchstone_gct_df) <- sel_Touchstone_gct_df[1,]
  sel_Touchstone_gct_df <- sel_Touchstone_gct_df[-1, ]
  
  ### selecting the right cell line
  Touchstone_gct_df_cell <- sel_Touchstone_gct_df[, (colnames(sel_Touchstone_gct_df) %in% c(in_cell_id, 'na'))]
  header_vector <- paste(as.character(Touchstone_gct_df_cell[1,]), as.character(colnames(Touchstone_gct_df_cell)), sep = '_')
  header_vector[1] <- 'name'
  colnames(Touchstone_gct_df_cell) <- header_vector
  Touchstone_gct_df_cell <- Touchstone_gct_df_cell[-1,]
  
  ### Separating the statins and non-statin compounds
  Touchstone_gct_df_cp_statin <- Touchstone_gct_df_cell[which(Touchstone_gct_df_cell$name %in% in_statin_vector),]
  Touchstone_gct_df_cp_nonstatin <- Touchstone_gct_df_cell[-which(Touchstone_gct_df_cell$name %in% in_statin_vector),]
  
  ### getting the tau scores for statin compounds
  statin_Tau_vector <- as.numeric(unlist(Touchstone_gct_df_cp_statin[2:length(colnames(Touchstone_gct_df_cp_statin))]))
  statin_Tau_vector <- statin_Tau_vector[!is.na(statin_Tau_vector)]
  statin_Tau_df <- data.frame(Tau_score = statin_Tau_vector, Type = 'statin')
  
  ### getting the tau scores for non-statin compounds
  nonstatin_Tau_vector <- as.numeric(unlist(Touchstone_gct_df_cp_nonstatin[2:length(colnames(Touchstone_gct_df_cp_nonstatin))]))
  nonstatin_Tau_df <- data.frame(Tau_score = nonstatin_Tau_vector, Type = 'non-statin')
  
  ### stacking the two dfs 
  Tau_df_cp_ann <- rbind(statin_Tau_df, nonstatin_Tau_df)
  
  ### plotting the density plot
  curr_plot <- ggplot(Tau_df_cp_ann, aes(x = Tau_score, fill = Type)) + 
    geom_density(alpha = 0.5)+
    theme_bw()+
    xlim(-100, 100) +
    ggtitle(in_cell_id)+
    theme(plot.title = element_text(size = 10),
          legend.key.size = unit(0.2, 'cm'), #change legend key size
          legend.title = element_blank(), #change legend title font size
          legend.text = element_text(size=10), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  return(curr_plot)
}

### plotting for all 9 statins
plot1 <- drawing_density('A375', statin_vector)
plot2 <- drawing_density('A549', statin_vector)
plot3 <- drawing_density('HA1E', statin_vector)
plot4 <- drawing_density('HCC515', statin_vector)
plot5 <- drawing_density('HEPG2', statin_vector)
plot6 <- drawing_density('HT29', statin_vector)
plot7 <- drawing_density('MCF7', statin_vector)
plot8 <- drawing_density('PC3', statin_vector)
plot9 <- drawing_density('VCAP', statin_vector)

total_plot <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, labels = "AUTO", nrow = 3, ncol = 3)
total_plot
out_fig_file <- paste(c(out_fig_dir, 'Touchstone_statin_vs_non_statin_ggplot_density.pdf'), collapse = '')
ggsave(
  out_fig_file,
  plot = total_plot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 297,
  height = 210,
  units = c("mm"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL)


### plotting for the six statins studied in this study
selected_statins <- c("atorvastatin", "lovastatin", "mevastatin", "fluvastatin", "simvastatin", "rosuvastatin")
plot1 <- drawing_density('A375', selected_statins)
plot2 <- drawing_density('A549', selected_statins)
plot3 <- drawing_density('HA1E', selected_statins)
plot4 <- drawing_density('HCC515', selected_statins)
plot5 <- drawing_density('HEPG2', selected_statins)
plot6 <- drawing_density('HT29', selected_statins)
plot7 <- drawing_density('MCF7', selected_statins)
plot8 <- drawing_density('PC3', selected_statins)
plot9 <- drawing_density('VCAP', selected_statins)

selected_plot <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, labels = "AUTO", nrow = 3, ncol = 3)
total_plot
out_fig_file <- paste(c(out_fig_dir, 'Touchstone_statin_vs_non_statin_selected_ggplot_density.pdf'), collapse = '')
ggsave(
  out_fig_file,
  plot = selected_plot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 297,
  height = 210,
  units = c("mm"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL)
