### module load R/4.1.0+sf

library(ggplot2)
library(stringr)
library(dplyr)
# library(egg)
library(ggpubr)
library(ggforce)

in_file_list <-  c('/home/uqjjian3/Statin_SMR_checkpoints_10052022/HMGCR_rs12916_HMGCR_eQTLGEN_SMR_Result_collated.txt', '/home/uqjjian3/Statin_SMR_checkpoints_10052022/ITGAL_rs11574938_ITGAL_eQTLGEN_SMR_Result_collated.txt', '/home/uqjjian3/Statin_SMR_checkpoints_10052022/HDAC2_rs9481408_HDAC2_eQTLGEN_SMR_Result_collated.txt', '/home/uqjjian3/Statin_SMR_checkpoints_10052022/HMGCR_rs17671591_HMGCR_Brain_psychENCODE_SMR_Result_collated.txt', '/home/uqjjian3/Statin_SMR_checkpoints_10052022/PCSK9_rs12117661_PCSK9_Blood_GTEX_SMR_Result_collated.txt')

sig_thresh <- 0.05/(3*29) 
all_selected_traits <- c('baso_Vuckovic_2020_N', 'baso_p_Vuckovic_2020_N', 'eo_Vuckovic_2020_N', 'eo_p_Vuckovic_2020_N', 'hlr_Vuckovic_2020_N', 'hlr_p_Vuckovic_2020_N', 'irf_Vuckovic_2020_N', 'lymph_Vuckovic_2020_N', 'lymph_p_Vuckovic_2020_N', 'mpv_Vuckovic_2020_N', 'mrv_Vuckovic_2020_N', 'mono_Vuckovic_2020_N', 'mono_p_Vuckovic_2020_N', 'neut_Vuckovic_2020_N', 'neut_p_Vuckovic_2020_N', 'plt_Vuckovic_2020_N', 'pct_Vuckovic_2020_N', 'pdw_Vuckovic_2020_N', 'rbc_Vuckovic_2020_N', 'rdw_cv_Vuckovic_2020_N', 'ret_Vuckovic_2020_N', 'ret_p_Vuckovic_2020_N', 'wbc_Vuckovic_2020_N', 'IL6_ahola_olli', 'CRP_Han_GCST009777', 'jointGwasMc_HDL', 'jointGwasMc_LDL', 'jointGwasMc_TG', 'CAD_VDH_UKB_ebiaGCST005194', 'T2D_Xue_ebiaGCST006867', 'BMI_UKB_ukbb19953', 'PGC_UKB_depression_no23andme', 'neuroticism_Nagel_2018', 'depressed_affect_Nagel_2018_NatGenet', 'worry_Nagel_2018_NatGenet')

all_selected_traits_reform <- c('Basophil count', 'Basophil percentage of white cells', 'Eosinophil count', 'Eosinophil percentage of white cells', 'High light scatter reticulocyte count', 'High light scatter reticulocyte percentage of red cells', 'Immature fraction of reticulocytes', 'Lymphocyte count', 'Lymphocyte percentage of white cells', 'Mean platelet volume', 'Mean reticulocyte volume', 'Monocyte count', 'Monocyte percentage of white cells', 'Neutrophil count', 'Neutrophil percentage of white cells', 'Platelet count', 'Platelet crit', 'Platelet distribution width', 'Red blood cell count', 'Red cell distribution width', 'Reticulocyte count', 'Reticulocyte fraction of red cells', 'White blood cell count', 'Blood IL6 levels', 'Serum CRP levels', 'HDL-C', 'LDL-C', 'TG', 'CAD', 'T2D', 'BMI', 'Depression',  'Neuroticism', 'Depressed affect', 'Worry')

all_control_selected_traits_reform <- c('HDL-C', 'LDL-C', 'TG', 'CAD', 'T2D', 'BMI')

#############################################
### for the eQTLGEN data
#############################################
for (in_file in in_file_list){
  print(in_file)
  
  ### Importing the df
  smr_df <- read.table(in_file, header = TRUE, sep = '\t')
  
  ### Filtering for existing traits
  matching_index <- c()
  i <- 1
  while (i <= length(all_selected_traits)){
    if (all_selected_traits[i] %in% smr_df$trait) {
      matching_index <- c(matching_index, i)
    }
    i <- i + 1
  }
  
  selected_traits <- all_selected_traits[matching_index]
  selected_traits_reform <- all_selected_traits_reform[matching_index]
  
  control_selected_traits_reform <- c()
  for (trait in all_control_selected_traits_reform){
    if (trait %in% selected_traits_reform) {
      control_selected_traits_reform <- c(control_selected_traits_reform, trait)
    }
  }
  
  print(length(selected_traits))
  print(length(selected_traits_reform))

  ### renaming traits
  smr_df_sel <- smr_df[smr_df$trait %in% selected_traits, ]
  smr_df_sel <- smr_df_sel[order(factor(smr_df_sel$trait, levels = selected_traits)),]
  smr_df_sel$trait <- selected_traits_reform

  ### Calculating error bars
  smr_df_sel$b_SMR_perSDdecrease <- as.numeric(smr_df_sel$b_SMR_perSDdecrease)
  smr_df_sel$se_SMR <- as.numeric(smr_df_sel$se_SMR)
  smr_df_sel$Upper_CI_beta <- smr_df_sel$b_SMR_perSDdecrease + smr_df_sel$se_SMR * 1.96
  smr_df_sel$Lower_CI_beta <- smr_df_sel$b_SMR_perSDdecrease - smr_df_sel$se_SMR * 1.96

  ### Separating control traits
  smr_df_sel$trait_type <- ifelse(smr_df_sel$trait %in% control_selected_traits_reform, 'Proof-of-principle traits', 'Outcome traits of interest')
  smr_df_sel$trait_type <- factor(smr_df_sel$trait_type, levels = c('Proof-of-principle traits', 'Outcome traits of interest'))

  ### Annotating the significance for control
  smr_df_sel_mul_sig <- subset(smr_df_sel, smr_df_sel$p_SMR < sig_thresh)
  smr_df_sel_nom_sig <- subset(smr_df_sel, smr_df_sel$p_SMR < 0.05)
  smr_df_sel_sig_HEIDI <- subset(smr_df_sel, smr_df_sel$p_HEIDI < 0.01)


  ### Drawing figures
  curr_fig <- ggplot(smr_df_sel, aes(x = b_SMR_perSDdecrease, y = factor(trait, levels=rev(selected_traits_reform))), alpha = 0.5) +
    geom_vline(aes(xintercept = 0), size = 0.4, linetype = 'dashed') +
    geom_errorbarh(aes(xmax = Upper_CI_beta, xmin = Lower_CI_beta), size = .5, height =
                     .2, color = 'gray50') +
    geom_point(size = 1.5, color = 'black') +
    geom_point(data=smr_df_sel_nom_sig, colour="gold1", size = 1.5) +
    geom_point(data=smr_df_sel_mul_sig, colour="red", size = 1.5) +
    geom_point(data=smr_df_sel_sig_HEIDI, shape = 4, size = 4 ) +
    theme_bw()+
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    ylab('') +
    xlab('Beta (unit change)') +
    ggforce::facet_col(.~trait_type, scales = 'free_y', space = 'free')

  out_fig_file <- str_replace(in_file, '.txt', '_beta_multiple_testing.pdf')
  ggsave(
    out_fig_file,
    plot = curr_fig,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 150,
    height = 150,
    units = c('mm'),
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
  )

}

