### module load R/4.1.0+sf

library(stringr)
library(dplyr)
library(data.table)

in_file_list <-  c('/home/uqjjian3/Statin_SMR_checkpoints_10052022/HMGCR_rs12916_HMGCR_eQTLGEN_SMR_Result_collated.txt', '/home/uqjjian3/Statin_SMR_checkpoints_10052022/ITGAL_rs11574938_ITGAL_eQTLGEN_SMR_Result_collated.txt', '/home/uqjjian3/Statin_SMR_checkpoints_10052022/HDAC2_rs9481408_HDAC2_eQTLGEN_SMR_Result_collated.txt', '/home/uqjjian3/Statin_SMR_checkpoints_10052022/HMGCR_rs17671591_HMGCR_Brain_psychENCODE_SMR_Result_collated.txt', '/home/uqjjian3/Statin_SMR_checkpoints_10052022/PCSK9_rs12117661_PCSK9_Blood_GTEX_SMR_Result_collated.txt')

sig_thresh <- 0.05/(3*29) 
all_selected_traits <- c('baso_Vuckovic_2020_N', 'baso_p_Vuckovic_2020_N', 'eo_Vuckovic_2020_N', 'eo_p_Vuckovic_2020_N', 'hlr_Vuckovic_2020_N', 'hlr_p_Vuckovic_2020_N', 'irf_Vuckovic_2020_N', 'lymph_Vuckovic_2020_N', 'lymph_p_Vuckovic_2020_N', 'mpv_Vuckovic_2020_N', 'mrv_Vuckovic_2020_N', 'mono_Vuckovic_2020_N', 'mono_p_Vuckovic_2020_N', 'neut_Vuckovic_2020_N', 'neut_p_Vuckovic_2020_N', 'plt_Vuckovic_2020_N', 'pct_Vuckovic_2020_N', 'pdw_Vuckovic_2020_N', 'rbc_Vuckovic_2020_N', 'rdw_cv_Vuckovic_2020_N', 'ret_Vuckovic_2020_N', 'ret_p_Vuckovic_2020_N', 'wbc_Vuckovic_2020_N', 'IL6_ahola_olli', 'CRP_Han_GCST009777', 'jointGwasMc_HDL', 'jointGwasMc_LDL', 'jointGwasMc_TG', 'CAD_VDH_UKB_ebiaGCST005194', 'T2D_Xue_ebiaGCST006867', 'BMI_UKB_ukbb19953', 'PGC_UKB_depression_no23andme', 'neuroticism_Nagel_2018', 'depressed_affect_Nagel_2018_NatGenet', 'worry_Nagel_2018_NatGenet')

all_selected_traits_reform <- c('Basophil count', 'Basophil percentage of white cells', 'Eosinophil count', 'Eosinophil percentage of white cells', 'High light scatter reticulocyte count', 'High light scatter reticulocyte percentage of red cells', 'Immature fraction of reticulocytes', 'Lymphocyte count', 'Lymphocyte percentage of white cells', 'Mean platelet volume', 'Mean reticulocyte volume', 'Monocyte count', 'Monocyte percentage of white cells', 'Neutrophil count', 'Neutrophil percentage of white cells', 'Platelet count', 'Platelet crit', 'Platelet distribution width', 'Red blood cell count', 'Red cell distribution width', 'Reticulocyte count', 'Reticulocyte fraction of red cells', 'White blood cell count', 'Blood IL6 levels', 'Serum CRP levels', 'HDL-C', 'LDL-C', 'TG', 'CAD', 'T2D', 'BMI', 'Depression',  'Neuroticism', 'Depressed affect', 'Worry')

all_control_selected_traits_reform <- c('HDL-C', 'LDL-C', 'TG', 'CAD', 'T2D', 'BMI')

#############################################
### for the SMR data
#############################################
for (in_file in in_file_list){
  print(in_file)
  
  ### Importing the df
  smr_df_original <- fread(in_file)
  smr_df <- smr_df_original[, c('trait', 'Gene', 'targetSNP', 'A1', 'A2', 'Freq', 'b_GWAS', 'se_GWAS', 'p_GWAS', 'b_eQTL', 'se_eQTL', 'p_eQTL', 'b_SMR', 'se_SMR', 'p_SMR', 'p_HEIDI', 'nsnp_HEIDI', 'b_SMR_perSDdecrease')]
  colnames(smr_df) <- c('Trait', 'Gene', 'Genetic_instrument', 'Effect_allele', 'Other_allele', 'Effect_allele_freq', 'b_GWAS', 'se_GWAS', 'p_GWAS', 'b_eQTL', 'se_eQTL', 'p_eQTL', 'b_SMR', 'se_SMR', 'p_SMR', 'p_HEIDI', 'nsnp_HEIDI', 'b_SMR_perSDdecrease')
  smr_df$F_stat <- (as.numeric(smr_df$b_eQTL)^2)/(as.numeric(smr_df$se_eQTL)^2)
  smr_df <- smr_df[, c('Trait', 'Gene', 'Genetic_instrument', 'Effect_allele', 'Other_allele', 'Effect_allele_freq', 'b_GWAS', 'se_GWAS', 'p_GWAS', 'b_eQTL', 'se_eQTL', 'p_eQTL', 'F_stat', 'b_SMR', 'se_SMR', 'p_SMR', 'p_HEIDI', 'nsnp_HEIDI', 'b_SMR_perSDdecrease')]
  
  ### Filtering for selected traits
  matching_index <- c()
  i <- 1
  while (i <= length(all_selected_traits)){
    if (all_selected_traits[i] %in% smr_df$Trait) {
      matching_index <- c(matching_index, i)
    }
    i <- i + 1
  }
  
  selected_traits <- all_selected_traits[matching_index]
  selected_traits_reform <- all_selected_traits_reform[matching_index]
  
  print(length(selected_traits))
  print(length(selected_traits_reform))
  
  ### renaming traits
  smr_df_sel <- smr_df[smr_df$Trait %in% selected_traits, ]
  smr_df_sel <- smr_df_sel[order(factor(smr_df_sel$Trait, levels = selected_traits)),]
  smr_df_sel$Trait <- selected_traits_reform
  
  ### Calculating error bars
  smr_df_sel$b_SMR_perSDdecrease <- as.numeric(smr_df_sel$b_SMR_perSDdecrease)
  smr_df_sel$se_SMR <- as.numeric(smr_df_sel$se_SMR)
  smr_df_sel$Upper_CI_beta <- smr_df_sel$b_SMR_perSDdecrease + smr_df_sel$se_SMR * 1.96
  smr_df_sel$Lower_CI_beta <- smr_df_sel$b_SMR_perSDdecrease - smr_df_sel$se_SMR * 1.96
  
  smr_df_sel_outfile <- str_replace(in_file, '.txt', '_organised.txt')
  write.table(smr_df_sel, smr_df_sel_outfile, sep = '\t', row.names = FALSE)

}

