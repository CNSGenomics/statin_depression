### This script is used to perform chi-square test to test for the enrichment of antidepressants amongst compounds that showed high connectivity scores (average Tau > 90) to statins. The connectivity profiles were computed by submitting statin signatures (n most up-regulated and n down-regualted landmark genes). 
### If a drug has been profiled more than once by CMap (with more than one compound ID), we only retain the ID profiled in the TouchStone collection (See: https://clue.io/data/TS).

library(data.table)
library(stringr)

##########################################################################
### This script is written to perform chi square test to test for enrichment of antidepressants in compounds that show high Tau scores (Tau > 90).
##########################################################################
enrichment_dir <- './6_Chisqr_antidepressant_enrichment/'

##########################################################################
### Getting a list of antidepressant names
##########################################################################
psych_Drug_ATC_file <- './Data/ATC_drugs/ATC_classification_psyc_drugs.txt'
psych_Drug_ATC_df <- fread(psych_Drug_ATC_file, header = FALSE, col.names = c('ATC_code', 'name'))
antidepressant_ATC_df <- psych_Drug_ATC_df[psych_Drug_ATC_df$ATC_code %in% str_subset(psych_Drug_ATC_df$ATC_code, 'N06A'), ]

##########################################################################
### Getting a list of compound IDs in the TouchStone collection
##########################################################################
TS_pert_info_file <- './Data/LINCS/GSE92742_Broad_LINCS_pert_info.txt'
TS_pert_info_df <- fread(TS_pert_info_file)

### filtering for compounds
cp_pert_info_df <- TS_pert_info_df[TS_pert_info_df$pert_type == 'trt_cp', ]

### filtering for compound IDs that are in the TouchStone collection
cp_TS_df <- cp_pert_info_df[cp_pert_info_df$is_touchstone == '1', ]
cp_TS_ID <- cp_TS_df$pert_id

##########################################################################
### Writing a function to perform chi square test to test for enrichment of antidepressants in compounds that show high Tau scores (Tau > 90).
##########################################################################
performing_chi <- function(in_tau_file, filename) { # the in_tau_file is the compound connectivity profiles computed by CLUE
  ### Importing the tau file, which contains the CP tau scores for the statins 
  in_tau_df <- fread(in_tau_file, data.table = FALSE)

  ### If a drug (e.g. "maprotline") has more than one compound ID in the CMap database (e.g. different brands), we retain only the ID in the TouchStone collected of well profiled compounds to avoid bias in results
  ### finding a list of compounds that are duplicates
  duplicate_cp <- unique(in_tau_df$Compound_name[duplicated(in_tau_df$Compound_name)])
  
  ### marking their duplicate status in the dataframe
  in_tau_df$duplicate_status <- ifelse(in_tau_df$Compound_name %in% duplicate_cp, '1', '0')
  
  ### removing compounds that are duplicates AND are not in the TouchStone collection
  in_tau_df$cp_removed <- ifelse(((in_tau_df$duplicate_status == '1') & (!in_tau_df$Compound_ID %in% cp_TS_ID)), 'removed', 'retained') ### removing duplicates that are not in the TouchStone
  filtered_tau_df <- in_tau_df[in_tau_df$cp_removed == 'retained', ]
  
  ### saving the TS tau profile
  filtered_tau_df$duplicate_status <- NULL
  filtered_tau_df$cp_removed <- NULL
  fil_tau_file <- str_replace(in_tau_file, '.txt', '_TS_filtered.txt')
  print(fil_tau_file)
  write.table(filtered_tau_df, fil_tau_file, sep = '\t', quote = FALSE, row.names = FALSE)
  
  ### Finding antidepressants
  antidepressants_tau_df <- filtered_tau_df[filtered_tau_df$Compound_name %in% antidepressant_ATC_df$name, ]
  
  ### Performing the Chi-sqr test
  ### Counting the number of compounds with and without average Tau > 90
  cp_high_Tau <- filtered_tau_df[filtered_tau_df$Averaged_tau > 90, ]
  count_high_Tau_cp <- nrow(cp_high_Tau)
  cp_nonhigh_Tau <- filtered_tau_df[filtered_tau_df$Averaged_tau <= 90, ]
  count_nonhigh_Tau_cp <- nrow(cp_nonhigh_Tau)
  
  ### Counting the number of antidepressants with and without Tau > 90
  antidepressants_TS_high_Tau <- antidepressants_tau_df[antidepressants_tau_df$Averaged_tau > 90, ]
  count_high_Tau_antidep <- nrow(antidepressants_TS_high_Tau)
  antidepressants_TS_nonhigh_Tau <- antidepressants_tau_df[antidepressants_tau_df$Averaged_tau <= 90, ]
  count_nonhigh_Tau_antidep <- nrow(antidepressants_TS_nonhigh_Tau)

  ### constructing a contingency table
  contingency_table <- data.frame(antidepressant = c(count_high_Tau_antidep, count_nonhigh_Tau_antidep), non_antidepressant = c((count_high_Tau_cp-count_high_Tau_antidep), (count_nonhigh_Tau_cp-count_nonhigh_Tau_antidep)))
  row.names(contingency_table) <- c('high_tau', 'non_high_tau')
  
  ### running the chi square test
  test_result <- chisq.test(contingency_table) 
  print(test_result)
  
  ### saving the results
  test_result_file <- paste(c(enrichment_dir, filename, '_antidepressants_enrichment_ChiSqr_TouchStone.txt'), collapse = '')
  sink(test_result_file) # Start writing to an output file
  print(filename)
  print(contingency_table)
  print(test_result)
  print('Warnings')
  print(warnings())
  sink() # Stop writing to the file	
}


performing_chi('./5_CLUE_query/CLUE_results_formatted/Statin_HA1E_10uM_24h_50_19Jan2023_HA1E_cp_formatted.txt', 'Statin_HA1E_10uM_24h_50_HA1E')
performing_chi('./5_CLUE_query/CLUE_results_formatted/Statin_HA1E_10uM_24h_100_19Jan2023_HA1E_cp_formatted.txt', 'Statin_HA1E_10uM_24h_100_HA1E')
performing_chi('./5_CLUE_query/CLUE_results_formatted/Statin_HA1E_10uM_24h_150_19Jan2023_HA1E_cp_formatted.txt', 'Statin_HA1E_10uM_24h_150_HA1E')
performing_chi('./5_CLUE_query/CLUE_results_formatted/Statin_NPC_10uM_24h_50_19Jan2023_summary_cp_formatted.txt', 'Statin_NPC_10uM_24h_50_summary')

