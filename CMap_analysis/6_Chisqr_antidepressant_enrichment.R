### This script is used to perform chi-square test to test for the enrichment of antidepressants amongst compounds that showed high connectivity scores (average Tau > 90) to statins. The connectivity profiles were computed by submitting statin signatures (n most up-regulated and n down-regualted landmark genes). 

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
### Writing a function to perform chi square test to test for enrichment of antidepressants in compounds that show high Tau scores (Tau > 90).
##########################################################################
performing_chi <- function(in_tau_file, filename) { # the in_tau_file is the compound connectivity profiles computed by CLUE
  ### Importing the tau file, which contains the CP tau scores for the statins 
  in_tau_df <- fread(in_tau_file, data.table = FALSE)
  
  ### Finding antidepressants
  antidepressants_tau_df <- in_tau_df[in_tau_df$Compound_name %in% antidepressant_ATC_df$name, ]
  
  ### Performing the Chi-sqr test
  ### Counting the number of compounds with and without average Tau > 90
  cp_high_Tau <- in_tau_df[in_tau_df$Averaged_tau > 90, ]
  count_high_Tau_cp <- nrow(cp_high_Tau)
  cp_nonhigh_Tau <- in_tau_df[in_tau_df$Averaged_tau <= 90, ]
  count_nonhigh_Tau_cp <- nrow(cp_nonhigh_Tau)
  
  ### Counting the number of antidepressants with and without Tau > 90
  antidepressants_TS_high_Tau <- antidepressants_tau_df[antidepressants_tau_df$Averaged_tau > 90, ]
  count_high_Tau_antidep <- nrow(antidepressants_TS_high_Tau)
  antidepressants_TS_nonhigh_Tau <- antidepressants_tau_df[antidepressants_tau_df$Averaged_tau <= 90, ]
  count_nonhigh_Tau_antidep <- nrow(antidepressants_TS_nonhigh_Tau)
  
  ### Counting the total number of unique antidepressants
  uniq_count_high_Tau_antidep <- length(unique(antidepressants_TS_high_Tau$Compound_name))
  uniq_count_nonhigh_Tau_antidep <- length(unique(antidepressants_TS_nonhigh_Tau$Compound_name))
  
  ### constructing a contingency table
  contingency_table <- data.frame(antidepressant = c(count_high_Tau_antidep, count_nonhigh_Tau_antidep), non_antidepressant = c((count_high_Tau_cp-count_high_Tau_antidep), (count_nonhigh_Tau_cp-count_nonhigh_Tau_antidep)))
  row.names(contingency_table) <- c('high_tau', 'non_high_tau')
  
  test_result <- chisq.test(contingency_table) 
  
  test_result_file <- paste(c(enrichment_dir, filename, '_antidepressants_enrichment_ChiSqr.txt'), collapse = '')
  sink(test_result_file) # Start writing to an output file
  print(filename)
  print(contingency_table)
  print(test_result)
  print('Number of unique antidepressants with average Tau > 90')
  print(uniq_count_high_Tau_antidep)
  print('Number of unique antidepressants with average Tau <= 90')
  print(uniq_count_nonhigh_Tau_antidep)
  print('Warnings')
  print(warnings())
  sink() # Stop writing to the file	
}


performing_chi('./Statin_HA1E_10uM_24h_50_HA1E_cp_formatted.txt', 'Statin_HA1E_10uM_24h_50_HA1E')
performing_chi('./Statin_HA1E_10uM_24h_100_HA1E_cp_formatted.txt', 'Statin_HA1E_10uM_24h_100_HA1E')
performing_chi('./Statin_HA1E_10uM_24h_150_HA1E_cp_formatted.txt', 'Statin_HA1E_10uM_24h_150_HA1E')
performing_chi('./Statin_NPC_10uM_24h_50_summary_cp_formatted.txt', 'Statin_NPC_10uM_24h_50_summary')

