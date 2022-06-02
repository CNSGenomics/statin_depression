library(GO.db) # https://bioconductor.org/packages/release/data/annotation/manuals/GO.db/man/GO.db.pdf
library(data.table)
library(ggplot2)

#############################################################
### Getting the level 2 terms of biological process (GO:0008150)
#############################################################
BP_root_term <- c('GO:0008150')
level_2_terms_df <- as.data.frame(GOBPCHILDREN[c('GO:0008150')])
level_2_terms_df <- level_2_terms_df[which(level_2_terms_df$RelationshipType == 'isa'),]
level_2_terms <- level_2_terms_df$go_id

#############################################################
### Annotating GO terms for genes perturbed in the same direction and opposite directions
#############################################################
Ind_gProfiler2_dir <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/10_Common_DEGs_statin_antidepressant/Individual_statin_vs_antidepressant/Individual_statin_vs_antidepressant_gProfiler2/'

Creating_annotation_file <- function(in_pattern){
  GO_BP_ID_vector <- c()
  all_pattern_file_list <- list.files(Ind_gProfiler2_dir, pattern = in_pattern, full.names = TRUE)
  for (curr_file in all_pattern_file_list){
    if (!file.size(curr_file) < 10) {
      curr_df <- read.csv(curr_file, sep = '\t')
      curr_df_BP <- curr_df[curr_df$source == 'GO:BP',]
      curr_df_BP_terms <- curr_df_BP$term_id
      GO_BP_ID_vector <- c(GO_BP_ID_vector, curr_df_BP_terms)
    }
  } 
  GO_BP_ID_vector_unique <- unique(GO_BP_ID_vector)
  
  statin_antidep_BP_df <- as.data.frame(GO_BP_ID_vector_unique)
  colnames(statin_antidep_BP_df) <- c('GO_BP_ID')
  statin_antidep_BP_df$Term <- c('NA')
  statin_antidep_BP_df$ancestor_level2 <- c('NA')
  for (BP_ID in statin_antidep_BP_df$GO_BP_ID){
    curr_term_ancestor_vector <- as.list(GOBPANCESTOR[c(BP_ID)])[[1]]
    annotated_ancestor <- c()
    for (curr_ancestor in c(curr_term_ancestor_vector, BP_ID)){
      if (curr_ancestor %in% level_2_terms){
        annotated_ancestor <- c(annotated_ancestor, curr_ancestor)
      }
    }
    annotated_ancestor_string <- paste(annotated_ancestor, collapse = '|')
    statin_antidep_BP_df$ancestor_level2 <- ifelse(statin_antidep_BP_df$GO_BP_ID == BP_ID, annotated_ancestor_string, statin_antidep_BP_df$ancestor_level2)
    
    Curr_term <- Term(as.list(GOTERM[c(BP_ID)])[[1]])
    statin_antidep_BP_df[which(statin_antidep_BP_df$GO_BP_ID == BP_ID), 'Term'] <- Curr_term
  }
  
  BP_sig_annotation_file <- paste(c('/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_vs_antidepressant_gProfiler2_HA1E/Individual_statin_vs_antidepressant_gProfiler2_', in_pattern, '_annotated_level2_ancestor.txt'), collapse = '')
  write.table(statin_antidep_BP_df, BP_sig_annotation_file, sep = '\t', quote = FALSE, row.names = FALSE)
}

#############################################################
### Counting the annotation frequency for each level 2 term
#############################################################
counting_frequency <- function(in_pattern){
  BP_sig_annotation_file <- paste(c('/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_vs_antidepressant_gProfiler2_HA1E/Individual_statin_vs_antidepressant_gProfiler2_', in_pattern, '_annotated_level2_ancestor.txt'), collapse = '')
  
  statin_antidep_BP_df <- fread(BP_sig_annotation_file) 
  statin_antidep_annotation_vector_unsplit <- statin_antidep_BP_df$ancestor_level2
  statin_antidep_annotation_vector_split <- unlist(strsplit(statin_antidep_annotation_vector_unsplit, '\\|'))
  statin_antidep_annotation_summary <- as.data.frame(table(statin_antidep_annotation_vector_split))
  colnames(statin_antidep_annotation_summary) <- c('GO_BP_ID', 'Frequency')
  statin_antidep_annotation_summary$Term <- c('NA')
  
  for (GO_BP in statin_antidep_annotation_summary$GO_BP_ID){
    Curr_term <- Term(as.list(GOTERM[c(GO_BP)])[[1]])
    statin_antidep_annotation_summary[which(statin_antidep_annotation_summary$GO_BP_ID == GO_BP), 'Term'] <- Curr_term
  }
  
  statin_antidep_annotation_summary_sort <- statin_antidep_annotation_summary[order(statin_antidep_annotation_summary$Frequency, decreasing = TRUE),]
  statin_antidep_annotation_summary_sort$Term <- ifelse(statin_antidep_annotation_summary_sort$Term == 'localization', 'localisation', statin_antidep_annotation_summary_sort$Term)
  
  statin_antidep_annotation_count_file <- paste(c('/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_vs_antidepressant_gProfiler2_HA1E/Individual_statin_vs_antidepressant_gProfiler2_', in_pattern, '_annotated_level2_ancestor_count.txt'), collapse = '')
  write.table(statin_antidep_annotation_summary_sort, statin_antidep_annotation_count_file, sep = '\t', quote = FALSE, row.names = FALSE)
}



Creating_annotation_file('same_direction')
counting_frequency('same_direction')

Creating_annotation_file('opposite_direction')
counting_frequency('opposite_direction')
