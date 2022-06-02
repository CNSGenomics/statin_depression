library(GO.db) # https://bioconductor.org/packages/release/data/annotation/manuals/GO.db/man/GO.db.pdf
library(data.table)
library(ggplot2)

#############################################################
### Getting the level 2 terms of biological process (GO:0002376)
#############################################################
BP_root_term <- c('GO:0002376')
immune_process_terms_df <- as.data.frame(GOBPCHILDREN[c('GO:0002376')])
immune_process_terms_df <- immune_process_terms_df[which(immune_process_terms_df$RelationshipType == 'isa'),]
immune_process_terms <- immune_process_terms_df$go_id

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
  statin_antidep_BP_df$ancestor_immune_process <- c('NA')
  for (BP_ID in statin_antidep_BP_df$GO_BP_ID){
    curr_term_ancestor_vector <- as.list(GOBPANCESTOR[c(BP_ID)])[[1]]
    annotated_ancestor <- c()
    for (curr_ancestor in c(curr_term_ancestor_vector, BP_ID)){
      if (curr_ancestor %in% immune_process_terms){
        annotated_ancestor <- c(annotated_ancestor, curr_ancestor)
      }
    }
    annotated_ancestor_string <- paste(annotated_ancestor, collapse = '|')
    statin_antidep_BP_df$ancestor_immune_process <- ifelse(statin_antidep_BP_df$GO_BP_ID == BP_ID, annotated_ancestor_string, statin_antidep_BP_df$ancestor_immune_process)
    
    Curr_term <- Term(as.list(GOTERM[c(BP_ID)])[[1]])
    statin_antidep_BP_df[which(statin_antidep_BP_df$GO_BP_ID == BP_ID), 'Term'] <- Curr_term
  }
  
  BP_sig_annotation_file <- paste(c('/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_vs_antidepressant_gProfiler2_HA1E/Individual_statin_vs_antidepressant_gProfiler2_', in_pattern, '_annotated_immune_process_ancestor.txt'), collapse = '')
  write.table(statin_antidep_BP_df, BP_sig_annotation_file, sep = '\t', quote = FALSE, row.names = FALSE)
}

Creating_annotation_file('same_direction')
Creating_annotation_file('opposite_direction')
#############################################################
### Counting the annotation frequency for each ancestor term
#############################################################
counting_frequency <- function(in_pattern){
  BP_sig_annotation_file <- paste(c('/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_vs_antidepressant_gProfiler2_HA1E/Individual_statin_vs_antidepressant_gProfiler2_', in_pattern, '_annotated_immune_process_ancestor.txt'), collapse = '')
  
  statin_antidep_BP_df <- fread(BP_sig_annotation_file) 
  statin_antidep_annotation_vector_unsplit <- statin_antidep_BP_df$ancestor_immune_process
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
  
  statin_antidep_annotation_count_file <- paste(c('/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_vs_antidepressant_gProfiler2_HA1E/Individual_statin_vs_antidepressant_gProfiler2_', in_pattern, '_annotated_immune_process_ancestor_count.txt'), collapse = '')
  write.table(statin_antidep_annotation_summary_sort, statin_antidep_annotation_count_file, sep = '\t', quote = FALSE, row.names = FALSE)
}

counting_frequency('same_direction')
counting_frequency('opposite_direction')

#############################################################
### Drawing the plots side-by-side
#############################################################
getting_df <- function(in_pattern){
  BP_sig_annotation_file <- paste(c('/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_vs_antidepressant_gProfiler2_HA1E/Individual_statin_vs_antidepressant_gProfiler2_', in_pattern, '_annotated_immune_process_ancestor.txt'), collapse = '')
  statin_antidep_BP_df <- fread(BP_sig_annotation_file) 
  total_GO_BP_count <- length(row.names(statin_antidep_BP_df))
  print(total_GO_BP_count)
  
  statin_antidep_annotation_count_file <- paste(c('/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_vs_antidepressant_gProfiler2_HA1E/Individual_statin_vs_antidepressant_gProfiler2_', in_pattern, '_annotated_immune_process_ancestor_count.txt'), collapse = '')
  statin_antidep_annotation_summary_sort <- fread(statin_antidep_annotation_count_file)
  statin_antidep_annotation_summary_sort$Perc <- statin_antidep_annotation_summary_sort$Frequency/total_GO_BP_count*100
  statin_antidep_annotation_summary_sort$label <- paste('(',statin_antidep_annotation_summary_sort$GO_BP_ID, ')', '\n', statin_antidep_annotation_summary_sort$Term, sep = '')
  statin_antidep_annotation_summary_sort$group <- in_pattern
  statin_antidep_annotation_summary_sort_sep <- statin_antidep_annotation_summary_sort[,c('label', 'Perc')]
  
  return(statin_antidep_annotation_summary_sort_sep)
}

### Merging the dataframes
same_df <- getting_df('same_direction')
opposite_df <- getting_df('opposite_direction')

merged_BP <- merge(same_df, opposite_df, by = 'label', suffixes = c('_same', '_opposite'), all = TRUE)
merged_BP[is.na(merged_BP)] <- 0

### same direction 
same_df_filled <- merged_BP[, c('label', 'Perc_same')]
colnames(same_df_filled) <- c('label', 'Perc')
same_df_filled$group <- 'same direction'

### opposite direction
opposite_df_filled <- merged_BP[, c('label', 'Perc_opposite')]
colnames(opposite_df_filled) <- c('label', 'Perc')
opposite_df_filled$group <- 'opposite direction'

### ordering the BPs based on the frequency in the opposite group
opposite_df_filled_ordered <- opposite_df_filled[order(opposite_df_filled$Perc, decreasing = TRUE), ]
BP_ordered <- opposite_df_filled_ordered$label

### merging to the final dataframe for plotting
merged_df <- rbind(same_df_filled, opposite_df_filled)

curr_fig <- ggplot(data=merged_df, aes(x=factor(label, levels = BP_ordered), y=Perc, fill=group)) +
  geom_bar(stat="identity", position=position_dodge())+
  ylab('Percentage of significant biological pathway (%)')+
  ylim(0, 10)+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = 'black', fill=NA),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_blank())

out_fig_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_vs_antidepressant_gProfiler2_HA1E/Individual_statin_vs_antidepressant_gProfiler2_immune_process_ancestor_count.pdf'
ggsave(
  out_fig_file,
  plot = curr_fig,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 210,
  height = 150,
  units = c("mm"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL)
