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
### Annotating statin significant terms
#############################################################
### Each signficant BP term may trace back to more than one level-2 ancestor, and thus counted twice
statin_sig_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/5_gProfiler2_analysis/gProfiler2_results_gSCS_files_all_genes/Statins_HA1E_GO_BP_sig_in_atleast1_up_down.csv'
statin_sig_df <- fread(statin_sig_file)
statin_term_vector <- statin_sig_df$term_id
statin_sig_df$ancestor_level2 <- c('NA')
statin_sig_df$Term <- c('NA')

for (statin_term in statin_term_vector){
  curr_term_ancestor_vector <- as.list(GOBPANCESTOR[c(statin_term)])[[1]]
  annotated_ancestor <- c()
  for (curr_ancestor in c(curr_term_ancestor_vector, statin_term)){
    if (curr_ancestor %in% level_2_terms){
      annotated_ancestor <- c(annotated_ancestor, curr_ancestor)
    }
  }
  annotated_ancestor_string <- paste(annotated_ancestor, collapse = '|')
  statin_sig_df$ancestor_level2 <- ifelse(statin_sig_df$term_id == statin_term, annotated_ancestor_string, statin_sig_df$ancestor_level2)
  
  Curr_term <- Term(as.list(GOTERM[c(statin_term)])[[1]])
  statin_sig_df[which(statin_sig_df$term_id == statin_term), 'Term'] <- Curr_term
}

statin_sig_annotation_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_gProfiler2_HA1E_all_genes/Statins_HA1E_GO_BP_sig_in_atleast1_up_down_annotated_level2_ancestor.txt'
write.table(statin_sig_df, statin_sig_annotation_file, sep = '\t', quote = FALSE, row.names = FALSE)

#############################################################
### Counting the annotation frequency for each level 2 term
#############################################################
statin_sig_annotation_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_gProfiler2_HA1E_all_genes/Statins_HA1E_GO_BP_sig_in_atleast1_up_down_annotated_level2_ancestor.txt'
statin_sig_df <- fread(statin_sig_annotation_file)
total_GO_BP_count <- length(row.names(statin_sig_df))
statin_sig_annotation_vector_unsplit <- statin_sig_df$ancestor_level2
statin_sig_annotation_vector_split <- unlist(strsplit(statin_sig_annotation_vector_unsplit, '\\|'))
statin_sig_annotation_summary <- as.data.frame(table(statin_sig_annotation_vector_split))
colnames(statin_sig_annotation_summary) <- c('GO_BP_ID', 'Frequency')
statin_sig_annotation_summary$Term <- c('NA')

for (GO_BP in statin_sig_annotation_summary$GO_BP_ID){
  Curr_term <- Term(as.list(GOTERM[c(GO_BP)])[[1]])
  statin_sig_annotation_summary[which(statin_sig_annotation_summary$GO_BP_ID == GO_BP), 'Term'] <- Curr_term
}

statin_sig_annotation_summary_sort <- statin_sig_annotation_summary[order(statin_sig_annotation_summary$Frequency, decreasing = TRUE),]
statin_sig_annotation_summary_sort$Term <- ifelse(statin_sig_annotation_summary_sort$Term == 'localization', 'localisation', statin_sig_annotation_summary_sort$Term)

statin_sig_annotation_count_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_gProfiler2_HA1E_all_genes/Statins_HA1E_GO_BP_sig_in_atleast1_up_down_annotated_level2_ancestor_count.txt'
write.table(statin_sig_annotation_summary_sort, statin_sig_annotation_count_file, sep = '\t', quote = FALSE, row.names = FALSE)

#############################################################
### Drawing the frequency of each level 2 term
#############################################################
statin_sig_annotation_count_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_gProfiler2_HA1E_all_genes/Statins_HA1E_GO_BP_sig_in_atleast1_up_down_annotated_level2_ancestor_count.txt'
statin_sig_annotation_summary_sort <- fread(statin_sig_annotation_count_file)
statin_sig_annotation_summary_sort$Perc <- statin_sig_annotation_summary_sort$Frequency/total_GO_BP_count*100
statin_sig_annotation_summary_sort$label <- paste('(',statin_sig_annotation_summary_sort$GO_BP_ID, ')', '\n', statin_sig_annotation_summary_sort$Term, sep = '')

curr_fig <- ggplot(data=statin_sig_annotation_summary_sort, aes(x=factor(label, levels = statin_sig_annotation_summary_sort$label), y=Perc)) +
  geom_bar(stat="identity", fill='steelblue', alpha = 0.5)+
  ggtitle(label = 'GO Biological processes enriched in differentially expressed genes')+
  ylab('Percentage of significant biological pathway (%)')+
  geom_text(aes(label = Frequency), vjust = 0) +  
  ylim(0, 100)+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = 'black', fill=NA),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_blank())

curr_fig

out_fig_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_gProfiler2_HA1E_all_genes/Statins_HA1E_GO_BP_sig_in_atleast1_up_down_annotated_level2_ancestor_count.pdf'
ggsave(
  out_fig_file,
  plot = curr_fig,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 210,
  height = 100,
  units = c("mm"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL)


