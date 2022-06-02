library(GO.db) # https://bioconductor.org/packages/release/data/annotation/manuals/GO.db/man/GO.db.pdf
library(data.table)
library(ggplot2)
library(GSA)

#############################################################
### Getting the child terms of primary metabolic process (GO:0044238)
#############################################################
primary_meta_process_term <- c('GO:0044238')
primary_meta_process_sub_terms_df <- as.data.frame(GOBPCHILDREN[primary_meta_process_term])
colnames(primary_meta_process_sub_terms_df) <- c('child', 'parent', 'RelationshipType')
primary_meta_process_sub_terms_df <- primary_meta_process_sub_terms_df[which(primary_meta_process_sub_terms_df$RelationshipType == 'isa'),]
primary_meta_process_sub_terms <- primary_meta_process_sub_terms_df$child
primary_meta_process_sub_terms_df$desc <- 'NA'

########################################################
### 2. Getting the GO BP gene entrez IDs
########################################################
gProfiler_BP_gmt <- '/Users/uqjjian3/Data/gProfiler/gprofiler_hsapiens.ENSG/hsapiens.GO_BP.ENSG.gmt'
gProfiler_list <- GSA.read.gmt(gProfiler_BP_gmt)
gProfiler_genesets<- gProfiler_list$genesets
gProfiler_geneset_names <- gProfiler_list$geneset.names
gProfiler_geneset_desc <- gProfiler_list$geneset.descriptions

for (in_GO_BP_ID in primary_meta_process_sub_terms) {
  in_GO_BP_index <- match(in_GO_BP_ID, gProfiler_geneset_names)
  in_GO_BP_desc <- gProfiler_geneset_desc[[in_GO_BP_index]]
  primary_meta_process_sub_terms_df[which(primary_meta_process_sub_terms_df$child == in_GO_BP_ID), 'desc'] <- in_GO_BP_desc
}

#############################################################
### Annotating statin significant terms
#############################################################
### Each signficant BP term may trace back to more than one ancestor ancestor, and thus counted twice
statin_sig_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/5_gProfiler2_analysis/gProfiler2_results_gSCS_files_all_genes/Statins_HA1E_GO_BP_sig_in_atleast1_up_down.csv'
statin_sig_df <- fread(statin_sig_file)
statin_term_vector <- statin_sig_df$term_id
statin_sig_df$ancestor_primary_meta <- c('NA')
statin_sig_df$child_Term <- c('NA')

for (statin_term in statin_term_vector){
  curr_term_ancestor_vector <- as.list(GOBPANCESTOR[c(statin_term)])[[1]]
  annotated_ancestor <- c()
  for (curr_ancestor in c(curr_term_ancestor_vector, statin_term)){
    if (curr_ancestor %in% primary_meta_process_sub_terms){
      annotated_ancestor <- c(annotated_ancestor, curr_ancestor)
    }
  }
  annotated_ancestor_string <- paste(annotated_ancestor, collapse = '|')
  statin_sig_df[which(statin_sig_df$term_id == statin_term), 'ancestor_primary_meta'] <- annotated_ancestor_string
  
  Curr_term <- Term(as.list(GOTERM[c(statin_term)])[[1]])
  statin_sig_df[which(statin_sig_df$term_id == statin_term), 'child_Term'] <- Curr_term
}

statin_sig_annotation_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_gProfiler2_HA1E_all_genes/Statins_HA1E_GO_BP_sig_in_atleast1_up_down_annotated_primary_metabolic_ancestor.txt'
write.table(statin_sig_df, statin_sig_annotation_file, sep = '\t', quote = FALSE, row.names = FALSE)

#############################################################
### Counting the annotation frequency for each ancestor term
#############################################################
statin_sig_annotation_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_gProfiler2_HA1E_all_genes/Statins_HA1E_GO_BP_sig_in_atleast1_up_down_annotated_primary_metabolic_ancestor.txt'
statin_sig_df <- fread(statin_sig_annotation_file)
total_GO_BP_count <- length(row.names(statin_sig_df))
statin_sig_annotation_vector_unsplit <- statin_sig_df$ancestor_primary_meta
statin_sig_annotation_vector_split <- unlist(strsplit(statin_sig_annotation_vector_unsplit, '\\|'))
statin_sig_annotation_summary <- as.data.frame(table(statin_sig_annotation_vector_split))
colnames(statin_sig_annotation_summary) <- c('primary_meta_GO_BP_ID', 'Frequency')
statin_sig_annotation_summary$Term <- c('NA')

for (GO_BP in statin_sig_annotation_summary$primary_meta_GO_BP_ID){
  Curr_term <- Term(as.list(GOTERM[c(GO_BP)])[[1]])
  statin_sig_annotation_summary[which(statin_sig_annotation_summary$primary_meta_GO_BP_ID == GO_BP), 'Term'] <- Curr_term
}

statin_sig_annotation_count_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/25_Annotating_GO_BP_levels/statin_gProfiler2_HA1E_all_genes/Statins_HA1E_GO_BP_sig_in_atleast1_up_down_annotated_primary_metabolic_ancestor_count.txt'
write.table(statin_sig_annotation_summary, statin_sig_annotation_count_file, sep = '\t', quote = FALSE, row.names = FALSE)
