### This script is used to annotate significantly enriched biological processes with their ancestors. As the results of gProfiler2 often contain redundant and overlapping GO:BP terms, we have grouped these GO:BP terms based on their ancestor terms.

library(GO.db) # https://bioconductor.org/packages/release/data/annotation/manuals/GO.db/man/GO.db.pdf
library(data.table)
library(ggplot2)
library(GSA)

#############################################################
### Getting list of pathways associated with each ancestor pathway
#############################################################
getting_child_terms <- function(in_ancestor_term) {
  child_term_df  <- as.data.frame(GOBPCHILDREN[c(in_ancestor_term)])
  child_term_df <- child_term_df[child_term_df$RelationshipType == 'isa', ]
  child_term <- child_term_df$go_id
  return(child_term)
}

root_child_terms <- getting_child_terms('GO:0008150')
immune_child_terms <- getting_child_terms('GO:0002376')
metabolic_child_terms <- getting_child_terms('GO:0044238')

#############################################################
### A function to check and annotate the ancestors
#############################################################
annotating_ancestor <- function(in_BP, in_ancestor_vector) {
  ### Getting the ancestors of the current term
  curr_BP_ancestor_vector <- as.list(GOBPANCESTOR[c(in_BP)])[[1]]
  
  ### Adding the current term to the vector in case that it is one of the ancestor terms
  curr_BP_ancestor_vector <- c(curr_BP_ancestor_vector, in_BP)
  
  ### Now looking for its ancestor terms 
  annotated_ancestor_vector <- intersect(in_ancestor_vector, curr_BP_ancestor_vector)
  annotated_ancestor_str <- paste(annotated_ancestor_vector, collapse = '|')
  return(annotated_ancestor_str)
}

#############################################################
### Annotating the pathways for their ancestors
#############################################################
Annotating_df <- function(in_file) {
  ### Importing df 
  in_df <- fread(in_file)
  if (length(colnames(in_df)) > 1){ # check if file is not empty
    
    ### Filtering for GO:BP terms
    BP_df <- in_df[in_df$source == 'GO:BP', ]
    
    ### Annotating the statin type, antidepressant type, cell type and the direction of DEGs
    fname <- str_split(in_file, '/')[[1]][length(str_split(in_file, '/')[[1]])]
    fname <- str_split(fname, '\\.')[[1]][1]
    curr_statin <- str_split(fname, '_')[[1]][1]
    curr_antidepressant <- str_split(fname, '_')[[1]][2]
    curr_cell <- str_split(fname, '_')[[1]][3]
    curr_direction <- str_split(fname, '_gProfiler2_')[[1]][length(str_split(fname, '_gProfiler2_')[[1]])]
    BP_df$Statin <- curr_statin
    BP_df$Antidepressant <- curr_antidepressant
    BP_df$Cell <- curr_cell
    BP_df$Direction_DEGs <- curr_direction
    
    ### Annotating the ancestor terms
    BP_df$BP_ancestor <- lapply(BP_df$term_id, annotating_ancestor, in_ancestor_vector=root_child_terms)
    BP_df$Immune_ancestor <- lapply(BP_df$term_id, annotating_ancestor, in_ancestor_vector=immune_child_terms)
    BP_df$Primary_metabolic_ancestor <- lapply(BP_df$term_id, annotating_ancestor, in_ancestor_vector=metabolic_child_terms)
    
    return(BP_df)
  }

}

#############################################################
### A function for annotating the files for each cell type
#############################################################
Annotating_dir <- function(in_dir, out_fname) {
  ### Getting the files
  same_direction_files <- list.files(in_dir, 'same_direction', full.names = TRUE)
  oppo_direction_files <- list.files(in_dir, 'opposite_direction', full.names = TRUE)
  all_files <- c(same_direction_files, oppo_direction_files)

  ### Annotating each file and collating them
  rm(collated_df)
  for (curr_file in all_files) {
    print(curr_file)
    curr_df <- Annotating_df(curr_file)

    if (exists('collated_df') == FALSE) {
      collated_df <- curr_df
    } else {
      collated_df <- rbind(collated_df, curr_df)
    }
  
  }

  
  outfile <- paste(c(in_dir, out_fname), collapse = '')
  fwrite(collated_df, outfile, sep = '\t', row.names = FALSE, quote = FALSE)
}

Annotating_dir('./8_Statins_antidepressant_DEG_gProfiler2/HA1E_statin_antidepressant/', 'Statin_antidepressants_HA1E_sharedDEGS_gprofiler2_annotated.txt')

Annotating_dir('./8_Statins_antidepressant_DEG_gProfiler2/NPC_statin_antidepressant/', 'Statin_antidepressants_NPC_sharedDEGS_gprofiler2_annotated.txt')
