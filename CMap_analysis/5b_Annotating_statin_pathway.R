### This script is used to annotate significantly enriched biological processes with their ancestors. As the results of gProfiler2 often contain redundant and overlapping GO:BP terms, we have grouped these GO:BP terms based on their ancestor terms.

library(GO.db) # https://bioconductor.org/packages/release/data/annotation/manuals/GO.db/man/GO.db.pdf
library(data.table)
library(ggplot2)
library(GSA)
library(stringr)

#############################################################
### Getting list of pathways associated with each ancestor pathway
#############################################################
getting_child_terms <- function(in_ancestor_term) {
  child_term_df  <- as.data.frame(GOBPCHILDREN[c(in_ancestor_term)])
  child_term_df <- child_term_df[child_term_df$RelationshipType == 'isa', ]
  child_term <- child_term_df$go_id
  return(child_term)
}

root_child_terms <- getting_child_terms('GO:0008150') # the master term for "biological processes"
immune_child_terms <- getting_child_terms('GO:0002376') # the master term for immune processes
metabolic_child_terms <- getting_child_terms('GO:0044238') # the master term for primary metabolic process

#############################################################
### A function to check and annotate the ancestors
#############################################################
annotating_ancestor <- function(in_BP, in_ancestor_vector) {
  if (in_BP == 'GO:0055114') {
    # It was found that "GO:0055114", which is an obsolete term, were not in the database, so it was manually annotated to a "metabolic process", based on its description on https://www.ebi.ac.uk/QuickGO/term/GO:0055114#:~:text=Definition%20(GO%3A0055114%20GONUTS%20page)&text=A%20metabolic%20process%20that%20results,of%20a%20proton%20or%20protons. 
    curr_BP_ancestor_vector <- c('GO:0008152')
  } else {
    ### Getting the ancestors of the current term
    curr_BP_ancestor_vector <- as.list(GOBPANCESTOR[c(in_BP)])[[1]]
    ### Adding the current term to the vector in case that it is one of the ancestor terms
    curr_BP_ancestor_vector <- c(curr_BP_ancestor_vector, in_BP)
  }
  ### Now looking at whether its ancestor terms overlap with the ancestor terms of interest
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
  if (length(colnames(in_df)) > 1){ # to confirm that file is not empty
    
    ### Filtering for GO:BP terms
    BP_df <- in_df[in_df$source == 'GO:BP', ]
    
    ### Annotating the statin type, cell type, the threshold for defining DEGs and the direction of DEGs
    fname <- str_split(in_file, '/')[[1]][length(str_split(in_file, '/')[[1]])]
    fname <- str_split(fname, '.txt')[[1]][1]
    curr_statin <- str_split(fname, '_')[[1]][1]
    curr_cell <- str_split(fname, '_')[[1]][2]
    curr_direction <- str_split(fname, '_gProfiler2_')[[1]][length(str_split(fname, '_gProfiler2_')[[1]])]
    curr_thresh <- str_split(fname, '_')[[1]][3]
    
    BP_df$Statin <- curr_statin
    BP_df$Cell <- curr_cell
    BP_df$Direction_DEGs <- curr_direction
    BP_df$thresh_DEGs <- curr_thresh
    
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
  up_direction_files <- list.files(in_dir, 'up.txt', full.names = TRUE)
  down_direction_files <- list.files(in_dir, 'down.txt', full.names = TRUE)
  all_files <- c(up_direction_files, down_direction_files)

  ### Annotating each file and collating them
  if (exists('collated_df')) {
    rm(collated_df)  
  }
  
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

Annotating_dir('./5_Statin_DEG_gProfiler/HA1E_statins/', 'statins_HA1E_DEGs_gprofiler2_annotated.txt')
Annotating_dir('./5_Statin_DEG_gProfiler/NPC_statins/', 'statins_NPC_DEGs_gprofiler2_annotated.txt')
