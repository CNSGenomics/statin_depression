library(cmapR)
library(dplyr)


###########################################################################
### 1. Getting the 24h signature IDs for Statins
###########################################################################
statin_24h_file <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/1_GSE92742_sig/GSE92742_Broad_LINCS_sig_info_statin_24h_HA1E_TS.txt'
statin_24h_df <- read.csv(statin_24h_file, sep = '\t')
statin_24h_df <- statin_24h_df %>%
  rowwise() %>%      # for each row
  mutate(new_id = paste(c(pert_iname, cell_id), collapse = "_")) %>%
  ungroup()

statin_24h_sig_ID_list <- statin_24h_df$sig_id
statin_24h_sig_ID_list_new <- statin_24h_df$new_id

###########################################################################
### 2. Getting the 24h signatures
###########################################################################
lincs_level5_gctx <- "/Users/uqjjian3/Data/LINCS/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx" # Retrieved from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
sig_dir <- '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/1_GSE92742_sig/'
lincs_level5_statin <- parse_gctx(lincs_level5_gctx, cid = statin_24h_sig_ID_list, matrix_only=TRUE)
lincs_level5_statin_zscores_m <- mat(lincs_level5_statin)
colnames(lincs_level5_statin_zscores_m) <- statin_24h_sig_ID_list_new
lincs_level5_statin_zscores_df <- as.data.frame(lincs_level5_statin_zscores_m)

all_z_file <- paste(c(sig_dir, 'GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_Statins_HA1E_all_genes.txt'), collapse = '')
write.table(lincs_level5_statin_zscores_df, all_z_file, sep = '\t', col.names=NA)

lincs_level5_statin_zscores_m_lm <-head(lincs_level5_statin_zscores_m, 978)
lincs_level5_statin_zscores_df_lm <- as.data.frame(lincs_level5_statin_zscores_m_lm)
lm_z_file <- paste(c(sig_dir, 'GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_Statins_HA1E_lm_genes.txt'), collapse = '')
write.table(lincs_level5_statin_zscores_df_lm, lm_z_file, sep = '\t', col.names=NA)

###########################################################################
### 3. Generating query signatures
###########################################################################
Getting_query <- function(in_full_sig_df, in_statin_name, in_cell_name, number_of_top_genes){
  in_sig_name <- paste(c(in_statin_name, in_cell_name), collapse = "_")
  print(in_sig_name)
  if(in_sig_name %in% colnames(in_full_sig_df))  {
    current_df <- select(in_full_sig_df, {{in_sig_name}})
    print(dim(current_df))
    
    colnames(current_df) <- c("Zval")
    current_query_up<-Getting_top_upreg(current_df, 50)
    print(current_query_up)
    current_query_down<-Getting_top_downreg(current_df, 50)
    print(current_query_down)
    
    current_up_gmt_list <- list('head' = in_sig_name, 'desc' = in_sig_name, 'len' = "50" , 'entry' = current_query_up)
    current_up_gmt_list_out <- list(in_sig_name = current_up_gmt_list)
    print(current_up_gmt_list_out)
    up_gmt_file <- paste(c(gmt_dir,in_sig_name,'_24h_10uM_up_query.gmt'), collapse = "")
    write_gmt(current_up_gmt_list_out, up_gmt_file)
    
    current_down_gmt_list <- list('head' = in_sig_name, 'desc' = in_sig_name, 'len' = "50" , 'entry' = current_query_down)
    current_down_gmt_list_out <- list(in_sig_name = current_down_gmt_list)
    print(current_down_gmt_list_out)
    down_gmt_file <- paste(c(gmt_dir,in_sig_name,'_24h_10uM_down_query.gmt'), collapse = "")
    write_gmt(current_down_gmt_list_out, down_gmt_file)
  } else {
    print("signature not found")
  }
}
Getting_top_upreg <- function(in_fil_df, number_of_top_genes) {
  Upreg_fil_df <- filter(in_fil_df, Zval > 0)
  print(paste(c('Total number of Upreg gene', length(row.names(Upreg_fil_df))), collapse = ' '))
  Upreg_fil_df_sorted <- Upreg_fil_df %>% arrange(desc(Zval))
  Top_upreg_genes = rownames(Upreg_fil_df_sorted)[1:number_of_top_genes] 
  return(c(Top_upreg_genes))
}
Getting_top_downreg <- function(in_fil_df, number_of_top_genes) {
  Downreg_fil_df <- filter(in_fil_df, Zval < 0)
  print(paste(c('Total number of Downreg gene', length(row.names(Downreg_fil_df))), collapse = ' '))
  Downreg_fil_df_sorted <- arrange(Downreg_fil_df,Zval)
  Top_downreg_genes = rownames(Downreg_fil_df_sorted)[1:number_of_top_genes]
  return(c(Top_downreg_genes))
}  

gmt_dir <-  '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/2_CMAP_query_Tau/Statins_HA1E_gmt_files/'
statins_names_list <- c("atorvastatin", "cerivastatin", "fluvastatin", "lovastatin", "mevastatin", "pitavastatin", "pravastatin", "rosuvastatin", "simvastatin")
Cell_id_TS <- c("HA1E")
for (query_statin in statins_names_list){
  for (Cell in Cell_id_TS) {
    Getting_query(lincs_level5_statin_zscores_df_lm, query_statin, Cell, 50)
  }
}

###########################################################################
### 4. Combining query for the online CLUE tool
###########################################################################
statins_HA1E_names_list <- c("atorvastatin", "fluvastatin", "lovastatin", "mevastatin", "pravastatin", "rosuvastatin", "simvastatin")
gmt_files_list <- list.files(gmt_dir)

for (Cell_line in Cell_id_TS){
  uplist <- list()
  downlist <- list()
  for (statin_name in statins_HA1E_names_list){
    sig_name <- paste(c(statin_name, "_", Cell_line), collapse = "")
    up_gmt_file <- paste(c(sig_name,"_24h_10uM_up_query.gmt"), collapse = "")
    down_gmt_file <- paste(c(sig_name, "_24h_10uM_down_query.gmt"), collapse = "")
    if (up_gmt_file %in% gmt_files_list){
      up_gmt_file <- paste(c(gmt_dir, up_gmt_file), collapse = "")
      down_gmt_file <- paste(c(gmt_dir, down_gmt_file), collapse = "")
      up_gmt <- parse_gmt(up_gmt_file)
      down_gmt <- parse_gmt(down_gmt_file)
      
      uplist <- c(uplist, sig_name = up_gmt)
      downlist <- c(downlist, sig_name = down_gmt)
    }
  }
  print(Cell_line)
  up_cell_line_gmt_file <- paste(c(gmt_dir,"statins_",Cell_line, "_24h_10uM_up_query_50.gmt"), collapse = "")
  write_gmt(uplist, up_cell_line_gmt_file)
  print(length(uplist))
  down_cell_line_gmt_file <- paste(c(gmt_dir,"statins_",Cell_line, "_24h_10uM_down_query_50.gmt"), collapse = "")
  write_gmt(downlist, down_cell_line_gmt_file)
  print(length(downlist))
}

###########################################################################
### 5. Submitted to only CLUE query, compared against v1.0 without the fastgutc tool
###########################################################################
### Gene expression (L1000)
### Touchstone
### Batch query
### 1.0
### Submission name: Statins_24h_10uM_HA1E_v1_batch
### Submission time: July 27th 2021 5:28pm