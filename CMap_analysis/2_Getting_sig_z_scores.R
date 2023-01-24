### This script is used to retrieve specific CMap gene expression signatures from drug-treated human cell lines, using signature IDs

library(cmapR)
library(dplyr)
library(stringr)
library(data.table)

z_score_dir <- './2_z_score_files/'

############################################################
### Getting the z-score signatures
############################################################
lincs_level5_gctx <- "./Data/LINCS/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"

getting_z_score <- function(in_sig_id_list, in_cp_sig_list, outfile) {
  ### Getting signatures
  lincs_level5_selected <- parse_gctx(lincs_level5_gctx, cid = in_sig_id_list, matrix_only=TRUE)
  lincs_level5_selected_df <- as.data.frame(mat(lincs_level5_selected))
  colnames(lincs_level5_selected_df) <- in_cp_sig_list
  
  ### Saving the file
  all_z_file <- paste(c(z_score_dir, outfile, '_all_genes.txt'), collapse = '')
  write.table(lincs_level5_selected_df, all_z_file, sep = '\t', col.names = NA, quote = FALSE)
  
  ### Getting the landmark genes only
  lincs_level5_selected_lm <-head (lincs_level5_selected_df, 978)
  lm_z_file <- paste(c(z_score_dir, outfile, '_lm_genes.txt'), collapse = '')
  write.table(lincs_level5_selected_lm, lm_z_file, sep = '\t', col.names=NA, quote = FALSE)
  
}

### Statins
getting_z_score(c('CPC004_HA1E_24H:BRD-K66296774-001-02-0:10', 'CPC003_HA1E_24H:BRD-K09416995-001-21-7:10', 'CPC001_HA1E_24H:BRD-K94441233-001-03-1:10', 'CPC002_HA1E_24H:BRD-K60511616-236-01-4:10', 'CPC004_HA1E_24H:BRD-K82941592-238-02-9:10', 'CPC005_HA1E_24H:BRD-A81772229-001-01-6:10', 'CPC006_HA1E_24H:BRD-U88459701-000-01-8:10'), c('fluvastatin_CPC004_HA1E_24H:BRD-K66296774-001-02-0:10', 'lovastatin_CPC003_HA1E_24H:BRD-K09416995-001-21-7:10', 'mevastatin_CPC001_HA1E_24H:BRD-K94441233-001-03-1:10', 'pravastatin_CPC002_HA1E_24H:BRD-K60511616-236-01-4:10', 'rosuvastatin_CPC004_HA1E_24H:BRD-K82941592-238-02-9:10', 'simvastatin_CPC005_HA1E_24H:BRD-A81772229-001-01-6:10', 'atorvastatin_CPC006_HA1E_24H:BRD-U88459701-000-01-8:10'), 'Statin_HA1E_10uM_24h')
getting_z_score(c('NMH001_NPC_24H:BRD-K69726342-001-02-6:10', 'CPC015_NPC_24H:BRD-K66296774-001-02-0:10', 'CPC015_NPC_24H:BRD-K09416995-001-21-7:10', 'CPC015_NPC_24H:BRD-K82941592-238-02-9:10', 'CPC016_NPC_24H:BRD-K94441233-001-03-1:10', 'CPC015_NPC_24H:BRD-K22134346-001-11-6:10'), c('atorvastatin_NMH001_NPC_24H:BRD-K69726342-001-02-6:10', 'fluvastatin_CPC015_NPC_24H:BRD-K66296774-001-02-0:10', 'lovastatin_CPC015_NPC_24H:BRD-K09416995-001-21-7:10', 'rosuvastatin_CPC015_NPC_24H:BRD-K82941592-238-02-9:10', 'mevastatin_CPC016_NPC_24H:BRD-K94441233-001-03-1:10', 'simvastatin_CPC015_NPC_24H:BRD-K22134346-001-11-6:10'), 'Statin_NPC_10uM_24h')

### Antidepressants
getting_z_score(c('CPC004_HA1E_24H:BRD-K60762818-003-15-3:10', 'CPC002_HA1E_24H:BRD-K91263825-001-03-6:10', 'CPC004_HA1E_24H:BRD-K37991163-003-06-8:10', 'CPC004_HA1E_24H:BRD-K82036761-003-07-0:10', 'CPC004_HA1E_24H:BRD-A19195498-050-09-1:10'), c('desipramine_CPC004_HA1E_24H:BRD-K60762818-003-15-3:10', 'nortriptyline_CPC002_HA1E_24H:BRD-K91263825-001-03-6:10', 'paroxetine_CPC004_HA1E_24H:BRD-K37991163-003-06-8:10', 'sertraline_CPC004_HA1E_24H:BRD-K82036761-003-07-0:10', 'trimipramine_CPC004_HA1E_24H:BRD-A19195498-050-09-1:10'), 'Antidepresssants_HA1E_10uM_24h')
getting_z_score(c('CPC017_NPC_24H:BRD-K60762818-003-15-3:10', 'CPC018_NPC_24H:BRD-K91263825-001-03-6:10', 'CPC015_NPC_24H:BRD-K37991163-003-06-8:10', 'CPC015_NPC_24H:BRD-K82036761-003-07-0:10', 'CPC016_NPC_24H:BRD-A19195498-050-09-1:10'), c('desipramine_CPC017_NPC_24H:BRD-K60762818-003-15-3:10', 'nortriptyline_CPC018_NPC_24H:BRD-K91263825-001-03-6:10', 'paroxetine_CPC015_NPC_24H:BRD-K37991163-003-06-8:10', 'sertraline_CPC015_NPC_24H:BRD-K82036761-003-07-0:10', 'trimipramine_CPC016_NPC_24H:BRD-A19195498-050-09-1:10'), 'Antidepressants_NPC_10uM_24h')

### Control drugs
getting_z_score(c('CPC012_HA1E_6H:BRD-K89626439-001-01-0:10', 'CPC014_HA1E_6H:BRD-K83988098-003-01-8:10'), c('sirolimus_CPC012_HA1E_6H:BRD-K89626439-001-01-0:10', 'alvespimycin_CPC014_HA1E_6H:BRD-K83988098-003-01-8:10'), 'Control_HA1E_10uM_6h')
