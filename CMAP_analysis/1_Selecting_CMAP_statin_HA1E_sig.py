import os
import pandas as pd

########################################################################################################
### This script is to get the sig info for statins in the HA1E cell line for 10um and 24h treatment
Master_dir = '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/'
sig_info_dir = Master_dir + '1_GSE92742_sig/'
########################################################################################################
GSE92742_sig_info_file = '/Users/uqjjian3/Data/LINCS/GSE92742_Broad_LINCS_sig_info.txt' # Retrieved from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
statins_names_list = ["atorvastatin", "cerivastatin", "fluvastatin", "lovastatin", "mevastatin", "pitavastatin", "pravastatin", "rosuvastatin", "simvastatin"]
statins_IDs_list = ['BRD-U88459701', 'BRD-K81169441', 'BRD-K66296774', 'BRD-K09416995', 'BRD-K94441233', 'BRD-K30097969', 'BRD-K60511616', 'BRD-K82941592', 'BRD-A81772229']

GSE92742_sig_info_df = pd.read_csv(GSE92742_sig_info_file, sep = '\t', header = 0, dtype = str)
GSE92742_sig_info_df_cp = GSE92742_sig_info_df[GSE92742_sig_info_df['pert_type'] == 'trt_cp'] ### Filter for compound
GSE92742_sig_info_df_10uM = GSE92742_sig_info_df_cp[GSE92742_sig_info_df_cp['pert_idose'] == '10 µM'] ### Filter for 10um treatment
GSE92742_sig_info_df_24h = GSE92742_sig_info_df_10uM[GSE92742_sig_info_df_10uM['pert_time'] == '24'] ### Filter for 24h timepoint
GSE92742_sig_info_df_cell = GSE92742_sig_info_df_24h[GSE92742_sig_info_df_24h['cell_id'] == 'HA1E'] ### Filter for cell line
GSE92742_sig_info_df_statin = GSE92742_sig_info_df_cell[GSE92742_sig_info_df_cell['pert_iname'].isin(statins_names_list)] ### Filter for statin name
GSE92742_sig_info_df_statin_ID = GSE92742_sig_info_df_statin[GSE92742_sig_info_df_statin['pert_id'].isin(statins_IDs_list)] ### Some statins were profiled using different pert-id. Filter for the ID in touchstone
GSE92742_sig_info_df_statin_ID = GSE92742_sig_info_df_statin_ID.sort_values(by =['pert_iname'])
print(GSE92742_sig_info_df_statin_ID)

GSE92742_statin_24h_file = sig_info_dir + 'GSE92742_Broad_LINCS_sig_info_statin_24h_HA1E_TS.txt'
GSE92742_sig_info_df_statin_ID.to_csv(GSE92742_statin_24h_file, sep = '\t', index = False)
