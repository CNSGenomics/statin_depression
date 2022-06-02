import pandas as pd

Master_dir = '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/'
sig_info_dir = Master_dir + '1_GSE92742_sig/'


all_z_file = sig_info_dir + 'GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_Statins_HA1E_all_genes.txt'

Statins_HA1E_list = ["atorvastatin", "fluvastatin", "lovastatin", "mevastatin", "pravastatin", "rosuvastatin", "simvastatin"]

###########################################################################
### 1. Preparing dataframes
###########################################################################
def preparing_df(in_file):
	z_scores_df = pd.read_csv(in_file, sep = '\t', index_col = 0)
	return(z_scores_df)

# ###########################################################################
# ### 2. Generating differentially expressed genes (DEGs)
# ###########################################################################
def Getting_top_upreg(in_fil_df, in_threshold):
	Upreg_fil_df = in_fil_df[in_fil_df['Zval'] > in_threshold]
	print('Total number of Up DE gene', len(Upreg_fil_df.index))
	upreg_fil_df_sorted = Upreg_fil_df.sort_values(by = ['Zval'], ascending = False)
	DEG_upreg_genes = upreg_fil_df_sorted.index.tolist()
	return(DEG_upreg_genes)

def Getting_top_downreg(in_fil_df, in_threshold):
	Downreg_fil_df = in_fil_df[in_fil_df['Zval'] < (-in_threshold)]
	print('Total number of Downreg DE gene', len(Downreg_fil_df.index))
	Downreg_fil_df_sorted = Downreg_fil_df.sort_values(by = ['Zval'])
	DEG_upreg_genes = Downreg_fil_df_sorted.index.tolist()
	return(DEG_upreg_genes)

def Getting_DEG(in_full_sig_df, in_statin_name, in_cell_name, out_deg_dir,in_threshold):
	in_sig_name = in_statin_name+'_'+in_cell_name
	current_df = in_full_sig_df[[in_sig_name]]
	current_df.columns = ['Zval']

	current_DEG_up = Getting_top_upreg(current_df, in_threshold)
	current_DEG_down = Getting_top_downreg(current_df, in_threshold)

	up_DEG_file = DEG_dir + out_deg_dir + '/' + in_sig_name + '_24h_' + out_deg_dir +'_up.txt'
	with open (up_DEG_file, 'w') as uo:
		for up in current_DEG_up:
			uo.write(str(up) + '\n')

	down_DEG_file = DEG_dir + out_deg_dir + '/' + in_sig_name + '_24h_' + out_deg_dir +'_down.txt'
	with open (down_DEG_file, 'w') as do:
		for down in current_DEG_down:
			do.write(str(down) + '\n')

all_z_df = preparing_df(all_z_file)
for statin_name in Statins_HA1E_list:
	print(statin_name)
	DEG_dir = Master_dir + '4_HA1E_24_DEGs_thresh1/'
	Getting_DEG(all_z_df, statin_name, 'HA1E', 'DEG_thresh1_all_genes', 1)
	DEG_dir = Master_dir + '4_HA1E_24_DEGs_thresh1.5/'
	Getting_DEG(all_z_df, statin_name, 'HA1E', 'DEG_thresh1.5_all_genes', 1.5)
	DEG_dir = Master_dir + '4_HA1E_24_DEGs_thresh2/'
	Getting_DEG(all_z_df, statin_name, 'HA1E', 'DEG_thresh2_all_genes', 2)