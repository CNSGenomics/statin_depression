import pandas as pd 

Master_dir = '/Users/uqjjian3/Statin_project/Statins_HA1E_CMAP_July27/'
CLUE_gct_dir =  Master_dir + '2_CMAP_query_Tau/CLUE_query_results/'
CLUE_txt_dir = Master_dir + '2_CMAP_query_Tau/Statins_HA1E_txt_files/'
original_gct_file = CLUE_gct_dir + 'Statins_24h_10uM_HA1E_v1_batch.gct'

Touchstone_cellline_list = ["A375", "A549", "HA1E", "HCC515", "HEPG2", "HT29", "MCF7", "PC3", "VCAP"]
HA1E_statin_list = ["atorvastatin", "fluvastatin", "lovastatin", "mevastatin", "pravastatin", "rosuvastatin", "simvastatin"]
statins_names_list = ["atorvastatin", "cerivastatin", "fluvastatin", "lovastatin", "mevastatin", "pitavastatin", "pravastatin", "rosuvastatin", "simvastatin"]
Psyc_drug_ids_TS = ['BRD-A05186015', 'BRD-A07106394', 'BRD-A10715913', 'BRD-A11990600', 'BRD-A14395271', 'BRD-A14966924', 'BRD-A16332958', 'BRD-A19195498', 'BRD-A19661776', 'BRD-A29644307', 'BRD-A31159102', 'BRD-A34309505', 'BRD-A43974499', 'BRD-A43974575', 'BRD-A44448661', 'BRD-A47598013', 'BRD-A51714012', 'BRD-A53077924', 'BRD-A59198242', 'BRD-A60197193', 'BRD-A62035778', 'BRD-A62428732', 'BRD-A64977602', 'BRD-A65280694', 'BRD-A80793822', 'BRD-A84481105', 'BRD-K00532621', 'BRD-K01292756', 'BRD-K02227374', 'BRD-K02265150', 'BRD-K02404261', 'BRD-K02867583', 'BRD-K03319035', 'BRD-K06980535', 'BRD-K07237224', 'BRD-K08619574', 'BRD-K10314788', 'BRD-K10995081', 'BRD-K12102668', 'BRD-K15409150', 'BRD-K16508793', 'BRD-K18779551', 'BRD-K18895904', 'BRD-K19352500', 'BRD-K19456237', 'BRD-K20141153', 'BRD-K26801045', 'BRD-K28761384', 'BRD-K29582115', 'BRD-K32398298', 'BRD-K35559145', 'BRD-K36116267', 'BRD-K36616567', 'BRD-K37289225', 'BRD-K37516142', 'BRD-K37688416', 'BRD-K37814297', 'BRD-K37991163', 'BRD-K38436528', 'BRD-K39915878', 'BRD-K42098891', 'BRD-K43786866', 'BRD-K44876623', 'BRD-K48722833', 'BRD-K50422030', 'BRD-K51671335', 'BRD-K52989797', 'BRD-K53318339', 'BRD-K53737926', 'BRD-K53857191', 'BRD-K54094468', 'BRD-K54759182', 'BRD-K55127134', 'BRD-K57432881', 'BRD-K57875380', 'BRD-K57930253', 'BRD-K59058766', 'BRD-K59273480', 'BRD-K59332007', 'BRD-K60762818', 'BRD-K63675182', 'BRD-K67783091', 'BRD-K68867920', 'BRD-K69116396', 'BRD-K70301876', 'BRD-K70358946', 'BRD-K70487031', 'BRD-K70778732', 'BRD-K70883034', 'BRD-K71103788', 'BRD-K72676686', 'BRD-K77947974', 'BRD-K78643075', 'BRD-K79145749', 'BRD-K79425933', 'BRD-K82036761', 'BRD-K82147103', 'BRD-K86595100', 'BRD-K87024524', 'BRD-K88611939', 'BRD-K89732114', 'BRD-K89997465', 'BRD-K90789829', 'BRD-K91263825', 'BRD-K92657060', 'BRD-K92984783', 'BRD-K93332168', 'BRD-K93461745', 'BRD-K97158071', 'BRD-K97530723'] ### NOTE: if a drug was profiled twice using different perturbagen ids in the dataset, the pert id found using /sig on the CLUE platform was used. 
Statin_pert_ids_TS = ["BRD-U88459701", "BRD-K81169441", "BRD-K66296774", "BRD-K09416995", "BRD-K94441233", "BRD-K30097969", "BRD-K60511616", "BRD-K82941592", "BRD-A81772229"]

###############################################################################
### 1. Getting the CMAP psyc drugs
###############################################################################
CMAP_psyc_drug_txt = "/Users/uqjjian3/Data/LINCS/GSE92742_Broad_LINCS_pert_ATC_psyc_drugs.txt"
CMAP_psyc_drug_list = []
ATC_code_dict = dict()
with open (CMAP_psyc_drug_txt, 'r') as a:
	lines = a.readlines()[1:]
	for line in lines:
		line = line.strip()
		fields = line.split('\t')
		if fields[0] in Psyc_drug_ids_TS:
			if fields[0] == "BRD-K51671335":
				ATC_code_dict[fields[0]] = "N05AL07" # /sig levosulpiride shows that it was mistakenly recorded as sulpiride
			else:
				ATC_code_dict[fields[0]] = fields[2]

print('Total number of CMAP psyc drugs', len(ATC_code_dict))

###############################################################################
### 2. Keeping only within cell line Tau Scores
###############################################################################
def getting_within_cellline_Tau(in_cell_line):
	print('Getting within cell line Tau')
	original_gct = pd.read_csv(original_gct_file, header=[3,4], sep='\t', dtype= str)
	original_gct = original_gct.dropna(how='all')
	same_cell_gct_df = original_gct[['cell_id', 'na', in_cell_line]]
	same_cell_gct_df.columns = same_cell_gct_df.columns.get_level_values(1)

	###renaming column names
	column_indices = [0,1,2,3,4,5,6,7,8]
	new_names = ['id', 'type', 'name', 'description', 'target', 'belongs_to', 'pc', 'median_tau_score', 'ts_pc']
	old_names = same_cell_gct_df.columns[column_indices]
	same_cell_gct_df.rename(columns=dict(zip(old_names, new_names)), inplace=True)
	same_cell_gct_df = same_cell_gct_df.drop(index = [0,1,2,3], columns = ['pc', 'belongs_to', 'median_tau_score', 'ts_pc'])

	to_strip_str = '_'+in_cell_line
	same_cell_gct_df.columns = same_cell_gct_df.columns.str.rstrip(to_strip_str)

	within_cell_txt_file = CLUE_txt_dir+"statins_"+in_cell_line+"_24h_within_cellline_Tau.txt"
	same_cell_gct_df.to_csv(within_cell_txt_file, sep = '\t', header = True, index = False)
	return

###############################################################################
### 3. Keeping only within cell line Tau Scores for psyc drugs
###############################################################################
def getting_psyc_drug_Tau(in_cell_line):
	print('Getting Tau with CMAP psyc drugs')
	within_cell_txt_file = CLUE_txt_dir+"statins_"+in_cell_line+"_24h_within_cellline_Tau.txt"
	within_cell_df = pd.read_csv(within_cell_txt_file, header = 0, sep = '\t')
	psyc_tau_df = within_cell_df[within_cell_df['id'].isin(Psyc_drug_ids_TS)]
	print(psyc_tau_df)
	psyc_tau_df = psyc_tau_df.dropna(subset = HA1E_statin_list)
	print(psyc_tau_df)
	psyc_tau_df['ATC_code'] = psyc_tau_df['id'].map(ATC_code_dict)
	psyc_tau_df.loc[psyc_tau_df.id == 'BRD-K51671335', 'name'] = "levosulpiride"
	print(psyc_tau_df.shape)

	psyc_tau_txt_file = CLUE_txt_dir+"statins_"+in_cell_line+"_24h_within_cellline_psyc_drug_Tau.txt"
	psyc_tau_df.to_csv(psyc_tau_txt_file, sep = '\t', header = True, index = False)
	return

###############################################################################
### 4. Keeping only within cell line Tau Scores for statins
###############################################################################
def getting_statin_drug_Tau(in_cell_line):
	print('Getting Tau with CMAP statin drugs')
	within_cell_txt_file = CLUE_txt_dir+"statins_"+in_cell_line+"_24h_within_cellline_Tau.txt"
	within_cell_df = pd.read_csv(within_cell_txt_file, header = 0, sep = '\t')
	statin_tau_df = within_cell_df[within_cell_df['id'].isin(Statin_pert_ids_TS)]
	profiled_statin_list = statin_tau_df.columns.values.tolist()[5:]
	statin_tau_df = statin_tau_df[statin_tau_df['name'].isin(profiled_statin_list)]
	statin_tau_df = statin_tau_df.sort_values(by = ['name'])

	statin_tau_txt_file = CLUE_txt_dir+"statins_"+in_cell_line+"_24h_within_cellline_statins_Tau.txt"
	statin_tau_df.to_csv(statin_tau_txt_file, sep = '\t', header = True, index = False)
	return

###############################################################################
### 5. Running functions
###############################################################################
TS_cellline = 'HA1E'
getting_within_cellline_Tau(TS_cellline)
getting_psyc_drug_Tau(TS_cellline)
getting_statin_drug_Tau(TS_cellline)

