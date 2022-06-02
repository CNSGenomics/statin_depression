
ATC_index = "/Users/uqjjian3/Data/ATC_drugs/ATC_classification_unofficial.txt" # the file was downloaded from https://bioportal.bioontology.org/ontologies/ATC. Version 2021AA
CMAP_GSE92742_pert_metrics_txt = "/Users/uqjjian3/Data/LINCS/GSE92742_Broad_LINCS_pert_metrics.txt" # Retrieved from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
Anno_outdir = "/Users/uqjjian3/Statin_project/CMAP_statin/Annotating_psychotic_drugs/"
Drugs_ID_list = ['N05A', 'N05B', 'N05C', 'N06A', 'N06B']

############################################################################################################
print('Getting ATC psyc drugs')
Drug_line_dict = dict()
Drug_line_list = []
Drug_name_list = []
with open (ATC_index, 'r') as a:
	lines = a.readlines()
	for line in lines:
		line = line.strip()
		fields = line.split('\t')
		ATC_code = fields[0]
		drug_name = fields[1]
		if len(ATC_code) == 7:
			if ATC_code[:4] in Drugs_ID_list:
				# print(drug_name)
				Drug_line_list.append(line)
				Drug_line_dict[drug_name] = ATC_code
				Drug_name_list.append(drug_name)

Drug_line_list = list(set(Drug_line_list))
print('Total number of ATC psyc drugs', len(Drug_line_list))

Psyc_drug_txt = "/Users/uqjjian3/Data/ATC_drugs/ATC_classification_unofficial_psyc_drugs.txt"
with open (Psyc_drug_txt, 'w') as ao:
	for line in Drug_line_list:
		ao.write(line +'\n')

############################################################################################################
print('Filtering CMAP perturbagens for ATC psyc drugs')
CMAP_psyc_drug_line_list = []
with open (CMAP_GSE92742_pert_metrics_txt, 'r') as b:
	lines = b.readlines()
	for line in lines[1:]:
		line = line.strip()
		fields = line.split('\t')

		pert_id = fields[0]
		pert_iname = fields[1]
		pert_type = fields[2]

		if pert_type == "trt_cp":
			if pert_iname in Drug_name_list:
				# print(pert_iname)
				ATC_code_current = Drug_line_dict[pert_iname]
				newline = '\t'.join([pert_id, pert_iname, ATC_code_current])
				CMAP_psyc_drug_line_list.append(newline)
print('Total number of ATC psyc drugs profiled by CMAP', len(CMAP_psyc_drug_line_list))

CMAP_psyc_drug_line_list = list(set(CMAP_psyc_drug_line_list))
CMAP_psyc_drug_txt = "/Users/uqjjian3/Data/LINCS/GSE92742_Broad_LINCS_pert_ATC_psyc_drugs.txt"
with open (CMAP_psyc_drug_txt, 'w') as bo:
	header = '\t'.join(['pert_id', 'pert_iname', 'ATC_code'])
	bo.write(header +'\n')
	for line in CMAP_psyc_drug_line_list:
		bo.write(line +'\n')



