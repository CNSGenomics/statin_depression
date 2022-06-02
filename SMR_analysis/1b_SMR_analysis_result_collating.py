import os
import math
import numpy as np

HDAC2_eQTLGEN_dir = '/home/uqjjian3/Statin_SMR_checkpoints_10052022/HDAC2_eQTLGEN/SMR_results/'
HMGCR_Brain_psychENCODE_dir = '/home/uqjjian3/Statin_SMR_checkpoints_10052022/HMGCR_Brain_psychENCODE/SMR_results/'
HMGCR_eQTLGEN_dir = '/home/uqjjian3/Statin_SMR_checkpoints_10052022/HMGCR_eQTLGEN/SMR_results/'
ITGAL_eQTLGEN_dir = '/home/uqjjian3/Statin_SMR_checkpoints_10052022/ITGAL_eQTLGEN/SMR_results/'
PCSK9_Blood_GTEX_dir = '/home/uqjjian3/Statin_SMR_checkpoints_10052022/PCSK9_Blood_GTEX/SMR_results/'

def collating_results(SMR_dir, in_gene, in_gene_ID, in_SNP):
	print(in_gene, in_SNP)

	gene_SNP_ID = in_gene + '_' + in_SNP + '_'
	in_label = SMR_dir.split('/')[4]

	smr_result_list = []
	header_line = '\t'.join(['trait', 'probeID', 'ProbeChr', 'Gene', 'Probe_bp', 'targetSNP', 'targetSNP_chr', 'targetSNP_bp', 'A1', 'A2', 'Freq', 'b_GWAS', 'se_GWAS', 'p_GWAS', 'b_eQTL', 'se_eQTL', 'p_eQTL', 'b_SMR', 'se_SMR', 'p_SMR', 'p_HEIDI', 'nsnp_HEIDI', 'b_SMR_perSDdecrease'])
	smr_result_list.append(header_line)

	for root, dirs, files in os.walk(SMR_dir):
		for file in files:
			if in_SNP in file and file.endswith('.smr'):
				full_path = os.path.join(root, file)

				trait = file.split(gene_SNP_ID)[-1]
				trait = trait.split('_peqtl_1')[0]
				trait = trait.replace('formatted_', '')

				with open (full_path, 'r') as a:
					lines = a.readlines()[1:]
					for line in lines:
						line = line.strip()
						fields = line.split()

						curr_gene_ID = fields[0]
						curr_SNP = fields[4]
						curr_bSMR = fields[16]
						b_SMR_perSDdecrease = str(-(float(curr_bSMR)))

						if curr_gene_ID.startswith(in_gene_ID):
							newline = '\t'.join([trait, line, b_SMR_perSDdecrease])
							smr_result_list.append(newline)


	out_collated_file = '/home/uqjjian3/Statin_SMR_checkpoints_10052022/' + in_gene + '_' + in_SNP + '_' + in_label + '_SMR_Result_collated.txt'
	with open (out_collated_file, 'w') as a:
		for line in smr_result_list:
			a.write(line + '\n')

	return


collating_results(HMGCR_eQTLGEN_dir, 'HMGCR', 'ENSG00000113161', 'rs12916')
collating_results(ITGAL_eQTLGEN_dir, 'ITGAL', 'ENSG00000005844', 'rs11574938')
collating_results(HDAC2_eQTLGEN_dir, 'HDAC2', 'ENSG00000196591', 'rs9481408')
collating_results(HMGCR_Brain_psychENCODE_dir, 'HMGCR', 'ENSG00000113161', 'rs17671591')
collating_results(PCSK9_Blood_GTEX_dir, 'PCSK9', 'ENSG00000169174', 'rs12117661')




