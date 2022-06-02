###Tinaroo running condition

#PBS -l select=1:ncpus=10:mem=20GB
#PBS -l walltime=24:00:00
#PBS -N SMR_analysis
#PBS -e /home/uqjjian3/jobs_sterr_stout/
#PBS -o /home/uqjjian3/jobs_sterr_stout/
#PBS -A UQ-IMB-CNSG

import os

file_list = ['/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_baso_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_lymph_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_neut_p_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_rdw_cv_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_baso_p_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_lymph_p_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_ret_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_eo_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_mpv_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_ret_p_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_eo_p_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_mrv_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_plt_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_wbc_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_hlr_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_mono_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_pct_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_hlr_p_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_mono_p_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_pdw_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_irf_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_neut_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_rbc_Vuckovic_2020_N', '/QRISdata/Q4319/Statin_project_Clara/data/inflam_ahola_olli_2016/formatted/formatted_IL6_ahola_olli', '/QRISdata/Q4319/Statin_project_Clara/data/CRP_Han_GCST009777_gwascatalogue/formatted_CRP_Han_GCST009777', '/QRISdata/Q4319/Statin_project_Clara/data/GLGC_lipids_gwas/formatted_jointGwasMc_HDL', '/QRISdata/Q4319/Statin_project_Clara/data/GLGC_lipids_gwas/formatted_jointGwasMc_LDL', '/QRISdata/Q4319/Statin_project_Clara/data/GLGC_lipids_gwas/formatted_jointGwasMc_TG', '/QRISdata/Q4319/Statin_project_Clara/data/openGWAS/formatted_CAD_VDH_UKB_ebiaGCST005194', '/QRISdata/Q4319/Statin_project_Clara/data/openGWAS/formatted_T2D_Xue_ebiaGCST006867', '/QRISdata/Q4319/Statin_project_Clara/data/openGWAS/formatted_BMI_UKB_ukbb19953', '/QRISdata/Q4319/data/MDD_Howard_2018_no23andme/formatted_PGC_UKB_depression_no23andme', '/QRISdata/Q4319/Statin_project_Clara/data/Neuroticism_Nagel_2018_gwascatalogue/formatted_neuroticism_Nagel_2018', '/QRISdata/Q4319/Statin_project_Clara/data/Neuroticism_Nagel_2018_gwascatalogue/formatted_depressed_affect_Nagel_2018_NatGenet', '/QRISdata/Q4319/Statin_project_Clara/data/Neuroticism_Nagel_2018_gwascatalogue/formatted_worry_Nagel_2018_NatGenet'] 

### Command
smr_Linux_cmd = '/home/uqjjian3/utils/SMR/smr_Linux_YW'

### eQTL files
eQTLGEN_eQTL_file = '/QRISdata/Q1032/eQTLgen/cis_eQTL_SMR/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense'
Blood_GTEX_eQTL_file = '/QRISdata/Q3985/GTExV8/besd_hg19/Whole_Blood.v8.eqtl_signifpairs_hg19'
Brain_psychENCODE_eQTL_file = '/QRISdata/Q1032/PsychENCODE/Gandal_PsychENCODE_eQTL_HCP100+gPCs20_QTLtools.txt'

### gene list file
gene_list_hg19 = '/QRISdata/Q1032/Statins_Psychiatric/Chenwen_analysis/data/glist-hg19'

### LD reference files
in_5_bfile = '/QRISdata/Q1032/UKB_ref_LD/ukbEURu_imp_chr5_v3_impQC_10k_mac1'
in_16_bfile = '/QRISdata/Q1032/UKB_ref_LD/ukbEURu_imp_chr16_v3_impQC_10k_mac1'
in_6_bfile = '/QRISdata/Q1032/UKB_ref_LD/ukbEURu_imp_chr6_v3_impQC_10k_mac1'
in_1_bfile = '/QRISdata/Q1032/UKB_ref_LD/ukbEURu_imp_chr1_v3_impQC_10k_mac1'

def rerunning_SMR(in_gene, in_bfile, in_SNP_names, in_probe, in_cis_eQTL_file, in_label):
	print(in_gene)
	in_outfile_dir = '/home/uqjjian3/Statin_SMR_checkpoints_10052022/' + in_label + '/SMR_results/'
	os.makedirs(in_outfile_dir)

	for curr_trait_full_path in file_list:
		trait = curr_trait_full_path.split('/')[-1]
		print(trait)

		for SNP in in_SNP_names:
			print(SNP)
			curr_outfile = in_outfile_dir + in_gene + '_' + SNP + '_' + trait + '_peqtl_1_SMR'
			cmd = smr_Linux_cmd + ' --bfile '+in_bfile+ ' --gwas-summary '+curr_trait_full_path+ ' --beqtl-summary '+in_cis_eQTL_file +' --target-snp ' + SNP + ' --peqtl-smr 1 --out '+curr_outfile + ' --thread-num 10  --diff-freq 1'
			os.system(cmd)

			print('Generating plot files')
			curr_outfile = in_outfile_dir + in_gene + '_' + in_probe + '_' + SNP + '_' + trait + '_peqtl_1'
			cmd = smr_Linux_cmd + ' --bfile '+in_bfile+ ' --gwas-summary '+curr_trait_full_path+ ' --beqtl-summary '+in_cis_eQTL_file +' --target-snp ' + SNP + ' --peqtl-smr 1 --out ' + curr_outfile + ' --plot --probe ' + in_probe + ' --probe-wind 500  --diff-freq 1 --gene-list ' + gene_list_hg19  + ' --thread-num 10'
			os.system(cmd)

	return

rerunning_SMR('HMGCR', in_5_bfile, ['rs12916'], 'ENSG00000113161', eQTLGEN_eQTL_file, 'HMGCR_eQTLGEN')
rerunning_SMR('ITGAL', in_16_bfile, ['rs11574938'], 'ENSG00000005844', eQTLGEN_eQTL_file, 'ITGAL_eQTLGEN')
rerunning_SMR('HDAC2', in_6_bfile, ['rs9481408'], 'ENSG00000196591', eQTLGEN_eQTL_file, 'HDAC2_eQTLGEN')
rerunning_SMR('HMGCR', in_5_bfile, ['rs17671591'], 'ENSG00000113161', Brain_psychENCODE_eQTL_file, 'HMGCR_Brain_psychENCODE')
rerunning_SMR('PCSK9', in_1_bfile, ['rs12117661'], 'ENSG00000169174', Blood_GTEX_eQTL_file, 'PCSK9_Blood_GTEX')




