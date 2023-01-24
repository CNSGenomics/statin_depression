### This script is used to run SMR analysis, to investigate the causal association between genetically predicted exposure (HMGCR, ITGAL, HDAC2 and PCSK9 inhibition), and various outcomes 

import os

file_list = ['~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_baso_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_lymph_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_neut_p_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_rdw_cv_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_baso_p_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_lymph_p_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_ret_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_eo_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_mpv_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_ret_p_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_eo_p_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_mrv_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_plt_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_wbc_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_hlr_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_mono_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_pct_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_hlr_p_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_mono_p_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_pdw_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_irf_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_neut_Vuckovic_2020_N', '~/data/Blood_Cell_Vuckovic_2020/gcta_formatted/formatted_rbc_Vuckovic_2020_N', '~/data/inflam_ahola_olli_2016/formatted/formatted_IL6_ahola_olli', '~/data/CRP_Han_GCST009777_gwascatalogue/fomatted_CRP_Han_GCST009777', '~/data/GLGC_lipids_gwas/formatted_jointGwasMc_HDL', '~/data/GLGC_lipids_gwas/formatted_jointGwasMc_LDL', '~/data/GLGC_lipids_gwas/formatted_jointGwasMc_TG', '~/data/openGWAS/formatted_CAD_VDH_UKB_ebiaGCST005194', '~/data/openGWAS/formatted_T2D_Xue_ebiaGCST006867', '~/data/openGWAS/formatted_BMI_UKB_ukbb19953', '~/data/MDD_Howard_2018_no23andme/formatted_PGC_UKB_depression_no23andme', '~/data/Neuroticism_Nagel_2018_gwascatalogue/formatted_neuroticism_Nagel_2018', '~/data/Neuroticism_Nagel_2018_gwascatalogue/formatted_depressed_affect_Nagel_2018_NatGenet', '~/data/Neuroticism_Nagel_2018_gwascatalogue/formatted_worry_Nagel_2018_NatGenet'] 

### Command
smr_Linux_cmd = '~/utils/SMR/smr_Linux_YW'

### eQTL files
eQTLGEN_eQTL_file = '~/eQTLgen/cis_eQTL_SMR/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense'
Blood_GTEX_eQTL_file = '~/GTExV8/besd_hg19/Whole_Blood.v8.eqtl_signifpairs_hg19'
Brain_psychENCODE_eQTL_file = '~/PsychENCODE/Gandal_PsychENCODE_eQTL_HCP100+gPCs20_QTLtools.txt'

### gene list file
gene_list_hg19 = '~/data/glist-hg19'

### LD reference files
in_5_bfile = '~/UKB_ref_LD/ukbEURu_imp_chr5_v3_impQC_10k_mac1'
in_16_bfile = '~/UKB_ref_LD/ukbEURu_imp_chr16_v3_impQC_10k_mac1'
in_6_bfile = '~/UKB_ref_LD/ukbEURu_imp_chr6_v3_impQC_10k_mac1'
in_1_bfile = '~/UKB_10K_LD/ukbEURu_imp_chr1_v3_impQC_10k_mac1'

def rerunning_SMR(in_gene, in_bfile, in_SNP_names, in_probe, in_cis_eQTL_file, in_label):
	print(in_gene)
	in_outfile_dir = './' + in_label + '/'
	os.makedirs(in_outfile_dir)

	for curr_trait_full_path in file_list:
		trait = curr_trait_full_path.split('/')[-1]
		print(trait)

		for SNP in in_SNP_names:
			print(SNP)
			curr_outfile = in_outfile_dir + in_gene + '_' + SNP + '_' + trait + '_peqtl_1_SMR'
			cmd = smr_Linux_cmd + ' --bfile '+in_bfile+ ' --gwas-summary '+curr_trait_full_path+ ' --beqtl-summary '+in_cis_eQTL_file +' --target-snp ' + SNP + ' --peqtl-smr 1 --out '+curr_outfile + ' --thread-num 20  --diff-freq 1'
			os.system(cmd)

	return

rerunning_SMR('HMGCR', in_5_bfile, ['rs12916'], 'ENSG00000113161', eQTLGEN_eQTL_file, 'HMGCR_eQTLGEN')
rerunning_SMR('ITGAL', in_16_bfile, ['rs11574938'], 'ENSG00000005844', eQTLGEN_eQTL_file, 'ITGAL_eQTLGEN')
rerunning_SMR('HDAC2', in_6_bfile, ['rs9481408'], 'ENSG00000196591', eQTLGEN_eQTL_file, 'HDAC2_eQTLGEN')
rerunning_SMR('HMGCR', in_5_bfile, ['rs17671591'], 'ENSG00000113161', Brain_psychENCODE_eQTL_file, 'HMGCR_Brain_psychENCODE')
rerunning_SMR('PCSK9', in_1_bfile, ['rs12117661'], 'ENSG00000169174', Blood_GTEX_eQTL_file, 'PCSK9_Blood_GTEX')




