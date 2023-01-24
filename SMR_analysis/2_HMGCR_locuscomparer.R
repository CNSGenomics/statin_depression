### The script is used to plot the locusCompare plots between HMGCR eQTL data and the LDL-C GWAS summary statistics

library(locuscomparer)
library(ggplot2)
library(stringr)

### The summary statistics for LDL-C
GWAS_file <- '~/data/GLGC_lipids_gwas/formatted_jointGwasMc_LDL'
GWAS_title <- 'LDL-C GWAS'

### The HMGCR eQTL file
eQTL_file <- '~/data/eQTLGEN_summary/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense_p_1_HMGCR.txt'
eQTL_title <- 'HMGCR eQTL (eQTLGen)'

### Plotting the locusCompare plot
out_pdf <- './LDLC_HMGCR_eQTL_locusCompare_rs12916.pdf'
curr_fig <-locuscompare(in_fn1=GWAS_file, in_fn2=eQTL_file, title1=GWAS_title, title2=eQTL_title, marker_col1= 'SNP', pval_col1='p', marker_col2='SNP', pval_col2='p', snp = 'rs12916')
ggsave(out_pdf, plot = curr_fig, device = NULL, path = NULL, scale = 1, width = 200, height = 180, units = c('mm'), dpi = 300, limitsize = TRUE, bg = NULL) 

