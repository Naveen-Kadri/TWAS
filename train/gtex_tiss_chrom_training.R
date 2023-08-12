##source("gtex_v7_nested_cv_elnet.R")
"%&%" <- function(a,b) paste(a,b, sep='')

chrom <- snakemake@wildcards[["chr"]]
snp_annot_file <- snakemake@input [["snp_annotation"]]
gene_annot_file <-snakemake@input [["gene_annotation"]]
genotype_file <- snakemake@input [["genotypes"]]
mpheno_file <-snakemake@input [["mpheno"]]
covariates_file <-snakemake@input [["covariates"]]

model_summary_file <- snakemake@output[["model_summary_file"]]
weights_file <- snakemake@output[["weights_file"]]
tiss_chr_summ_f<-snakemake@output[["tiss_chr_summ_f"]]
covariance_file<-snakemake@output[["covariance_file"]]


## chrom <- 25
## snp_annot_file <- '/cluster/work/pausch/naveen/TWAS/GENO/CHR25/snp_annotation.txt'
## gene_annot_file <- '/cluster/work/pausch/naveen/TWAS/splicing/testis/annotation.txt'
## genotype_file <- '/cluster/work/pausch/naveen/TWAS/GENO/CHR25/genotypes.txt'
## mpheno_file <- '/cluster/work/pausch/naveen/TWAS/splicing/testis/mpheno.txt'
## covariates_file <- '/cluster/work/pausch/naveen/TWAS/splicing/testis/covariates.txt'

## model_summary_file <- 'model_summary.txt'
## weights_file <- 'wts.txt'
## tiss_chr_summ_f<-'tsschr.txt'
## covariance_file<-'cov.txt'



## prefix removed
.libPaths('/cluster/work/pausch/naveen/Rlib/')
source ('~/TWAS/SNAKES/train/gtex_v7_nested_cv_elnet.R')

main(snp_annot_file, gene_annot_file, genotype_file, mpheno_file, covariates_file, as.numeric(chrom), null_testing=FALSE)
##main(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, as.numeric(chrom), prefix, null_testing=FALSE)


