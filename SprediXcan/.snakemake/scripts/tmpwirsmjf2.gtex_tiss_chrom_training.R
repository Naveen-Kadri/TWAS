
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('/cluster/work/pausch/naveen/TWAS/GENO/CHR24/snp_annotation.txt', '/cluster/work/pausch/naveen/TWAS/GENO/CHR24/genotypes.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/annotation.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/covariates.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/mpheno.txt', "snp_annotation" = '/cluster/work/pausch/naveen/TWAS/GENO/CHR24/snp_annotation.txt', "genotypes" = '/cluster/work/pausch/naveen/TWAS/GENO/CHR24/genotypes.txt', "gene_annotation" = '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/annotation.txt', "covariates" = '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/covariates.txt', "mpheno" = '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/mpheno.txt'),
    output = list('/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/covariance.txt', "model_summary_file" = '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/summary.txt', "weights_file" = '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/weights.txt', "tiss_chr_summ_f" = '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/tiss_summ.txt', "covariance_file" = '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/covariance.txt'),
    params = list(),
    wildcards = list('splicing', 'vas_d', '24', "mpheno" = 'splicing', "tissue" = 'vas_d', "chr" = '24'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', 'walltime', "mem_mb" = 6000, "mem_mib" = 5723, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/scratch/257441532.tmpdir', "walltime" = '48:00'),
    config = list("molecular_pheno" = list("expression" = '/cluster/work/pausch/naveen/TWAS/TPM/{tissue}.tsv', "splicing" = '/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/{tissue}.bed.gz'), "covariates" = list("expression" = '/cluster/work/pausch/xena/eQTL/cis_all/{tissue}/{tissue}_covariates.txt', "splicing" = '/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/covariates_{tissue}.txt'), "vcfs" = list("twas" = '/cluster/work/pausch/naveen/ASE/VCF/chr_{chr}_beagle5.2.vcf.gz', "gwas" = '/cluster/work/pausch/naveen/SNPDATA/TOSEQ/CHR{chr}/combined.vcf.gz'), "gtf" = '/cluster/work/pausch/naveen/ASE/INPUTS/Bos_taurus.ARS-UCD1.2.104.gtf.gz', "OUT_DIR" = '/cluster/work/pausch/naveen/TWAS', "old_gwas_files" = '/cluster/work/pausch/xena/gwas/gcta/results/{{trait}}/BV/sire/CHR{chr}/result.mlma', "gwas_files" = '/cluster/work/pausch/xena/gwas/new/{inheri}/{{trait}}/CHR{chr}/result.mlma.gz'),
    rule = 'train',
    bench_iteration = as.numeric(NA),
    scriptdir = '/cluster/home/nkadri/TWAS/SNAKES/train',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
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


