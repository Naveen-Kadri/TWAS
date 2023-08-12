
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
    input = list('/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/covariates_vas_d.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR25/mpheno.txt', "mpheno" = '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR25/mpheno.txt'),
    output = list('/cluster/work/pausch/naveen/TWAS/splicing/vas_d/covariates.txt', "outfile" = '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/covariates.txt'),
    params = list(),
    wildcards = list('splicing', 'vas_d', "mpheno" = 'splicing', "tissue" = 'vas_d'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', 'walltime', "mem_mb" = 2000, "mem_mib" = 1908, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/scratch/257441117.tmpdir', "walltime" = '00:20'),
    config = list("molecular_pheno" = list("expression" = '/cluster/work/pausch/naveen/TWAS/TPM/{tissue}.tsv', "splicing" = '/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/{tissue}.bed.gz'), "covariates" = list("expression" = '/cluster/work/pausch/xena/eQTL/cis_all/{tissue}/{tissue}_covariates.txt', "splicing" = '/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/covariates_{tissue}.txt'), "vcfs" = list("twas" = '/cluster/work/pausch/naveen/ASE/VCF/chr_{chr}_beagle5.2.vcf.gz', "gwas" = '/cluster/work/pausch/naveen/SNPDATA/TOSEQ/CHR{chr}/combined.vcf.gz'), "gtf" = '/cluster/work/pausch/naveen/ASE/INPUTS/Bos_taurus.ARS-UCD1.2.104.gtf.gz', "OUT_DIR" = '/cluster/work/pausch/naveen/TWAS', "old_gwas_files" = '/cluster/work/pausch/xena/gwas/gcta/results/{{trait}}/BV/sire/CHR{chr}/result.mlma', "gwas_files" = '/cluster/work/pausch/xena/gwas/new/{inheri}/{{trait}}/CHR{chr}/result.mlma.gz'),
    rule = 'format_covariates',
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
## covariatesFile <- "~/TWAS/SNAKES/train/testis_covariates.txt"
## geneExpFile <- "/cluster/work/pausch/naveen/TWAS/EXPR/testis/expr_trans.txt"

##covariatesFile <- snakemake@input [["covariates"]]
covariatesFile <- snakemake@input [[1]]
mphenoFile <- snakemake@input[["mpheno"]]
outfile <- snakemake@output[["outfile"]]



cat ('covariate file is:', covariatesFile, '\n' )
cat ('reading the covriate file\n')
covariates <- read.table (covariatesFile, header=T)
##remove the first column and use it as the row ids
rnames <- covariates [,1]
covariates <- covariates [, -c(1)]
rownames (covariates) <-  rnames

##is the colnames of covariates is in same order as the rownames of gene exp?
cat ('reading the molecular phenotype file \n')
mpheno <- read.table (mphenoFile, head=T, row.names=1,sep="\t", stringsAsFactor=F)

#same order as in the phenotype file
covariates <- covariates [, rownames(mpheno) ]
if ( sum (colnames (covariates) == rownames (mpheno)) == nrow (mpheno) ) {
    write.table (covariates,  file = outfile, sep = "\t",row.names=TRUE)
}else {
    stop ("order of ids in covariantes and exp files do not match")
    #sum (colnames (covariates) == rownames (mpheno)) == nrow (mpheno) 
}


#ordering columns is easy!
#x<-data.frame (matrix (1:100, nrow=10))
#colnames (x) <- LETTERS [1:10]
#x [, rev (LETTERS [1:10])]


## ##peer_factors = read.csv(file = "./output/peer_out/X.csv", header = FALSE)
## covariatess = read.csv(file = covariatesFile, header = TRUE)
## gene_exp_transpose <- read.table (file = geneExpFile, head=T, row.names=T)


## #Set the column names for the PEER factors (covariates) as the subject IDs
## colnames(peer_factors) = rownames(gene_exp_transpose)

## # write out a covariates matrix
## ##write.table(peer_factors, file = "./output/covariates.txt", sep = "\t",
## #row.names = TRUE)
## write.table(peer_factors, file = outfile, sep = "\t",
##             row.names = TRUE)
