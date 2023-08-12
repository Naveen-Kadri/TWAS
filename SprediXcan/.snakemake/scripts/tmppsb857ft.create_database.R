
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
    input = list('/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR1/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR2/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR3/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR4/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR5/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR6/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR7/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR8/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR9/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR10/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR11/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR12/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR13/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR14/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR15/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR16/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR17/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR18/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR19/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR20/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR21/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR22/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR23/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR25/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR26/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR27/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR28/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR29/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR1/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR2/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR3/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR4/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR5/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR6/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR7/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR8/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR9/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR10/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR11/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR12/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR13/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR14/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR15/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR16/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR17/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR18/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR19/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR20/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR21/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR22/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR23/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR25/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR26/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR27/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR28/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR29/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR1/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR2/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR3/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR4/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR5/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR6/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR7/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR8/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR9/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR10/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR11/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR12/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR13/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR14/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR15/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR16/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR17/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR18/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR19/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR20/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR21/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR22/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR23/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR25/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR26/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR27/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR28/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR29/tiss_summ.txt', "model_summary_files" = c('/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR1/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR2/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR3/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR4/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR5/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR6/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR7/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR8/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR9/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR10/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR11/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR12/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR13/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR14/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR15/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR16/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR17/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR18/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR19/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR20/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR21/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR22/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR23/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR25/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR26/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR27/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR28/summary.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR29/summary.txt'), "weights_files" = c('/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR1/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR2/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR3/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR4/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR5/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR6/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR7/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR8/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR9/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR10/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR11/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR12/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR13/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR14/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR15/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR16/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR17/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR18/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR19/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR20/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR21/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR22/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR23/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR25/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR26/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR27/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR28/weights.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR29/weights.txt'), "tiss_chr_summ_files" = c('/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR1/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR2/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR3/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR4/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR5/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR6/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR7/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR8/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR9/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR10/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR11/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR12/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR13/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR14/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR15/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR16/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR17/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR18/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR19/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR20/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR21/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR22/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR23/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR24/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR25/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR26/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR27/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR28/tiss_summ.txt', '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR29/tiss_summ.txt')),
    output = list('/cluster/work/pausch/naveen/TWAS/splicing/vas_d/prediction_models.db', "database" = '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/prediction_models.db'),
    params = list(29, "nchr" = 29),
    wildcards = list('splicing', 'vas_d', "mpheno" = 'splicing', "tissue" = 'vas_d'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', 'walltime', "mem_mb" = 16000, "mem_mib" = 15259, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/scratch/257474238.tmpdir', "walltime" = '04:00'),
    config = list("molecular_pheno" = list("expression" = '/cluster/work/pausch/naveen/TWAS/TPM/{tissue}.tsv', "splicing" = '/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/{tissue}.bed.gz'), "covariates" = list("expression" = '/cluster/work/pausch/xena/eQTL/cis_all/{tissue}/{tissue}_covariates.txt', "splicing" = '/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/covariates_{tissue}.txt'), "vcfs" = list("twas" = '/cluster/work/pausch/naveen/ASE/VCF/chr_{chr}_beagle5.2.vcf.gz', "gwas" = '/cluster/work/pausch/naveen/SNPDATA/TOSEQ/CHR{chr}/combined.vcf.gz'), "gtf" = '/cluster/work/pausch/naveen/ASE/INPUTS/Bos_taurus.ARS-UCD1.2.104.gtf.gz', "OUT_DIR" = '/cluster/work/pausch/naveen/TWAS', "old_gwas_files" = '/cluster/work/pausch/xena/gwas/gcta/results/{{trait}}/BV/sire/CHR{chr}/result.mlma', "gwas_files" = '/cluster/work/pausch/xena/gwas/new/{inheri}/{{trait}}/CHR{chr}/result.mlma.gz'),
    rule = 'create_database',
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
##use r/4.2.2

##dependencies=TRUE
##also install all dependencies with the myconda.yaml
## install.packages ('DBI',lib='/cluster/work/pausch/naveen/Rlib/')
## install.packages ('DBI', lib='/cluster/work/pausch/naveen/Rlib/')
## install.packages ('RSQLite', lib='/cluster/work/pausch/naveen/Rlib/')
## install.packages ('xml2', lib='/cluster/work/pausch/naveen/Rlib/',dependencies=TRUE)
## install.packages ('tidyverse', lib='/cluster/work/pausch/naveen/Rlib/',dependencies=TRUE)

.libPaths("/cluster/work/pausch/naveen/Rlib/")
library (DBI)
library (RSQLite)
library ('tidyverse')


model_summary_files <-  snakemake@input[["model_summary_files"]]
tiss_summary_files <- snakemake@input [["tiss_chr_summ_files"]]
weight_files <- snakemake@input [["weights_files"]]
database <- snakemake@output[["database"]]
cat ('the data base is',database, '\n')
nchr <- snakemake@params[["nchr"]]
"%&%" <- function(a,b) paste(a,b, sep='')

## model_summary_files <- paste0("/cluster/work/pausch/naveen/TWAS/TRAIN/testis/CHR", 1:3, "/summary.txt")
## tiss_summary_files <- paste0("/cluster/work/pausch/naveen/TWAS/TRAIN/testis/CHR", 1:3, "/tiss_summ.txt")
## weight_files <- paste0("/cluster/work/pausch/naveen/TWAS/TRAIN/testis/CHR", 1:3, "/weights.txt")
## nchr<-3
## database<-"test.db"

driver <- dbDriver('SQLite')
model_summaries <- read.table(model_summary_files[1],header = T, stringsAsFactors = F)
tiss_summaries <- read.table(tiss_summary_files [1], header = T, stringsAsFactors = F)
n_samples <- tiss_summaries$n_samples
  
for (i in 2:nchr) {
    cat ('reading summary for chr',i,'\n')
    model_summaries <- rbind(model_summaries,read.table(model_summary_files [i], header = T, stringsAsFactors = F))
    cat ('reading tiss summary for chr',i, '\n')
    tiss_summaries <- rbind(tiss_summaries,read.table(tiss_summary_files[i], header = T, stringsAsFactors = F))
}

colnames (model_summaries)
#model_summaries <- rename(model_summaries, gene = gene_id)
colnames (model_summaries)  [which (colnames (model_summaries) == 'gene_id')] = "gene"


# Create a database connection
##conn <- dbConnect(drv = driver, './dbs/gtex_v7_models.db')
conn <- dbConnect(drv = driver, database)

dbWriteTable(conn, 'model_summaries', model_summaries, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX gene_model_summary ON model_summaries (gene)")

# Weights Table -----
##weights <- read.table('./weights/Model_training_chr1_weights.txt', header = T,stringsAsFactors = F)
weights <- read.table(weight_files [1], header = T,stringsAsFactors = F)


for (i in 2:nchr) {
  weights <- rbind(weights,
              read.table(weight_files [i], header = T, stringsAsFactors = F))
  
}
  
##weights <- rename(weights, gene = gene_id)
colnames (weights)  [which (colnames (weights) == 'gene_id')] = "gene"
dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbExecute(conn, "CREATE INDEX weights_gene ON weights (gene)")
dbExecute(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")

# Sample_info Table ----
sample_info <- data.frame(n_samples = n_samples, population = 'BSW') # Provide the population info
dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)
  
# Construction Table ----
construction <- tiss_summaries %>%
                    select(chrom, cv_seed) %>%
                    rename(chromosome = chrom)

dbWriteTable(conn, 'construction', construction, overwrite = TRUE)
dbDisconnect(conn)




##check the database
##mydb <- dbConnect(RSQLite::SQLite(), "test.db")
##mywts <- dbReadTable(mydb, 'weights')
##dim (mywts)
##table (sapply (strsplit(mywts$rsid, "_"), "[", 1))

