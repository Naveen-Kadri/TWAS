
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
    input = list('/cluster/work/pausch/naveen/TWAS/splicing/vas_d/prediction_models.db', "database" = '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/prediction_models.db'),
    output = list('/cluster/work/pausch/naveen/TWAS/splicing/vas_d/filtered_signif.db', "database" = '/cluster/work/pausch/naveen/TWAS/splicing/vas_d/filtered_signif.db'),
    params = list(),
    wildcards = list('splicing', 'vas_d', "mpheno" = 'splicing', "tissue" = 'vas_d'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', 'walltime', "mem_mb" = 16000, "mem_mib" = 15259, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/scratch/257474244.tmpdir', "walltime" = '04:00'),
    config = list("molecular_pheno" = list("expression" = '/cluster/work/pausch/naveen/TWAS/TPM/{tissue}.tsv', "splicing" = '/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/{tissue}.bed.gz'), "covariates" = list("expression" = '/cluster/work/pausch/xena/eQTL/cis_all/{tissue}/{tissue}_covariates.txt', "splicing" = '/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/covariates_{tissue}.txt'), "vcfs" = list("twas" = '/cluster/work/pausch/naveen/ASE/VCF/chr_{chr}_beagle5.2.vcf.gz', "gwas" = '/cluster/work/pausch/naveen/SNPDATA/TOSEQ/CHR{chr}/combined.vcf.gz'), "gtf" = '/cluster/work/pausch/naveen/ASE/INPUTS/Bos_taurus.ARS-UCD1.2.104.gtf.gz', "OUT_DIR" = '/cluster/work/pausch/naveen/TWAS', "old_gwas_files" = '/cluster/work/pausch/xena/gwas/gcta/results/{{trait}}/BV/sire/CHR{chr}/result.mlma', "gwas_files" = '/cluster/work/pausch/xena/gwas/new/{inheri}/{{trait}}/CHR{chr}/result.mlma.gz'),
    rule = 'filter_database',
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
#unfiltered_db <- './dbs/gtex_v7_models.db'
#filtered_db <- './dbs/gtex_v7_models_filtered_signif.db'


unfiltered_db <- snakemake@input[["database"]]
filtered_db <- snakemake@output[["database"]]
.libPaths("/cluster/work/pausch/naveen/Rlib/")
library ('DBI')
library ('RSQLite')
library ('tidyverse')

driver <- dbDriver("SQLite")
in_conn <- dbConnect(driver, unfiltered_db)
out_conn <- dbConnect(driver, filtered_db)
model_summaries <- dbGetQuery(in_conn, 'select * from model_summaries where zscore_pval < 0.05 and rho_avg > 0.1')
model_summaries <- model_summaries %>% 
                    rename(pred.perf.R2 = rho_avg_squared, genename = gene_name, pred.perf.pval = zscore_pval, n.snps.in.model = n_snps_in_model)
model_summaries$pred.perf.qval <- NA
dbWriteTable(out_conn, 'extra', model_summaries, overwrite = TRUE)
construction <- dbGetQuery(in_conn, 'select * from construction')
dbWriteTable(out_conn, 'construction', construction, overwrite = TRUE)
sample_info <- dbGetQuery(in_conn, 'select * from sample_info')
dbWriteTable(out_conn, 'sample_info', sample_info, overwrite = TRUE)
weights <- dbGetQuery(in_conn, 'select * from weights')
weights <- weights %>% filter(gene %in% model_summaries$gene) %>% rename(eff_allele = alt, ref_allele = ref, weight = beta)
dbWriteTable(out_conn, 'weights', weights, overwrite = TRUE)
dbExecute(out_conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbExecute(out_conn, "CREATE INDEX weights_gene ON weights (gene)")
dbExecute(out_conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
dbExecute(out_conn, "CREATE INDEX gene_model_summary ON extra (gene)")
dbDisconnect(in_conn)
dbDisconnect(out_conn)
