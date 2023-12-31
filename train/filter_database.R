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
