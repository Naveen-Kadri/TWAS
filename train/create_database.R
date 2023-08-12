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

