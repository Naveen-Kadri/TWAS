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
