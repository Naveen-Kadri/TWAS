cnames  <- c('pos', 'maf', 'add', 'nonadd')
colors  <- c('gray', 'black')
fontsize=1.5


plot_file <- snakemake@output[["plot_file"]]
gwas_files <- snakemake@input[["gwas"]]
twas_files <- snakemake@input[["twas"]]
maf_thresh <- as.numeric(snakemake@params[["maf_thresh"]])
thresh_line  <-snakemake@params[["thresh_line"]]
inheri <- snakemake@wildcards[["inheri"]]

cat ('inheri is' , inheri, '\n')



if (inheri == 'additive') {
    pcol <- 4
}else{
    pcol <-5
}



cat ('plotting to the file ', plot_file, "\n")
tiff (plot_file, height=12,width=20,units="in",res=300)
maxis  <- c()
res <- data.frame ()
nvar <- 0


#gwas_files <- paste0('/cluster/work/pausch/naveen/RECOMBINATION/REFALT/IMPUTE/GWAS/GRR/bv/female/forplot_chr',1:29,'.txt.gz')
##twas_files = paste0('/cluster/work/pausch/naveen/TWAS/expression/PrediXcan/head/testis/CHR',1:29,'/result_with_cords.txt.gz')
sizes <-rep (0, length (gwas_files) )
for ( k in 1:length (gwas_files) ) {
    myres  <- matrix (scan(gwas_files [k]), ncol=4, byrow=T)
    ##add chr col
    #myres <- data.frame (cbind (k, myres))
    myres <- data.frame (myres [order (myres[,2]),])
    
    colnames (myres) <- cnames
    nvar <- nvar + nrow (myres)
    #filter only on the gwas results
    myres <- myres [ myres$maf > maf_thresh ,]


    inf<- read.table (twas_files [k], header=TRUE,sep="\t")
    inf$mid <- (inf$start + inf$end ) /2
    inf<-  inf[, c( 'mid', 'pvalue')]
    
    all <- merge (myres, inf, by.x='pos', by.y='mid', all=TRUE)
    all <- data.frame (chr=k, all)
    colnames (all) [6] <- 'twas'
    
    res <- rbind (res,all)
    sizes [k] <- max (all$pos)
    cat (k, "\n")
}



cat ("number of variants tested, passing maf thresh", nvar/1e6, nrow(res)/1e6, "mi.\n" )
chrcol <-1
poscol <- 2

mycol = rep (NA, nrow (res))
chrs  <- unique (res [, chrcol]) ; chrs
evens <- chrs [1:length(chrs)%%2==0]
odds <- chrs [1:length(chrs)%%2==1]
mycol [which (res[,chrcol] %in% evens) ]  <-  colors [1]
mycol [which (res[,chrcol] %in% odds) ]  <-  colors [2]


cat (sizes)
##continuous positions
toadd <- c (0, cumsum (sizes))
toadd <- toadd [-c(length(toadd))]
newpos <- res[, poscol] + toadd [res [, chrcol]]
maxi <- -log10(min (res[,pcol], res$twas,na.rm=TRUE))
ylims <- c(-1,1) *maxi ; ylims
cat ('ylims = ', ylims , '\n')

par (mar = c (4,6,4,4), cex.axis=1.5, cex.lab=1.5 )
plot (newpos, -log10(res[, pcol]), col=mycol, xaxt="n", ylab=expression (-log [10](italic("P"))), ylim=ylims, cex.main=1 , cex.lab=fontsize, cex.axis=1, xlab="", pch=20, cex=1.5,yaxt="n", main='')
points (newpos, -1*-log10(res$twas), col=mycol, xaxt="n", ylab=expression (-log [10](italic("P"))), ylim=ylims, cex.main=1 , cex.lab=fontsize, cex.axis=1, xlab="", pch=20, cex=1.5,yaxt="n", main='')
abline (h =0, lty=2)

yats <- seq(floor (ylims [1] /10) *10, ceiling(ylims[2]/10) *10, length.out=11)
axis (side=2, at=yats, labels=abs(yats),las=1)
abline ( h = c(-1,1) *thresh_line, lty=2, col='red')

##vertical lines separting the chromosomes
abline (v= sizes + toadd,col='gray', lty=2 )


ats  <- c ()
for (chr in chrs) {
    myats <- mean (newpos [res [, chrcol] == chr  ])
    ##myats  <-  mean (range (res [res[,chrcol] == chr,"pos"]))
    ats  <- c (ats, myats)
    cat (chr, "\n")
}


text (ats, rep(maxi, length(ats)), labels=chrs, cex=1.5)
dev.off ()

#ylims <- c(-1,1) * 48 ; ylims
#yats <- seq(floor (ylims [1] /10) *10, ceiling(ylims[2]/10) *10, length.out=11) ;yats

