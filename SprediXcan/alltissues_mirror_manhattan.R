cnames  <- c('pos', 'maf', 'add', 'nonadd')
colors  <- c('gray', 'darkgray')

fontsize=1.5
tissues <- c('testis', 'epi_h', 'vas_d')
tissue_colors <- c("#B8DAF3", "#CCE085", "#D4C8E3")
tissue_names <- c("Testis", "Epididymis", "Vas deferens")

## gwas_files = paste0('/cluster/work/pausch/xena/gwas/new/plots/fbkall/forplot_chr',1:29,'.txt.gz')
## twas_file <- '/cluster/work/pausch/naveen/TWAS/splicing/PrediXcan/nonadditive/fbkall/result_with_cords.txt.gz'
## #twas_file <- '/cluster/work/pausch/naveen/TWAS/splicing/PrediXcan/fbkall/result_with_cords.txt.gz'
## inheri = 'nonadditive'
## maf_thresh <- 0.005
## thresh_line <- 6


plot_file <- snakemake@output[["plot_file"]]
gwas_files <- snakemake@input[["gwas"]]
##gwas_files <- gwas_files [1:6]
twas_file <- snakemake@input[["twas"]]
maf_thresh <- as.numeric(snakemake@params[["maf_thresh"]])
thresh_line  <-snakemake@params[["thresh_line"]]
inheri <- snakemake@wildcards[["inheri"]]
tissues <- snakemake@params [["tissues"]]
mpheno <- snakemake@wildcards [["mpheno"]]
cat ('inheri is' , inheri, '\n')

if (inheri == 'additive') {
    pcol <- 4
}else{
    pcol <-5
}



cat ('plotting to the file ', plot_file, "\n")
tiff(plot_file, height=12,width=26,units="in",res=300)
maxis  <- c()
res <- data.frame ()
nvar <- 0


twas<- read.table (twas_file, header=TRUE,sep="\t")
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



    ##TWAS file
    inf <- twas [twas$chr == k ,]
    inf$mid <- (inf$start + inf$end ) /2
    inf<-  inf[, c( 'mid', "gene", "gene_name", tissues)]
    
    all <- merge (myres, inf, by.x='pos', by.y='mid', all=TRUE)
    all <- data.frame (chr=k, all)
    all <- all [order (all$pos) ,]
    #colnames (all) [6] <- 'twas'
    
    res <- rbind (res,all)
    sizes [k] <- max (all$pos)
    cat (k, "\n")
}
##bck <- res


cat ('the minimum pvalues are', min(res[,pcol],na.rm=TRUE), min(res[, tissues],na.rm=TRUE), "\n" )

cat ("number of variants tested, passing maf thresh", nvar/1e6, nrow(res)/1e6, "mi.\n" )
chrcol <-1
poscol <- 2

mycol = rep (NA, nrow (res))
chrs  <- unique (res [, chrcol]) ; chrs
evens <- chrs [1:length(chrs)%%2==0]
odds <- chrs [1:length(chrs)%%2==1]
mycol [which (res[,chrcol] %in% evens) ]  <-  colors [1]
mycol [which (res[,chrcol] %in% odds) ]  <-  colors [2]
res$mycol <- mycol

cat (sizes)
##continuous positions
toadd <- c (0, cumsum (sizes))
toadd <- toadd [-c(length(toadd))]
res$newpos <- res[, poscol] + toadd [res [, chrcol]]




ylims <-1.1 * c ( -1* -log10(min( res[, tissues],na.rm=TRUE)),   -log10(min (res [, pcol],na.rm=TRUE)))

cat ('ylims = ', ylims , '\n')

par (mar = c (4,6,4,4), cex.axis=1.5, cex.lab=1.5 )


##quick plotting
##res <- res [c (which (is.na (res[,pcol])),which(res[,pcol] < 10^-6)) ,]


##dim (xres)

##HERE
xlims <- c(-25000000, max (res$newpos)*1.01 )

plot (res$newpos, -log10(res[, pcol]), col=res$mycol, xaxt="n", ylab=expression (-log [10](italic("P"))), ylim=ylims, cex.main=1 , cex.lab=fontsize, cex.axis=1, xlab="", pch=20, cex=2,yaxt="n", main='',xlim=xlims)

yats <- seq(floor (ylims [1] /10) *10, ceiling(ylims[2]/10) *10, length.out=11)
axis (side=2, at=yats, labels=abs(yats),las=1)


##vertical lines separting the chromosomes
abline (v= sizes + toadd,col='gray', lty=2 )

ats  <- c ()
for (chr in chrs) {
    myats <- mean (res [ res [,chrcol] == chr ,"newpos" ])
    ats  <- c (ats, myats)
    cat (chr, "\n")
}
text (ats, rep(ylims[2], length(ats)), labels=chrs, cex=1.5)
abline ( h = c(-1,1) *thresh_line, lty=5, col='orange')

##PEAKS <- list ()
PEAKS <- data.frame ()
for (ti in 1:length (tissues) ) {
##for (ti in 1) {
    mytissue <- res [, c("chr","newpos", "gene", "gene_name", tissues [ti], "mycol") ]
    #mytissue$col <- mycol
    mytissue [which(mytissue [,5] < 1e-6), 'mycol' ] <- tissue_colors [ti]
    points (mytissue$newpos, -1*-log10(mytissue[,5]), col=mytissue$mycol, xaxt="n", ylab=expression (-log [10](italic("P"))), ylim=ylims, cex.main=1 , cex.lab=fontsize, cex.axis=1, xlab="", pch=20, cex=2,yaxt="n", main='')
    ft <- mytissue [!is.na (mytissue [,5]) ,]
    peaks <- data.frame ()

    for (chr in 1:29) {
        mychr <- ft [ft$chr == chr,]
        mychr <- mychr [-log10(mychr[,5])>thresh_line ,]
        ##mychr <- mychr [!duplicated (mychr$gene_id) ,]
        cat (chr, nrow (mychr) )


        if (nrow(mychr) >0) {
            #mycut<-quantile(-log10(mychr[,5]), 0.95)
            #mypeak<-mychr [-log10(mychr[,5]) > mycut, ]
            mychr <- mychr [order (mychr[,5]) ,]
            mypeak <- head (mychr,2)
            peaks <- rbind (peaks, mypeak)
        }
    }

    if (nrow (peaks) > 0 ) {
        colnames (peaks) [5] <- "pvalue"
        peaks$tissue <- tissues [ti]
        
        PEAKS <- rbind (PEAKS, peaks)
    }
}


##EXTRACT gene name 
if (nrow (PEAKS) > 0  ) { 
    if (mpheno == "splicing") {
        gene_ids <- sapply (strsplit (PEAKS$gene_name, "_"), "[", 1)
        gene_names<-sapply (strsplit (PEAKS$gene_name, "_"), "[", 2)
        gene_names [gene_names == "NA"] <- gene_ids [gene_names =="NA"]
        PEAKS$gene_name <- gene_names
    }

    PEAKS <- PEAKS [order (PEAKS$newpos) ,]
    noname <- which(is.na(PEAKS$gene_name) )
    ##set missing gene_names to gene id if mpheno == expr
    PEAKS [noname, "gene_name"] <- PEAKS [noname, "gene"] 

    text (PEAKS$newpos, -1.02*-log10(PEAKS$pvalue), labels = PEAKS$gene_name,col=PEAKS$mycol, pos= c(4,1,2), font=2)
    ##text (PEAKS$newpos, -1*-log10(PEAKS$pvalue), labels = PEAKS$gene_name,col=PEAKS$mycol, pos= c(4,1,2), font=2)
}


##attempt to add arrows and text
## arrows (PEAKS$newpos, -1*-log10(PEAKS$pvalue),  PEAKS$newpos*1.1, -1*-log10(PEAKS$pvalue), angle=15, col = PEAKS$mycol,length=0.1 )
## text (1.1*PEAKS$newpos, -1*-log10(PEAKS$pvalue), labels = PEAKS$gene_name,col=PEAKS$mycol,  font=2,pos=4)


legend ('bottomright', col=tissue_colors, legend=tissue_names, pch=20, pt.cex=2, bty="n",cex=1.5)

abline (h =0, lty=2)
dev.off ()


##PEAKS [PEAKS$tissue == "epi_h" ,]






