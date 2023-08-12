
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
    input = list('/cluster/work/pausch/naveen/TWAS/expression/SprediXcan/nonadditive/con/result_with_cords.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr1.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr2.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr3.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr4.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr5.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr6.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr7.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr8.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr9.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr10.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr11.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr12.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr13.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr14.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr15.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr16.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr17.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr18.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr19.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr20.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr21.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr22.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr23.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr24.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr25.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr26.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr27.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr28.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr29.txt.gz', "twas" = '/cluster/work/pausch/naveen/TWAS/expression/SprediXcan/nonadditive/con/result_with_cords.txt.gz', "gwas" = c('/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr1.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr2.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr3.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr4.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr5.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr6.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr7.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr8.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr9.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr10.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr11.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr12.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr13.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr14.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr15.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr16.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr17.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr18.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr19.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr20.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr21.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr22.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr23.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr24.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr25.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr26.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr27.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr28.txt.gz', '/cluster/work/pausch/xena/gwas/new/plots/con/forplot_chr29.txt.gz')),
    output = list('/cluster/work/pausch/naveen/TWAS/expression/SprediXcan/nonadditive/con/manhattan_expression_con_nonadditive.tiff', "plot_file" = '/cluster/work/pausch/naveen/TWAS/expression/SprediXcan/nonadditive/con/manhattan_expression_con_nonadditive.tiff'),
    params = list(6, 0.005, c('testis', 'epi_h', 'vas_d'), "thresh_line" = 6, "maf_thresh" = 0.005, "tissues" = c('testis', 'epi_h', 'vas_d')),
    wildcards = list('expression', 'nonadditive', 'con', "mpheno" = 'expression', "inheri" = 'nonadditive', "trait" = 'con'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', 'walltime', "mem_mb" = 8000, "mem_mib" = 7630, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/scratch/257320170.tmpdir', "walltime" = '01:00'),
    config = list("molecular_pheno" = list("expression" = '/cluster/work/pausch/naveen/TWAS/TPM/{tissue}.tsv', "splicing" = '/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/{tissue}.bed.gz'), "covariates" = list("expression" = '/cluster/work/pausch/xena/eQTL/cis_all/{tissue}/{tissue}_covariates.txt', "splicing" = '/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/covariates_{tissue}.txt'), "vcfs" = list("twas" = '/cluster/work/pausch/naveen/ASE/VCF/chr_{chr}_beagle5.2.vcf.gz', "gwas" = '/cluster/work/pausch/naveen/SNPDATA/TOSEQ/CHR{chr}/combined.vcf.gz'), "gtf" = '/cluster/work/pausch/naveen/ASE/INPUTS/Bos_taurus.ARS-UCD1.2.104.gtf.gz', "OUT_DIR" = '/cluster/work/pausch/naveen/TWAS', "old_gwas_files" = '/cluster/work/pausch/xena/gwas/gcta/results/{{trait}}/BV/sire/CHR{chr}/result.mlma', "gwas_files" = '/cluster/work/pausch/xena/gwas/new/{inheri}/{{trait}}/CHR{chr}/result.mlma.gz'),
    rule = 'manhattan2',
    bench_iteration = as.numeric(NA),
    scriptdir = '/cluster/home/nkadri/TWAS/SNAKES/spredi',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
cnames  <- c('pos', 'maf', 'add', 'nonadd')
colors  <- c('gray', 'darkgray')

fontsize=1.5
tissues <- c('testis', 'epi_h', 'vas_d')
tissue_colors <- c("#B8DAF3", "#CCE085", "#D4C8E3")


## gwas_files = paste0('/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr',1:29,'.txt.gz')
## gwas_files <- gwas_files [1:6]
## twas_file <- '/cluster/work/pausch/naveen/TWAS/expression/SprediXcan/nonadditive/tail/result_with_cords.txt.gz'
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
cat ('inheri is' , inheri, '\n')

if (inheri == 'additive') {
    pcol <- 4
}else{
    pcol <-5
}



cat ('plotting to the file ', plot_file, "\n")
tiff (plot_file, height=12,width=24,units="in",res=300)
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
#bck <- res


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


ylims <- c ( -1* -log10(min( res[, tissues],na.rm=TRUE)),   -log10(min (res [, pcol],na.rm=TRUE)))

cat ('ylims = ', ylims , '\n')

par (mar = c (4,6,4,4), cex.axis=1.5, cex.lab=1.5 )


##quick plotting
##res <- res [c (which (is.na (res[,pcol])),which(res[,pcol] < 10^-6)) ,]


##dim (xres)
plot (res$newpos, -log10(res[, pcol]), col=res$mycol, xaxt="n", ylab=expression (-log [10](italic("P"))), ylim=ylims, cex.main=1 , cex.lab=fontsize, cex.axis=1, xlab="", pch=20, cex=2,yaxt="n", main='')

for (ti in 1:length (tissues) ) {
    mytissue <- res [, c("chr","newpos", "gene", "gene_name", tissues [ti], "mycol") ]
    #mytissue$col <- mycol
    mytissue [which(mytissue [,5] < 1e-6), 'mycol' ] <- tissue_colors [ti]
    points (mytissue$newpos, -1*-log10(mytissue[,5]), col=mytissue$mycol, xaxt="n", ylab=expression (-log [10](italic("P"))), ylim=ylims, cex.main=1 , cex.lab=fontsize, cex.axis=1, xlab="", pch=20, cex=2,yaxt="n", main='')
    ft <- mytissue [!is.na (mytissue [,5]) ,]
    peaks <- data.frame ()

    for (chr in 1:29) {
        mychr <- ft [ft$chr == chr,]
        mychr <- mychr [-log10(mychr[,5])>thresh_line ,]
        if (nrow(mychr) >0) {
            mycut<-quantile(-log10(mychr[,5]), 0.95)
            mypeak<-mychr [-log10(mychr[,5]) > mycut, ]
            peaks <- rbind (peaks, mypeak)
        }
    }

    
    ##when gene names are missing
    peaks [is.na (peaks$gene_name), "gene_name"] <- peaks [is.na (peaks$gene_name), "gene"]
    if (nrow (peaks) >0 ) {
        text (peaks$newpos, -1*-log10(peaks[,5]), pos=c(2,4), peaks$gene_name, col=peaks$col, font=2)
    }
}



abline (h =0, lty=2)

yats <- seq(floor (ylims [1] /10) *10, ceiling(ylims[2]/10) *10, length.out=11)
axis (side=2, at=yats, labels=abs(yats),las=1)
abline ( h = c(-1,1) *thresh_line, lty=2, col='orange')

##vertical lines separting the chromosomes
abline (v= sizes + toadd,col='gray', lty=2 )


ats  <- c ()
for (chr in chrs) {
    myats <- mean (res [ res [,chrcol] == chr ,"newpos" ])
    ats  <- c (ats, myats)
    cat (chr, "\n")
}

text (ats, rep(ylims[2], length(ats)), labels=chrs, cex=1.5)
dev.off ()







