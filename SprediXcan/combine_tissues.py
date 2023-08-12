import gzip
from collections import defaultdict
# tissues = ['testis', 'epi_h', 'vas_d']
# infiles = [f'/cluster/work/pausch/naveen/TWAS/expression/SprediXcan/additive/mot/{tissue}/result_with_cords.txt.gz' for tissue in tissues ]
# outfile = "todel.txt"

tissues = snakemake.params.tissues
infiles = snakemake.input.infiles
outfile = snakemake.output.outfile

out = gzip.open (outfile, "wt")

pvalues = defaultdict (lambda : ["NA"]*3)
header = "\t".join(['gene', 'gene_name', 'chr', 'start', 'end'] + tissues )
out.write (f'{header}\n')
for t,infile in enumerate(infiles):
    with gzip.open (infile, "rt") as inf:
        for lnum,line in enumerate(inf):
            if lnum ==0:
                continue
            spl = line.rstrip().split()
            mykey="\t".join ( [ spl [0], spl[1] ]  + spl[12:15] )
            pvalues [mykey] [t] = spl [4]



for cord in pvalues:
    mycord = pvalues [cord]
    out.write (f'{cord}\t{mycord[0]}\t{mycord[1]}\t{mycord[2]}\n')


out.close()
