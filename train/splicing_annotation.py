import gzip
from collections import defaultdict

gene_annotation = snakemake.input.gene_annotation
splicing = snakemake.input.splicing
outfile = snakemake.output.outfile
#gene_annotation = '/cluster/work/pausch/naveen/TWAS/gene_annotation.txt'
#splicing = '/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/testis.bed.gz'
#outfile = 'test.txt'


out = open(outfile, 'w')


stable = defaultdict(list)
with gzip.open(splicing, 'rt') as inf:
    for line in inf:
        spl = line.rstrip().split()
        stable[spl[4]].append(spl[3])


# chr	gene_id	gene_name	start	end	gene_type
with open(gene_annotation) as inf:
    for lnum, line in enumerate(inf):
        if lnum == 0:
            out.write(f'{line}')
        else:
            spl = line.rstrip().split()
            # combination of gene_id and gene_name in the third column - this shuld be unique so add lnum
            spl[2] = "_".join([spl[1], spl[2], str(lnum)])
            myclusters = stable[spl[1]]
            for mycluster in myclusters:
                # R rownames cannot start with numbers
                spl[1] = "G" + mycluster.replace(":", "_")
                tw = "\t".join(spl)
                out.write(f'{tw}\n')
