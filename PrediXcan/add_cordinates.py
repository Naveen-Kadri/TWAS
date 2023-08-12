import gzip
annotation = snakemake.input.annotation
result = snakemake.input.result
outfile = snakemake.output.outfile

out = gzip.open(outfile, 'wt')

# chr	gene_id	gene_name	start	end	gene_type
cords = {}
with open(annotation) as inf:
    for lnum, line in enumerate(inf):
        spl = line.rstrip().split()
        # print(spl)
        # if lnum > 3:
        #     exit()
        cords[spl[1]] = '\t'.join([spl[2], spl[0], spl[3], spl[4]])


#['gene', 'effect', 'se', 'zscore', 'pvalue', 'n_samples', 'status']
with open(result) as inf:
    for lnum, line in enumerate(inf):
        spl = line.rstrip().split()
        if lnum == 0:
            header = '\t'.join(spl + ['gene_name', 'chr', 'start', 'end'])
            out.write(f'{header}\n')
        else:
            tw = '\t'.join(spl + [cords[spl[0]]])
            out.write(f'{tw}\n')
