infile = snakemake.input.mpheno
outfile = snakemake.output.mpheno


#infile = '/cluster/work/pausch/naveen/TWAS/splicing/PrediXcan/head/testis/CHR6/predicted_mpheno_id_order_fixed.txt'
#outfile = 'xx'

out = open(outfile, 'w')

sums = dict()
with open(infile, "rt") as inf:
    for lnum, line in enumerate(inf):
        spl = line.rstrip().split()
        if lnum == 0:
            header = spl
        else:
            info = dict(zip(header[2:], spl[2:]))
            for mygene in info:
                sums[mygene] = sums.get(mygene, 0) + float(info[mygene])

to_remove = {}
with open(infile, "rt") as inf:
    for lnum, line in enumerate(inf):
        spl = line.rstrip().split()
        if lnum == 0:
            for myindex, el in enumerate(spl):
                if el in sums and sums[el] == 0.0:
                    to_remove[myindex] = 1

        tw = []
        for xindex, el in enumerate(spl):
            if xindex not in to_remove:
                tw.append(el)
        tw = "\t".join(tw)
        out.write(f'{tw}\n')

out.close()
