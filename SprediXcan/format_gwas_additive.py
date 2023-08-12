import gzip
# infiles = [
# y
# outfile = 'todel.txt.gz'

infiles = snakemake.input.infiles
outfile = snakemake.output.outfile


header = "\t" .join(['SNP', 'A1', 'A2', 'Freq', 'BETA', 'SE', 'P'])
out = open(outfile, "w")
out.write(f'{header}\n')

for infile in infiles:
    with gzip.open(infile, "rt") as inf:
        for lnum, line in enumerate(inf):
            if lnum == 0:
                continue
            if "nan" not in line:
                spl = line.rstrip().split()
                snp = "_".join([spl[1], spl[3], spl[4]])
                tw = "\t".join([snp] + spl[3:])
                out.write(f'{tw}\n')


out.close()
