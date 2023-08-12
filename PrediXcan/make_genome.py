import gzip

infiles = snakemake.input.infiles
outfile = snakemake.output.outfile
out = gzip.open (outfile, "wt")
for fnum,infile in enumerate (infiles):
    with gzip.open (infile, "rt") as inf:
        for i, line in enumerate (inf):
            if i==0:
                if fnum==0:
                    out.write (line)
            else:
                out.write (line)    
