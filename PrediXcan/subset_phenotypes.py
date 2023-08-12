import gzip
pheno = snakemake.input.pheno
vcf = snakemake.input.vcf
outfile = snakemake.output.pheno
out = open(outfile, "w")

vcf_ids = {}
print('reading vcf')
with gzip.open(vcf, "rt") as inf:
    for line in inf:
        print(line)
        if line[0:6] == "#CHROM":
            ids = line.rstrip().split("\t")[9:]
            print(f'Number of samples in the vcf ::  {len (ids)}')
            # vcf_ids = {myid: 1 for myid in ids} ##this does not work
            for myid in ids:
                vcf_ids[myid] = 1
            break

towrite = {}
print('readingt the phenotype file')
with open(pheno, "rt") as inf:
    for lnum, line in enumerate(inf):
        if lnum == 0:
            out.write(f'{line}')
        else:
            spl = line.rstrip().split()
            if spl[0] in vcf_ids:
                out.write(f'{line}')
