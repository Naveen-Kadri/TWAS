from collections import defaultdict
import gzip
'''
get additive and non additive (combination of recessive and dominat)
'''

recessive_files = snakemake.input.recessive
dominant_files = snakemake.input.dominant
outfile = snakemake.output.outfile

#outfile = 'xx.txt'
# recessive_files = [
# f'/cluster/work/pausch/xena/gwas/new/recessive/vol/CHR{xchr}/result.mlma.gz' for xchr in range(25, 26)]
# dominant_files = [
# f'/cluster/work/pausch/xena/gwas/new/dominant/vol/CHR{xchr}/result.mlma.gz' for xchr in range(25, 26)]


# Chr     SNP     bp      A1      A2      Freq    b       se      p
header = "\t" .join(['SNP', 'A1', 'A2', 'Freq', 'BETA', 'SE', 'P'])
out = open(outfile, "w")
out.write(f'{header}\n')


pvalues = {}
infos = {}
for infile in recessive_files:
    with gzip.open(infile, "rt") as inf:
        for lnum, line in enumerate(inf):
            if lnum == 0:
                continue
            if "nan" not in line:
                spl = line.rstrip().split()
                snp = "_".join([spl[1], spl[3], spl[4]])
                pvalues[snp] = float(spl[-1])


print('reading first set of files')
for infile in dominant_files:
    with gzip.open(infile, "rt") as inf:
        for lnum, line in enumerate(inf):
            if lnum == 0:
                continue
            if "nan" not in line:
                spl = line.rstrip().split()
                tw = "\t".join([spl[1]] + spl[3:8])
                snp = "_".join([spl[1], spl[3], spl[4]])
                tw = "\t".join([snp] + spl[3:8])
                if snp in pvalues:
                    pvalue = min(pvalues[snp], float(spl[-1]))
                else:
                    pvalue = spl[-1]
                out.write(f"{tw}\t{pvalue}\n")

out.close()
