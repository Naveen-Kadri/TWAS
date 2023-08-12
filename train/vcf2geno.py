'''
Get genotypes and snp_annotation files from vcf
genotype file : VarID, sampleids...
snp annotation file :  chromosome, pos, varID, ref_vcf, alt_vcf, rsid
'''

import gzip
# vcf = '/cluster/work/pausch/naveen/ASE/VCF/chr_25_beagle5.2.vcf.gz'
# genotypefile = 'geno'
# annotationfile = 'annot'

vcf = snakemake.input.vcf
genotypefile = snakemake.output.genotypes
annotationfile = snakemake.output.annotation

geno = open(genotypefile, "w")
annot = open(annotationfile, "w")

recode = {
    "0|0": "0",
    "0|1": "1",
    "1|0": "1",
    "1|1": "2"
}

annot_header = "\t".join(
    ["chromosome", "pos", "varID", "ref_vcf", "alt_vcf", "rsid"])


with gzip.open(vcf, "rt") as inf:
    for line in inf:
        if line[0:2] != "##":
            spl = line.rstrip().split()
            if line[0:6] == "#CHROM":
                ids = spl[9:]
                header = "\t".join(["varID"] + ids)
                geno.write(f'{header}\n')
                annot.write(f'{annot_header}\n')
            else:
                gts = spl[9:]
                gts = [recode.get(gt) for gt in gts]
                tw = "\t".join([spl[2]] + gts)
                geno.write(f'{tw}\n')
                # varid for rsid
                tw = "\t".join(spl[:5] + [spl[2]])
                annot.write(f'{tw}\n')
