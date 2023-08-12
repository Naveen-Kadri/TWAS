
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/cluster/work/pausch/naveen/mambaforge/envs/snakemake/lib/python3.11/site-packages', '/cluster/home/nkadri/.cache/snakemake/snakemake/source-cache/runtime-cache/tmpdq78hpxe/file/cluster/home/nkadri/TWAS/SNAKES/spredi', '/cluster/home/nkadri/TWAS/SNAKES/spredi']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\xee\x19\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8cG/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR1/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR2/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR3/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR4/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR5/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR6/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR7/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR8/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR9/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR10/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR11/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR12/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR13/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR14/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR15/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR16/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR17/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR18/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR19/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR20/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR21/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR22/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR23/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR24/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR25/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR26/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR27/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR28/result.mlma.gz\x94\x8cH/cluster/work/pausch/xena/gwas/new/recessive/fbkall/CHR29/result.mlma.gz\x94\x8cF/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR1/result.mlma.gz\x94\x8cF/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR2/result.mlma.gz\x94\x8cF/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR3/result.mlma.gz\x94\x8cF/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR4/result.mlma.gz\x94\x8cF/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR5/result.mlma.gz\x94\x8cF/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR6/result.mlma.gz\x94\x8cF/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR7/result.mlma.gz\x94\x8cF/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR8/result.mlma.gz\x94\x8cF/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR9/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR10/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR11/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR12/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR13/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR14/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR15/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR16/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR17/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR18/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR19/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR20/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR21/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR22/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR23/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR24/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR25/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR26/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR27/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR28/result.mlma.gz\x94\x8cG/cluster/work/pausch/xena/gwas/new/dominant/fbkall/CHR29/result.mlma.gz\x94e}\x94(\x8c\x06_names\x94}\x94(\x8c\trecessive\x94K\x00K\x1d\x86\x94\x8c\x08dominant\x94K\x1dK:\x86\x94u\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94ehM\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(hS)}\x94\x8c\x05_name\x94hMsNt\x94bhNhQhS\x85\x94R\x94(hS)}\x94hWhNsNt\x94bhGh\x06\x8c\tNamedlist\x94\x93\x94)\x81\x94(h\nh\x0bh\x0ch\rh\x0eh\x0fh\x10h\x11h\x12h\x13h\x14h\x15h\x16h\x17h\x18h\x19h\x1ah\x1bh\x1ch\x1dh\x1eh\x1fh h!h"h#h$h%h&e}\x94(hE}\x94hK]\x94(hMhNehMhQhS\x85\x94R\x94(hS)}\x94hWhMsNt\x94bhNhQhS\x85\x94R\x94(hS)}\x94hWhNsNt\x94bubhIh^)\x81\x94(h\'h(h)h*h+h,h-h.h/h0h1h2h3h4h5h6h7h8h9h:h;h<h=h>h?h@hAhBhCe}\x94(hE}\x94hK]\x94(hMhNehMhQhS\x85\x94R\x94(hS)}\x94hWhMsNt\x94bhNhQhS\x85\x94R\x94(hS)}\x94hWhNsNt\x94bubub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8cA/cluster/work/pausch/naveen/TWAS/GWAS/nonadditive/fbkall/gwas.txt\x94a}\x94(hE}\x94\x8c\x07outfile\x94K\x00N\x86\x94shK]\x94(hMhNehMhQhS\x85\x94R\x94(hS)}\x94hWhMsNt\x94bhNhQhS\x85\x94R\x94(hS)}\x94hWhNsNt\x94bh~h{ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(hE}\x94hK]\x94(hMhNehMhQhS\x85\x94R\x94(hS)}\x94hWhMsNt\x94bhNhQhS\x85\x94R\x94(hS)}\x94hWhNsNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94\x8c\x06fbkall\x94a}\x94(hE}\x94\x8c\x05trait\x94K\x00N\x86\x94shK]\x94(hMhNehMhQhS\x85\x94R\x94(hS)}\x94hWhMsNt\x94bhNhQhS\x85\x94R\x94(hS)}\x94hWhNsNt\x94bh\x9fh\x9cub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01M@\x1fM\xce\x1dM\xe8\x03M\xba\x03\x8c\x19/scratch/257492448.tmpdir\x94\x8c\x0504:00\x94e}\x94(hE}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06mem_mb\x94K\x02N\x86\x94\x8c\x07mem_mib\x94K\x03N\x86\x94\x8c\x07disk_mb\x94K\x04N\x86\x94\x8c\x08disk_mib\x94K\x05N\x86\x94\x8c\x06tmpdir\x94K\x06N\x86\x94\x8c\x08walltime\x94K\x07N\x86\x94uhK]\x94(hMhNehMhQhS\x85\x94R\x94(hS)}\x94hWhMsNt\x94bhNhQhS\x85\x94R\x94(hS)}\x94hWhNsNt\x94bh\xb3K\x01h\xb5K\x01h\xb7M@\x1fh\xb9M\xce\x1dh\xbbM\xe8\x03h\xbdM\xba\x03h\xbfh\xafh\xc1h\xb0ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(hE}\x94hK]\x94(hMhNehMhQhS\x85\x94R\x94(hS)}\x94hWhMsNt\x94bhNhQhS\x85\x94R\x94(hS)}\x94hWhNsNt\x94bub\x8c\x06config\x94}\x94(\x8c\x0fmolecular_pheno\x94}\x94(\x8c\nexpression\x94\x8c1/cluster/work/pausch/naveen/TWAS/TPM/{tissue}.tsv\x94\x8c\x08splicing\x94\x8c;/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/{tissue}.bed.gz\x94u\x8c\ncovariates\x94}\x94(\x8c\nexpression\x94\x8cG/cluster/work/pausch/xena/eQTL/cis_all/{tissue}/{tissue}_covariates.txt\x94\x8c\x08splicing\x94\x8cC/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/covariates_{tissue}.txt\x94u\x8c\x04vcfs\x94}\x94(\x8c\x04twas\x94\x8c>/cluster/work/pausch/naveen/ASE/VCF/chr_{chr}_beagle5.2.vcf.gz\x94\x8c\x04gwas\x94\x8cB/cluster/work/pausch/naveen/SNPDATA/TOSEQ/CHR{chr}/combined.vcf.gz\x94u\x8c\x03gtf\x94\x8cG/cluster/work/pausch/naveen/ASE/INPUTS/Bos_taurus.ARS-UCD1.2.104.gtf.gz\x94\x8c\x07OUT_DIR\x94\x8c /cluster/work/pausch/naveen/TWAS\x94\x8c\x0eold_gwas_files\x94\x8cR/cluster/work/pausch/xena/gwas/gcta/results/{{trait}}/BV/sire/CHR{chr}/result.mlma\x94\x8c\ngwas_files\x94\x8cM/cluster/work/pausch/xena/gwas/new/{inheri}/{{trait}}/CHR{chr}/result.mlma.gz\x94u\x8c\x04rule\x94\x8c\x17format_gwas_nonadditive\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c\'/cluster/home/nkadri/TWAS/SNAKES/spredi\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/cluster/home/nkadri/TWAS/SNAKES/spredi/format_gwas_nonadditive.py';
######## snakemake preamble end #########
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
