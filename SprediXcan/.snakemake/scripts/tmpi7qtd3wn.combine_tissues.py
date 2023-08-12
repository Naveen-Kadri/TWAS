
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/cluster/work/pausch/naveen/mambaforge/envs/snakemake/lib/python3.11/site-packages', '/cluster/home/nkadri/.cache/snakemake/snakemake/source-cache/runtime-cache/tmpnvmja_75/file/cluster/home/nkadri/TWAS/SNAKES/spredi', '/cluster/home/nkadri/TWAS/SNAKES/spredi']); import pickle; snakemake = pickle.loads(b"\x80\x04\x95\x10\n\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8ce/cluster/work/pausch/naveen/TWAS/splicing/SprediXcan/nonadditive/head/testis/result_with_cords.txt.gz\x94\x8cd/cluster/work/pausch/naveen/TWAS/splicing/SprediXcan/nonadditive/head/epi_h/result_with_cords.txt.gz\x94\x8cd/cluster/work/pausch/naveen/TWAS/splicing/SprediXcan/nonadditive/head/vas_d/result_with_cords.txt.gz\x94e}\x94(\x8c\x06_names\x94}\x94\x8c\x07infiles\x94K\x00K\x03\x86\x94s\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x14\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x1a)}\x94\x8c\x05_name\x94h\x14sNt\x94bh\x15h\x18h\x1a\x85\x94R\x94(h\x1a)}\x94h\x1eh\x15sNt\x94bh\x10h\x06\x8c\tNamedlist\x94\x93\x94)\x81\x94(h\nh\x0bh\x0ce}\x94(h\x0e}\x94h\x12]\x94(h\x14h\x15eh\x14h\x18h\x1a\x85\x94R\x94(h\x1a)}\x94h\x1eh\x14sNt\x94bh\x15h\x18h\x1a\x85\x94R\x94(h\x1a)}\x94h\x1eh\x15sNt\x94bubub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8c^/cluster/work/pausch/naveen/TWAS/splicing/SprediXcan/nonadditive/head/result_with_cords.txt.gz\x94a}\x94(h\x0e}\x94\x8c\x07outfile\x94K\x00N\x86\x94sh\x12]\x94(h\x14h\x15eh\x14h\x18h\x1a\x85\x94R\x94(h\x1a)}\x94h\x1eh\x14sNt\x94bh\x15h\x18h\x1a\x85\x94R\x94(h\x1a)}\x94h\x1eh\x15sNt\x94bh9h6ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94]\x94(\x8c\x06testis\x94\x8c\x05epi_h\x94\x8c\x05vas_d\x94ea}\x94(h\x0e}\x94\x8c\x07tissues\x94K\x00N\x86\x94sh\x12]\x94(h\x14h\x15eh\x14h\x18h\x1a\x85\x94R\x94(h\x1a)}\x94h\x1eh\x14sNt\x94bh\x15h\x18h\x1a\x85\x94R\x94(h\x1a)}\x94h\x1eh\x15sNt\x94bhNhHub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94(\x8c\x08splicing\x94\x8c\x0bnonadditive\x94\x8c\x04head\x94e}\x94(h\x0e}\x94(\x8c\x06mpheno\x94K\x00N\x86\x94\x8c\x06inheri\x94K\x01N\x86\x94\x8c\x05trait\x94K\x02N\x86\x94uh\x12]\x94(h\x14h\x15eh\x14h\x18h\x1a\x85\x94R\x94(h\x1a)}\x94h\x1eh\x14sNt\x94bh\x15h\x18h\x1a\x85\x94R\x94(h\x1a)}\x94h\x1eh\x15sNt\x94b\x8c\x06mpheno\x94h]hdh^hfh_ub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01M\xa0\x0fM\xe7\x0eM\xe8\x03M\xba\x03\x8c\x19/scratch/257474375.tmpdir\x94\x8c\x0500:20\x94e}\x94(h\x0e}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06mem_mb\x94K\x02N\x86\x94\x8c\x07mem_mib\x94K\x03N\x86\x94\x8c\x07disk_mb\x94K\x04N\x86\x94\x8c\x08disk_mib\x94K\x05N\x86\x94\x8c\x06tmpdir\x94K\x06N\x86\x94\x8c\x08walltime\x94K\x07N\x86\x94uh\x12]\x94(h\x14h\x15eh\x14h\x18h\x1a\x85\x94R\x94(h\x1a)}\x94h\x1eh\x14sNt\x94bh\x15h\x18h\x1a\x85\x94R\x94(h\x1a)}\x94h\x1eh\x15sNt\x94bh{K\x01h}K\x01h\x7fM\xa0\x0fh\x81M\xe7\x0eh\x83M\xe8\x03h\x85M\xba\x03h\x87hwh\x89hxub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0e}\x94h\x12]\x94(h\x14h\x15eh\x14h\x18h\x1a\x85\x94R\x94(h\x1a)}\x94h\x1eh\x14sNt\x94bh\x15h\x18h\x1a\x85\x94R\x94(h\x1a)}\x94h\x1eh\x15sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x0fmolecular_pheno\x94}\x94(\x8c\nexpression\x94\x8c1/cluster/work/pausch/naveen/TWAS/TPM/{tissue}.tsv\x94\x8c\x08splicing\x94\x8c;/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/{tissue}.bed.gz\x94u\x8c\ncovariates\x94}\x94(\x8c\nexpression\x94\x8cG/cluster/work/pausch/xena/eQTL/cis_all/{tissue}/{tissue}_covariates.txt\x94\x8c\x08splicing\x94\x8cC/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/covariates_{tissue}.txt\x94u\x8c\x04vcfs\x94}\x94(\x8c\x04twas\x94\x8c>/cluster/work/pausch/naveen/ASE/VCF/chr_{chr}_beagle5.2.vcf.gz\x94\x8c\x04gwas\x94\x8cB/cluster/work/pausch/naveen/SNPDATA/TOSEQ/CHR{chr}/combined.vcf.gz\x94u\x8c\x03gtf\x94\x8cG/cluster/work/pausch/naveen/ASE/INPUTS/Bos_taurus.ARS-UCD1.2.104.gtf.gz\x94\x8c\x07OUT_DIR\x94\x8c /cluster/work/pausch/naveen/TWAS\x94\x8c\x0eold_gwas_files\x94\x8cR/cluster/work/pausch/xena/gwas/gcta/results/{{trait}}/BV/sire/CHR{chr}/result.mlma\x94\x8c\ngwas_files\x94\x8cM/cluster/work/pausch/xena/gwas/new/{inheri}/{{trait}}/CHR{chr}/result.mlma.gz\x94u\x8c\x04rule\x94\x8c\x0fcombine_tissues\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c'/cluster/home/nkadri/TWAS/SNAKES/spredi\x94ub."); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/cluster/home/nkadri/TWAS/SNAKES/spredi/combine_tissues.py';
######## snakemake preamble end #########
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
