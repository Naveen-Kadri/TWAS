
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/cluster/work/pausch/naveen/mambaforge/envs/snakemake/lib/python3.11/site-packages', '/cluster/home/nkadri/.cache/snakemake/snakemake/source-cache/runtime-cache/tmpcw7774a9/file/cluster/home/nkadri/TWAS/SNAKES/train', '/cluster/home/nkadri/TWAS/SNAKES/train']); import pickle; snakemake = pickle.loads(b"\x80\x04\x95U\x08\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8c8/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/vas_d.bed.gz\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x10\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x16)}\x94\x8c\x05_name\x94h\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8c@/cluster/work/pausch/naveen/TWAS/splicing/vas_d/CHR19/mpheno.txt\x94a}\x94(h\x0c}\x94\x8c\x07outfile\x94K\x00N\x86\x94sh\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bh'h$ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94(\x8c\x08splicing\x94\x8c\x05vas_d\x94\x8c\x0219\x94e}\x94(h\x0c}\x94(\x8c\x06mpheno\x94K\x00N\x86\x94\x8c\x06tissue\x94K\x01N\x86\x94\x8c\x03chr\x94K\x02N\x86\x94uh\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94b\x8c\x06mpheno\x94hEhLhF\x8c\x03chr\x94hGub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01M\xe8\x03M\xba\x03M\xe8\x03M\xba\x03\x8c\x19/scratch/257441038.tmpdir\x94\x8c\x0500:20\x94e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06mem_mb\x94K\x02N\x86\x94\x8c\x07mem_mib\x94K\x03N\x86\x94\x8c\x07disk_mb\x94K\x04N\x86\x94\x8c\x08disk_mib\x94K\x05N\x86\x94\x8c\x06tmpdir\x94K\x06N\x86\x94\x8c\x08walltime\x94K\x07N\x86\x94uh\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bhdK\x01hfK\x01hhM\xe8\x03hjM\xba\x03hlM\xe8\x03hnM\xba\x03hph`hrhaub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x0fmolecular_pheno\x94}\x94(\x8c\nexpression\x94\x8c1/cluster/work/pausch/naveen/TWAS/TPM/{tissue}.tsv\x94\x8c\x08splicing\x94\x8c;/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/{tissue}.bed.gz\x94u\x8c\ncovariates\x94}\x94(\x8c\nexpression\x94\x8cG/cluster/work/pausch/xena/eQTL/cis_all/{tissue}/{tissue}_covariates.txt\x94\x8c\x08splicing\x94\x8cC/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/covariates_{tissue}.txt\x94u\x8c\x04vcfs\x94}\x94(\x8c\x04twas\x94\x8c>/cluster/work/pausch/naveen/ASE/VCF/chr_{chr}_beagle5.2.vcf.gz\x94\x8c\x04gwas\x94\x8cB/cluster/work/pausch/naveen/SNPDATA/TOSEQ/CHR{chr}/combined.vcf.gz\x94u\x8c\x03gtf\x94\x8cG/cluster/work/pausch/naveen/ASE/INPUTS/Bos_taurus.ARS-UCD1.2.104.gtf.gz\x94\x8c\x07OUT_DIR\x94\x8c /cluster/work/pausch/naveen/TWAS\x94\x8c\x0eold_gwas_files\x94\x8cR/cluster/work/pausch/xena/gwas/gcta/results/{{trait}}/BV/sire/CHR{chr}/result.mlma\x94\x8c\ngwas_files\x94\x8cM/cluster/work/pausch/xena/gwas/new/{inheri}/{{trait}}/CHR{chr}/result.mlma.gz\x94u\x8c\x04rule\x94\x8c\x1asplit_and_transpose_mpheno\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c&/cluster/home/nkadri/TWAS/SNAKES/train\x94ub."); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/cluster/home/nkadri/TWAS/SNAKES/train/split_and_transpose_mpheno.py';
######## snakemake preamble end #########
import gzip
from collections import defaultdict

infile = snakemake.input[0]
outfile = snakemake.output.outfile
mychr = snakemake.wildcards.chr
mpheno = snakemake.wildcards.mpheno

#infile = '/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/testis.bed.gz'
#outfile = 'todel.txt'
#mychr = '25'

info = defaultdict(dict)

if 'gz' in infile:
    inf = gzip.open(infile, "rt")
else:
    inf = open(infile)

features = []
for lnum, line in enumerate(inf):
    spl = line.rstrip().split()
    if lnum == 0:
        header = spl
        print(spl[:10])
        ids = spl[6:]
    else:
        if spl[0] == mychr:
            phenos = spl[6:]
            features.append(spl[3])
            for myid, mypheno in zip(ids, phenos):
                info[spl[3]][myid] = mypheno


inf.close()
print(f'number of features:{len(features)}')


##transpose and print
with open(outfile, 'w') as out:
    if mpheno == 'splicing':
        # in our splicing data the feature names start with numbers -- this is used as rownames for the R data frame and R does not like it
        # also R does not like '.' in the names so that is replaced with '_'
        Gfeatures = ["G"+feature.replace(":", "_") for feature in features]
        tw = "\t".join(Gfeatures)
    else:
        tw = "\t".join(features)
    out.write(f'{tw}\n')
    for myid in ids:  # important to keep the id order same across files
        out.write(f'{myid}')
        for myfeature in features:
            out.write(f'\t{info[myfeature][myid]}')
        out.write('\n')
