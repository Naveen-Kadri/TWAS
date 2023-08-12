
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/cluster/work/pausch/naveen/mambaforge/envs/snakemake/lib/python3.11/site-packages', '/cluster/home/nkadri/.cache/snakemake/snakemake/source-cache/runtime-cache/tmp340t8ny4/file/cluster/home/nkadri/TWAS/SNAKES/spredi', '/cluster/home/nkadri/TWAS/SNAKES/spredi']); import pickle; snakemake = pickle.loads(b"\x80\x04\x95\xf9\x08\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8cU/cluster/work/pausch/naveen/TWAS/splicing/SprediXcan/nonadditive/vol/vas_d/result.txt\x94\x8c>/cluster/work/pausch/naveen/TWAS/splicing/vas_d/annotation.txt\x94e}\x94(\x8c\x06_names\x94}\x94(\x8c\x06result\x94K\x00N\x86\x94\x8c\nannotation\x94K\x01N\x86\x94u\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x15\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x1b)}\x94\x8c\x05_name\x94h\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bh\x0fh\nh\x11h\x0bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8cc/cluster/work/pausch/naveen/TWAS/splicing/SprediXcan/nonadditive/vol/vas_d/result_with_cords.txt.gz\x94a}\x94(h\r}\x94\x8c\x07outfile\x94K\x00N\x86\x94sh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bh,h)ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\r}\x94h\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94(\x8c\x08splicing\x94\x8c\x0bnonadditive\x94\x8c\x03vol\x94\x8c\x05vas_d\x94e}\x94(h\r}\x94(\x8c\x06mpheno\x94K\x00N\x86\x94\x8c\x06inheri\x94K\x01N\x86\x94\x8c\x05trait\x94K\x02N\x86\x94\x8c\x06tissue\x94K\x03N\x86\x94uh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94b\x8c\x06mpheno\x94hJhRhKhThLhVhMub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01M\xe8\x03M\xba\x03M\xe8\x03M\xba\x03\x8c\x19/scratch/257492675.tmpdir\x94e}\x94(h\r}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06mem_mb\x94K\x02N\x86\x94\x8c\x07mem_mib\x94K\x03N\x86\x94\x8c\x07disk_mb\x94K\x04N\x86\x94\x8c\x08disk_mib\x94K\x05N\x86\x94\x8c\x06tmpdir\x94K\x06N\x86\x94uh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bhjK\x01hlK\x01hnM\xe8\x03hpM\xba\x03hrM\xe8\x03htM\xba\x03hvhgub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\r}\x94h\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x0fmolecular_pheno\x94}\x94(\x8c\nexpression\x94\x8c1/cluster/work/pausch/naveen/TWAS/TPM/{tissue}.tsv\x94\x8c\x08splicing\x94\x8c;/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/{tissue}.bed.gz\x94u\x8c\ncovariates\x94}\x94(\x8c\nexpression\x94\x8cG/cluster/work/pausch/xena/eQTL/cis_all/{tissue}/{tissue}_covariates.txt\x94\x8c\x08splicing\x94\x8cC/cluster/work/pausch/naveen/TWAS/EXC_RATIOS/covariates_{tissue}.txt\x94u\x8c\x04vcfs\x94}\x94(\x8c\x04twas\x94\x8c>/cluster/work/pausch/naveen/ASE/VCF/chr_{chr}_beagle5.2.vcf.gz\x94\x8c\x04gwas\x94\x8cB/cluster/work/pausch/naveen/SNPDATA/TOSEQ/CHR{chr}/combined.vcf.gz\x94u\x8c\x03gtf\x94\x8cG/cluster/work/pausch/naveen/ASE/INPUTS/Bos_taurus.ARS-UCD1.2.104.gtf.gz\x94\x8c\x07OUT_DIR\x94\x8c /cluster/work/pausch/naveen/TWAS\x94\x8c\x0eold_gwas_files\x94\x8cR/cluster/work/pausch/xena/gwas/gcta/results/{{trait}}/BV/sire/CHR{chr}/result.mlma\x94\x8c\ngwas_files\x94\x8cM/cluster/work/pausch/xena/gwas/new/{inheri}/{{trait}}/CHR{chr}/result.mlma.gz\x94u\x8c\x04rule\x94\x8c\x0eadd_cordinates\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c'/cluster/home/nkadri/TWAS/SNAKES/spredi\x94ub."); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/cluster/home/nkadri/TWAS/SNAKES/spredi/add_cordinates.py';
######## snakemake preamble end #########
import gzip


#annotation = '/cluster/work/pausch/naveen/TWAS/splicing/testis/annotation.txt'
#result = '/cluster/work/pausch/naveen/TWAS/splicing/SprediXcan/additive/tail/testis/result.txt'
#outfile = 'todel.gz'

annotation = snakemake.input.annotation
result = snakemake.input.result
outfile = snakemake.output.outfile

out = gzip.open(outfile, 'wt')

# chr	gene_id	gene_name	start	end	gene_type
cords = {}
with open(annotation) as inf:
    for lnum, line in enumerate(inf):
        spl = line.rstrip().split()
        cords[spl[1]] = '\t'.join([spl[0], spl[3], spl[4]])

# gene,gene_name,zscore,effect_size,pvalue,var_g,pred_perf_r2,pred_perf_pval,pred_perf_qval,n_snps_used,n_snps_in_cov,n_snps_in_model
with open(result) as inf:
    for lnum, line in enumerate(inf):
        spl = line.rstrip().split(",")
        if lnum == 0:
            header = '\t'.join(spl + ['chr', 'start', 'end'])
            out.write(f'{header}\n')
        else:
            tw = '\t'.join(spl + [cords[spl[0]]])
            out.write(f'{tw}\n')
