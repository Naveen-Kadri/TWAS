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
