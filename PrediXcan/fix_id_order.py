# expr = '/cluster/work/pausch/naveen/TWAS/PREDICT/head/testis/CHR6/expression.txt'
# pheno = '/cluster/work/pausch/naveen/TWAS/PREDICT/head/phenotypes.txt'
# out_pheno_file = 'pheno.txt'
# out_expr_file = 'expr.txt'

pheno = snakemake.input.pheno
mpheno = snakemake.input.mpheno
out_pheno_file = snakemake.output.pheno
out_expr_file = snakemake.output.mpheno

out_pheno = open(out_pheno_file, "w")
out_expr = open(out_expr_file, "w")


phenod = {}
id_order = []
print('reading the phenotype file')
with open(pheno, "rt") as inf:
    for lnum, line in enumerate(inf):
        if lnum == 0:
            out_pheno.write(f'{line}')
        else:
            spl = line.rstrip().split()
            # some formatting error in the last line  of the file
            if len(spl) == 3:
                phenod[spl[0]] = line.rstrip()
            else:
                continue
print('reading the expr file')
with open(mpheno, "rt") as inf:
    for lnum, line in enumerate(inf):
        if lnum == 0:
            out_expr.write(f'{line}')
        else:
            myid = line.rstrip().split()[0]
            if myid in phenod:
                id_order.append(myid)
                out_expr.write(f'{line}')
print('writing the phenotype file')
for myid in id_order:
    out_pheno.write(f'{phenod[myid]}\n')
