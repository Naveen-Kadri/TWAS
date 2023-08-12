'''
there are some parameters hard-coded into the training models
like the maf .. may be put it in the params ?
ref : https://github.com/hakyimlab/PredictDB-Tutorial
24-10-2022
'''
configfile: '/cluster/home/nkadri/TWAS/SNAKES/config.yaml'
OUT_DIR = '/cluster/work/pausch/naveen/TWAS'
# FILES
# 10 peer factors, 3 PCs, RIN, and age
# covariates = '/cluster/work/pausch/xena/eQTL/cis_all/{tissue}/{tissue}_covariates.txt'
# tpmg = '/cluster/work/pausch/naveen/TWAS/TPM/{tissue}.tsv'
gtf = "/cluster/work/pausch/naveen/ASE/INPUTS/Bos_taurus.ARS-UCD1.2.104.gtf.gz"
# vcf = "/cluster/work/pausch/naveen/ASE/VCF/chr_{chr}_beagle5.2.vcf.gz"

# new file with no / less stringent filter for HWE
# y

# WILDCARDS
chromosomes = range(1, 30)
tissues = ['testis', 'epi_h', 'vas_d']
#mphenos = ['splicing', 'expression']
mphenos = ['expression', 'splicing']


# TEST
#mphenos = ['expression']
#tissues = ['testis']
#chromosomes = [25]

# rule all:
#     input:
#         expand(OUT_DIR +
#                '/{mpheno}/{tissue}/CHR{chr}/covariance.txt', mpheno=mphenos, tissue=tissues, chr=chromosomes),
#         expand(OUT_DIR + '/{mpheno}/{tissue}/covariances.txt.gz',
#                mpheno=mphenos, tissue=tissues),
#         expand(OUT_DIR + '/{mpheno}/{tissue}/filtered_signif.db',
#                mpheno=mphenos, tissue=tissues)


rule expression_annotation:
    ''' We will have multiple splicing phenotypes per gene [splicing clusters]
    The program should be able to to associate the the name of the cluster to genomic cordinates
    So we need another annotation file
    '''
    input:
        gtf = config['gtf']
    output:
        outfile = OUT_DIR + '/expression/{tissue}/annotation.txt'
    script:
        "parse_gtf.py"

rule splicing_annotation:
    ''' The ids in the splicing phenotypes need to be associated to their genomic regions
    We have multiple phenotypes per gene!
    and should be done for each tissue separtely.. splicing might be tissue dependent
    '''
    input:
        gene_annotation = rules.expression_annotation.output.outfile,
        splicing = config['molecular_pheno']['splicing']
    output:
        outfile = OUT_DIR + '/splicing/{tissue}/annotation.txt'
    script:
        'splicing_annotation.py'


rule vcf2geno:
    ''' genotype format as required for the training of the model
    headers : VarID, sampleids...
    The info. on varID are in the snp_annotation file
    '''
    input:
        vcf = config['vcfs']['twas']
    output:
        genotypes = OUT_DIR + '/GENO/CHR{chr}/genotypes.txt',
        annotation = OUT_DIR + '/GENO/CHR{chr}/snp_annotation.txt'
    resources:
        mem_mb = 4000,
        walltime = "04:00"
    script:
        "vcf2geno.py"

rule split_and_transpose_mpheno:
    '''splicing data is huge.. better to split it... the file is also transposed as required for the Metaxscan
    In the splicing data the cluster names start with numbers ... these are used as rownames in R data frames and R does not like it
    So 'G' is added to the feature name .. same is done when the splicing_annotation rule.
    '''
    input:
        lambda wc: config['molecular_pheno'][wc.mpheno]
    output:
        outfile = OUT_DIR + '/{mpheno}/{tissue}/CHR{chr}/mpheno.txt'
    resources:
        mem_mb = 1000,
        walltime = "00:20"
    script:
        'split_and_transpose_mpheno.py'


rule format_covariates:
    '''
    Check the order of sample ids in the mpheno (chr25 across all chromosomes the id order is kept same) and the covariates data
    This is probably handled in elasic net script but to be on the safer side the order is checked.
    Set rownames to TRUE where rownames are the type of covariates
    '''
    input:
        lambda wildcards: config['covariates'][wildcards.mpheno],
        mpheno = OUT_DIR + '/{mpheno}/{tissue}/CHR25/mpheno.txt'
    output:
        outfile = OUT_DIR + '/{mpheno}/{tissue}/covariates.txt'
    resources:
        mem_mb = 2000,
        walltime = "00:20"
    script:
        'format_covariates.R'

'''
rsids are used in the code - trying with VarIDs as rsids.
Outputs 2 files to summary, 1 file  to covariances and 1 file  to weights directory
by default - code changed to pass file names explicity
'''

walltime = {
    "expression": "10:00",
    "splicing": "48:00"
}


def get_walltime(wc):
    return walltime[wc.mpheno]


rule train:
    input:
        snp_annotation = rules.vcf2geno.output.annotation,
        genotypes = rules.vcf2geno.output.genotypes,
        gene_annotation = OUT_DIR + '/{mpheno}/{tissue}/annotation.txt',
        covariates = OUT_DIR + '/{mpheno}/{tissue}/covariates.txt',
        mpheno = OUT_DIR + '/{mpheno}/{tissue}/CHR{chr}/mpheno.txt'
    output:
        model_summary_file = OUT_DIR +
        '/{mpheno}/{tissue}/CHR{chr}/summary.txt',
        weights_file = OUT_DIR + '/{mpheno}/{tissue}/CHR{chr}/weights.txt',
        tiss_chr_summ_f = OUT_DIR +
        '/{mpheno}/{tissue}/CHR{chr}/tiss_summ.txt',
        covariance_file = OUT_DIR +
        '/{mpheno}/{tissue}/CHR{chr}/covariance.txt'
    resources:
        mem_mb = 6000,
        walltime = get_walltime
    script:
        'gtex_tiss_chrom_training.R'


rule create_database:
    ''' All chromosomes are put together now.. can this also be done chromosomewise ?'''
    input:
        model_summary_files = expand(
            OUT_DIR + '/{{mpheno}}/{{tissue}}/CHR{chr}/summary.txt', chr=chromosomes),
        weights_files = expand(
            OUT_DIR + '/{{mpheno}}/{{tissue}}/CHR{chr}/weights.txt', chr=chromosomes),
        tiss_chr_summ_files = expand(
            OUT_DIR + '/{{mpheno}}/{{tissue}}/CHR{chr}/tiss_summ.txt', chr=chromosomes)
    params:
        nchr = max(chromosomes)
    output:
        database = OUT_DIR + '/{mpheno}/{tissue}/prediction_models.db'
    resources:
        mem_mb = 16000,
        walltime = "04:00"
    script:
        'create_database.R'

rule filter_database:
    ''' Parameters for filtering are hard-coded to zscore_pval < 0.05 and rho_avg > 0.1'''
    input:
        database = rules.create_database.output.database
    output:
        database = OUT_DIR + '/{mpheno}/{tissue}/filtered_signif.db'
    resources:
        mem_mb = 16000,
        walltime = "04:00"
    script:
        'filter_database.R'

rule cat_covariances:
    input:
        infiles = expand(
            OUT_DIR + '/{{mpheno}}/{{tissue}}/CHR{chr}/covariance.txt', chr=chromosomes)
    output:
        outfile = OUT_DIR + '/{mpheno}/{tissue}/covariances.txt.gz'
    shell:
        "awk '{{if(FNR==NR){{print}}else{{if(FNR>1){{print}} }}   }}' {input.infiles} |gzip > {output.outfile}"
