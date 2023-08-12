import os
'''
---IMPORTANT----
~/TWAS/MetXcan/software/metax/gwas/GWAS.py  -- MODIFIED NOT TO EXCLUDE SNPs with rsid
CHECK THE EFFECT OF CHANGING THE EFFECT ALLELE - MAKE SURE THE EFFECT ALLELE IS REALLY THE ALT. ALLELE



reference : https://github.com/hakyimlab/MetaXcan
example data (from the link in the repo) are here :/cluster/work/pausch/naveen/TWAS/example_data
we can already launch prediction once we have the summary stat from the GWAS
header for gwas file : SNP  A1  A2     FRQ    INFO    BETA      SE       P [what is INFO?]

'''

configfile: '/cluster/home/nkadri/TWAS/SNAKES/config.yaml'
include: '/cluster/home/nkadri/TWAS/SNAKES/train/train.smk'

# PROGRAMS
SPrediXcan = '/cluster/home/nkadri/TWAS/MetXcan/software/SPrediXcan.py '
# WILDACARDS
tissues = ['testis', 'epi_h', 'vas_d']
inheritances = ['additive', 'nonadditive']
traits = ['head', 'tail', 'vol', 'cnt', 'con',
          'mot', 'ntar', 'nrr', 'sbf', 'fbkall']
to_exclude = ['nrr', 'sbf']  # missing why?
traits = [trait for trait in traits if trait not in to_exclude]
print(f'analysis done for {traits}')

#tissues  =['epi_h', 'testis' ]
OUT_DIR = config['OUT_DIR']
rule all:
    input:
        expand( OUT_DIR + \
                '/{mpheno}/SprediXcan/{inheri}/{trait}/result_with_cords.txt.gz', mpheno=mphenos, inheri=inheritances, trait=traits),
        ##expand(OUT_DIR +
        #'/{mpheno}/SprediXcan/{inheri}/{trait}/{tissue}/manhattan_{mpheno}_{trait}_{tissue}_{inheri}.tiff', mpheno=mphenos, trait=traits, tissue=tissues, inheri=inheritances),
        expand( OUT_DIR + \
            '/{mpheno}/SprediXcan/{inheri}/{trait}/manhattan_{mpheno}_{trait}_{inheri}.tiff', mpheno=mphenos, trait=traits, inheri=inheritances),
        

rule format_gwas_additive:
    input:
        infiles = expand(config['gwas_files'],
                         chr=chromosomes, inheri='additive')
    output:
        outfile = OUT_DIR + '/GWAS/additive/{trait}/gwas.txt'
    script:
        'format_gwas_additive.py'

rule format_gwas_nonadditive:
    '''only positions common to recessive and  dominant files are considered'''
    input:
        recessive = expand(config['gwas_files'],
                           chr=chromosomes, inheri='recessive'),
        dominant = expand(config['gwas_files'],
                          chr=chromosomes, inheri='dominant')
    output:
        outfile = OUT_DIR + '/GWAS/nonadditive/{trait}/gwas.txt'
    resources:
        mem_mb = 8000,
        walltime = '04:00'
    script:
        'format_gwas_nonadditive.py'

rule SPrediXcan:
    ''' GWAS input can be a single file (gwas_file) or a directory (gwas_folder + gwas_file_pattern)
    Effect allele : how was the data coded ?
    '''
    input:
        db = rules.filter_database.output.database,
        gwas = OUT_DIR + '/GWAS/{inheri}/{trait}/gwas.txt',
        covariance = rules.cat_covariances.output.outfile
    output:
        result = OUT_DIR +
        '/{mpheno}/SprediXcan/{inheri}/{trait}/{tissue}/result.txt'
    resources:
        mem_mb = 36000,
        walltime = '02:00'
    shell:
        '''
        {SPrediXcan} \
        --model_db_path {input.db} \
        --covariance {input.covariance}  \
        --gwas_file {input.gwas} \
        --snp_column SNP \
        --effect_allele_column A2 \
        --non_effect_allele_column A1 \
        --beta_column BETA \
        --pvalue_column P \
        --output_file {output.result}
        '''
rule add_cordinates:
    input:
        result = rules.SPrediXcan.output.result,
        annotation = OUT_DIR + '/{mpheno}/{tissue}/annotation.txt'
    output:
        outfile = OUT_DIR + \
            '/{mpheno}/SprediXcan/{inheri}/{trait}/{tissue}/result_with_cords.txt.gz'
    script:
        'add_cordinates.py'

rule combine_tissues:
    input:
        infiles = expand (OUT_DIR + \
                          '/{{mpheno}}/SprediXcan/{{inheri}}/{{trait}}/{tissue}/result_with_cords.txt.gz', tissue=tissues)
    output:
        outfile = OUT_DIR + \
            '/{mpheno}/SprediXcan/{inheri}/{trait}/result_with_cords.txt.gz'
    resources:
        mem_mb=4000,
        walltime='00:20'
    params:
        tissues = tissues    
    script:
        "combine_tissues.py"
        
        
rule manhattan:
    '''
    gwas files have pvalues for both additive and non additive gwas
    '''
    input:
        twas = rules.add_cordinates.output.outfile,
        gwas = expand(
            '/cluster/work/pausch/xena/gwas/new/plots/{{trait}}/forplot_chr{chr}.txt.gz', chr=chromosomes)
    output:
        plot_file = OUT_DIR + \
            '/{mpheno}/SprediXcan/{inheri}/{trait}/{tissue}/manhattan_{mpheno}_{trait}_{tissue}_{inheri}.tiff'
    params:
        thresh_line = 6,
        maf_thresh = 0.005
    resources:
        mem_mb = 16000,
        walltime = '04:00'
    script:
        'manhattan_mirror.R'


rule manhattan2:
    '''
    gwas files have pvalues for both additive and non additive gwas
    '''
    input:
        twas = rules.combine_tissues.output.outfile,
        gwas = expand(
            '/cluster/work/pausch/xena/gwas/new/plots/{{trait}}/forplot_chr{chr}.txt.gz', chr=chromosomes)
    output:
        plot_file = OUT_DIR + \
            '/{mpheno}/SprediXcan/{inheri}/{trait}/manhattan_{mpheno}_{trait}_{inheri}.tiff'
    params:
        thresh_line = 6,
        maf_thresh = 0.005,
        tissues=tissues    
    resources:
        mem_mb = 8000,
        walltime = '01:00'
    script:
        'alltissues_mirror_manhattan.R'


        
