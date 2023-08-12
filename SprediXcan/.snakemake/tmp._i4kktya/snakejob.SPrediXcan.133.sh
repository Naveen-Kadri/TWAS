#!/bin/sh
# properties = {"type": "single", "rule": "SPrediXcan", "local": false, "input": ["/cluster/work/pausch/naveen/TWAS/expression/testis/filtered_signif.db", "/cluster/work/pausch/naveen/TWAS/GWAS/additive/vol/gwas.txt", "/cluster/work/pausch/naveen/TWAS/expression/testis/covariances.txt.gz"], "output": ["/cluster/work/pausch/naveen/TWAS/expression/SprediXcan/additive/vol/testis/result.txt"], "wildcards": {"mpheno": "expression", "inheri": "additive", "trait": "vol", "tissue": "testis"}, "params": {}, "log": [], "threads": 1, "resources": {"mem_mb": 36000, "mem_mib": 34333, "disk_mb": 2514, "disk_mib": 2398, "tmpdir": "<TBD>", "walltime": "02:00"}, "jobid": 133, "cluster": {}}
cd '/cluster/home/nkadri/TWAS/SNAKES/spredi' && /cluster/work/pausch/naveen/mambaforge/envs/snakemake/bin/python3.11 -m snakemake --snakefile '/cluster/home/nkadri/TWAS/SNAKES/spredi/spredi.smk' --target-jobs 'SPrediXcan:mpheno=expression,inheri=additive,trait=vol,tissue=testis' --allowed-rules 'SPrediXcan' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=36000' 'mem_mib=34333' 'disk_mb=2514' 'disk_mib=2398' --wait-for-files '/cluster/home/nkadri/TWAS/SNAKES/spredi/.snakemake/tmp._i4kktya' '/cluster/work/pausch/naveen/TWAS/expression/testis/filtered_signif.db' '/cluster/work/pausch/naveen/TWAS/GWAS/additive/vol/gwas.txt' '/cluster/work/pausch/naveen/TWAS/expression/testis/covariances.txt.gz' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'mtime' --skip-script-cleanup  --conda-frontend 'mamba' --use-envmodules  --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 10 --scheduler 'ilp' --scheduler-solver-path '/cluster/work/pausch/naveen/mambaforge/envs/snakemake/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && exit 0 || exit 1
