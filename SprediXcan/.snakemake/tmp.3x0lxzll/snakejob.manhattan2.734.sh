#!/bin/sh
# properties = {"type": "single", "rule": "manhattan2", "local": false, "input": ["/cluster/work/pausch/naveen/TWAS/splicing/SprediXcan/nonadditive/tail/result_with_cords.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr1.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr2.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr3.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr4.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr5.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr6.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr7.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr8.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr9.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr10.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr11.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr12.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr13.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr14.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr15.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr16.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr17.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr18.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr19.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr20.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr21.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr22.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr23.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr24.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr25.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr26.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr27.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr28.txt.gz", "/cluster/work/pausch/xena/gwas/new/plots/tail/forplot_chr29.txt.gz"], "output": ["/cluster/work/pausch/naveen/TWAS/splicing/SprediXcan/nonadditive/tail/manhattan_splicing_tail_nonadditive.tiff"], "wildcards": {"mpheno": "splicing", "inheri": "nonadditive", "trait": "tail"}, "params": {"thresh_line": 6, "maf_thresh": 0.005}, "log": [], "threads": 1, "resources": {"mem_mb": 8000, "mem_mib": 7630, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>", "walltime": "01:00"}, "jobid": 734, "cluster": {}}
cd '/cluster/home/nkadri/TWAS/SNAKES/spredi' && /cluster/work/pausch/naveen/mambaforge/envs/snakemake/bin/python3.11 -m snakemake --snakefile '/cluster/home/nkadri/TWAS/SNAKES/spredi/spredi.smk' --target-jobs 'manhattan2:mpheno=splicing,inheri=nonadditive,trait=tail' --allowed-rules 'manhattan2' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=8000' 'mem_mib=7630' 'disk_mb=1000' 'disk_mib=954' --wait-for-files-file '/cluster/home/nkadri/TWAS/SNAKES/spredi/.snakemake/tmp.3x0lxzll/snakejob.manhattan2.734.sh.waitforfilesfile.txt' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'mtime' --skip-script-cleanup  --conda-frontend 'mamba' --use-envmodules  --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 10 --scheduler 'ilp' --scheduler-solver-path '/cluster/work/pausch/naveen/mambaforge/envs/snakemake/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && exit 0 || exit 1
