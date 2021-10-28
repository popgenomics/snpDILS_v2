#!/usr/bin/bash
## launch DILS for 2 populations
## the provided argument is for --configfile, expecting the yaml file
module unload snakemake
module load pypy/2.7-5.10.0
module load snakemake/6.5.0
module load r/3.6.3
module load python/2.7
#binpath='/shared/mfs/data/home/croux/softwares/DILS/bin'
binpath='/shared/ifbstor1/projects/metapop/snpDILS/bin'
snakemake --snakefile ${binpath}/Snakefile_2pop -p --cores 400 -j 400 --configfile ${1} --cluster-config ${binpath}/cluster_2pop.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time}" --latency-wait 10

