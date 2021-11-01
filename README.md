# Table of contents
1. [Snakemake](#1---snakemake)   
	1. [info](#info)  
	2. [execution](#execution)  
	3. [dependencies](#dependencies)  
2. [Scripts in python](#2---python)  
	1. [scripts](#scripts)  
	2. [dependencies](#dependencies)  
3. [Scripts in R](#3---r)  
	1. [scripts](#scripts)  
	2. [dependencies](#dependencies)  
4. [Code in Java](#4---java)  
	1. [msms (by Gregory Ewing and Joachim Hermisson)](#msms)  
5. [External codes](#5---external)  
6. [Config files](#6---config-files)  
	1. [cluster.json](#clusterjson)  
	2. [config.yaml](#configyaml)  
7. [Workflow](#7---workflow)  
	1. [Two populations](#two-populations)  
8. [Rapid tuto](#8---rapid-tuto)  
	1. [Clone the git depo](#clone-the-git-depo)  
	2. [Small edits](#small-edits)  
	3. [Run analysis](#run-analysis)  
9. [fasta2dadi](#9---fasta2dadi)  

# 1 - snakemake  
## info  
**The entire workflow is based on snakemake.**  
https://snakemake.readthedocs.io/en/stable/  
  
## execution  
Please adapt the pathway to your system.  
```  
snakemake --snakefile /shared/mfs/data/home/croux/softwares/DILS/bin/Snakefile_2pop -p -j 50 --configfile /shared/home/croux/scratch/myProject/config.yaml --cluster-config /shared/home/croux/scratch/myProject/cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time}"  
```  
  
## dependencies  
Needs:  
1 Snakefile  
1 config.yaml file  
1 cluster.json file  
  
# 2 - python  
**Executables are all found in the bin/ subdirectory. The pathway of subdirectory has to be indicated in the Snakefiles. This is the only required file modification**  
## scripts  
dadi2ABC.py  
submit_simulations_gof_2pop.py  
submit_simulations_2pop.py  
priorgen_2pop.py  
priorgen_gof_2pop.py  
mscalc_2pop_SFS.py  
  
## dependencies  
**some scripts uses pypy as python interpreter**    
import sys  
import os  
import time  
from random import shuffle  
from random import randint  
from random import choice  
from numpy.random import uniform  
from numpy.random import binomial  
from numpy.random import beta  
from numpy.random import randint  
from numpy.random import choice  
from numpy import median  
   
# 3 - R  
## scripts  
**uses Rscript from /usr/bin or elsewhere**  
gof_2pop.R  
get_parameters_2pop.R  
model_comp_2pop_allModels.R  
collaborative_plot.R  
estimates_1pop_best.R  
estimates_2pop.R  
PCA.R  
estimates_2pop_best.R  
   
## dependencies  
abcrf  
nnet  
tidyverse  
data.table  
ggplot2  
ggpubr  
viridis  
plotly  
FactoMineR  
   
# 4 - Java
## msms  
### info  
```  
java -jar msms3.2rc-b163.jar --help  
```  
We use msms as simulator (https://www.mabs.at/ewing/msms/download.shtml)  
   
# 5 - external  
**pandoc** (https://pandoc.org/index.html)  
The Pandoc call requires in this workflow that **pdflatex** is pre-installed.  
  
# 6 - config files  
## cluster.json  
This file contains informations for **Slurm** about the submited jobs, in particular, the required resources (CPU, memory, duration).  
```
{
    "__default__" :
    {
        "node" : 1,
        "ntasks" : 1,
        "n" : 1,
	"cpusPerTask" : 1,
	"memPerCpu" : 3000,
	"time" : "04:00:00"
    },
    "dadi2ABC" :
    {
	"cpusPerTask" : 1,
	"time" : "01:00:00",
	"memPerCpu" : 3000
    },
    "modelComparison" :
    {
	"cpusPerTask" : 8,
	"time" : "03:00:00",
	"memPerCpu" : 5000
    },
    "estimation_best_model" :
    {
	"cpusPerTask" : 8,
	"time" : "24:00:00",
	"memPerCpu" : 5000
    },
    "estimation_best_model_2" :
    {
	"cpusPerTask" : 8,
	"time" : "03:00:00",
	"memPerCpu" : 2500
    },
    "estimation_best_model_3" :
    {
	"cpusPerTask" : 8,
	"time" : "03:00:00",
	"memPerCpu" : 2500
    },
    "estimation_best_model_4" :
    {
	"cpusPerTask" : 8,
	"time" : "03:00:00",
	"memPerCpu" : 2500
    },
    "estimation_best_model_5" :
    {
	"cpusPerTask" : 8,
	"time" : "03:00:00",
	"memPerCpu" : 2500
    },
    "PCA_SS" :
    {
	"cpusPerTask" : 2,
	"time" : "03:00:00",
	"memPerCpu" : 7000
    }
}
``` 
  
## config.yaml  
Configuration file used by Snakemake to adapt the workflow to a particular analysis. Contains information such as species names, genomic region (coding, noncoding), prior boundaries, etc...  
```  
mail_address: camille.roux.1983@gmail.com  
infile: /home/croux/Programmes/DILS_sfs/test/all_SNPs_subsampled.txt  
nspecies: 2  
nameA: ama  
nameB: chi  
nA: 10  
nB: 10  
nameOutgroup: NA  
useSFS: 1  
useOutgroup: 0  
config_yaml: /home/croux/Programmes/DILS_sfs/test/test.yaml  
timeStamp: i8U91jkfd  
population_growth: constant  
modeBarrier: bimodal  
N_min: 0  
N_max: 500000  
Tsplit_min: 0  
Tsplit_max: 1750000  
M_min: 1  
M_max: 40  
```  
   
# 7 - workflow  
## two populations  
![DAG (directed acyclic graph)](https://raw.githubusercontent.com/popgenomics/snpDILS/master/dag.png)  
  

# 8 - Rapid tuto  
## clone the Git depo  
```
git clone https://github.com/popgenomics/snpDILS
```
  
## Small edits  
Then, edit a **Snakefile_Xpop** in order to adapt the variable **binpath**, corresponding to the pathway to the **bin** directory of the cloned git depository.  
Similarly, edit a **DILS_Xpop.sh** in order to adapt the same variable.  
The goal of these two small edits is to set the full pathways of all scripts + programs located into the **/bin** subdirectory of the git depository.  

And finally, don't forget to adapt pathways of the two input files within the **.yaml** file, at the following variables:
```
infile: /shared/ifbstor1/projects/metapop/snpDILS_v2/test/all_SNPs_subsampled.txt
config_yaml: /shared/ifbstor1/projects/metapop/snpDILS_v2/test/test.yaml
```  
  
## Run analysis (example)  
To perform the **test** analysis, first move into the __snpDILS_v2/test__ sub-directory and then:   
```
../bin/DILS_2pop.sh test.yaml
```
Of course, your analysis can be performed in the directory of your choice, just, adapt the pathways into the **.yaml** in order to let **DILS** finding the pathway of the input files.  

# 9 - fasta2dadi  
We add a small python script that can be used with **pypy** instead of **python3** (but python3 works fine also, just slower).  
This script takes as __input file__ a __fasta file__ that can be used in **DILS**, and produces a SNP data table that can be used by DaDi, Moments or snpDILS.  
```
Ing	Out	Allele1	ama	chi	flo	mal	ros	txn	vul	zel	Allele2	ama	chi	flo	mal	ros	txn	vul	zel	Gene	Position
TGC	TGC	G	20	20	20	20	20	20	20	19	T	0	0	0	0	0	0	0	1	Hmel219002_6	6
AGT	AGT	G	18	20	20	20	16	20	19	20	C	0	0	0	0	4	0	1	0	Hmel219002_6	24
ACC	ACC	C	19	20	20	20	20	20	20	20	T	1	0	0	0	0	0	0	0	Hmel219002_6	61
CCA	CCA	C	20	5	0	20	20	0	20	14	T	0	15	20	0	0	20	0	6	Hmel219002_6	62
AAG	AAG	A	20	20	20	20	19	20	20	20	T	0	0	0	0	1	0	0	0	Hmel219002_6	82
GCT	GCT	C	20	15	20	20	20	20	20	20	T	0	5	0	0	0	0	0	0	Hmel219002_6	84
CGG	CGG	G	20	20	20	20	11	20	7	14	A	0	0	0	0	9	0	13	6	Hmel219002_6	95
CCC	CCC	C	20	20	20	20	16	20	19	20	T	0	0	0	0	4	0	1	0	Hmel219002_6	107
TCA	TCA	C	20	18	20	20	20	20	20	20	T	0	2	0	0	0	0	0	0	Hmel219002_6	119
```
  
To run it:  
```
pypy3 test.py -i mytilus_renamed.fas -o coding_noOut_mussels.out -r coding
```

Or more help here:
```
pypy3 test.py -h
```


