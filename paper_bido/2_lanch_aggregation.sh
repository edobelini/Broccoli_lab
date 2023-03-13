#!/bin/bash
#SBATCH --job-name=aggr
#SBATCH --account bellini.edoardo
#SBATCH --mem=70GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=50  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=ALL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=bellini.edoardo@hsr.it
#SBATCH --error="aggr.err"
#SBATCH --output="aggr.out"

echo "my job strart now" > aggr.log;
date >> aggr.log;

export PATH=/beegfs/scratch/ric.broccoli/ric.broccoli/cellranger-6.1.2:$PATH

cellranger aggr --id=aggr \
                --csv=libraries.csv \
                --normalize=mapped \
                --localcores=50 \
                --localmem=70

date >> aggr.log;
echo "all done!!" >> aggr.log