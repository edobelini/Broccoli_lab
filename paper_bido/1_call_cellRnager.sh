#!/bin/bash
#SBATCH --job-name=Cellranger
#SBATCH --account bellini.edoardo
#SBATCH --mem=70GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=50  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=ALL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=bellini.edoardo@hsr.it
#SBATCH --error="cellRanger.err"
#SBATCH --output="cellRanger.out"

echo "my job strart now" > cellRanger.log;
date >> cellRanger.log;

export PATH=/beegfs/scratch/ric.broccoli/ric.broccoli/cellranger-6.1.2:$PATH

for i in WT SNCA SNCA_IL10 IL10
    cellranger count --id=$i \
                 --transcriptome=/beegfs/scratch/ric.broccoli/ric.broccoli/mm10_SNCA_IL10 \
                 --fastqs=/beegfs/scratch/ric.broccoli/ric.broccoli/runs/220414_A00626_0432_AHCFNNDMXY_GS_Index/BroccoliV_1709_synaptosome_scRNA/$i \
                 --localcores=50 \
                 --localmem=70
done

date >> cellRanger.log;
echo "all done!!" >> cellRanger.log