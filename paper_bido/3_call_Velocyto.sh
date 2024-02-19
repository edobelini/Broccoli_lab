#!/bin/bash
#SBATCH --job-name=Velocyto
#SBATCH --account bellini.edoardo
#SBATCH --mem=50GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=70  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=ALL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=bellini.edoardo@hsr.it
#SBATCH --error="vel.err"
#SBATCH --output="vel.out"

echo "my job strart now" > vel.log;
date >> vel.log;
. /opt/common/tools/ric.cosr/miniconda3/bin/activate;

conda activate velocyto;
for i in WT SNCA SNCA_IL10 IL10
do
    velocyto run10x  /beegfs/scratch/ric.broccoli/ric.broccoli/sc_Bido_Il10/$i /beegfs/scratch/ric.broccoli/ric.broccoli/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf
done
date >> vel.log;
echo "all done!!" >> vel.log