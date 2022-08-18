#!/bin/bash #This bash script was used to calculated Contact Domains on mega maps of Neurons and NPCs D868
#SBATCH --job-name=arrowhead
#SBATCH --account ric.broccoli
#SBATCH --mem=200GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=40  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=ALL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=zaghi.mattia@hsr.it
#SBATCH --error="arrowhead.err"
#SBATCH --output="arrowhead.out"

echo "my job strart now" > myjob.log;
date >> myjob.log; 

juiceDir="/beegfs/scratch/ric.broccoli/ric.broccoli/juicer1.5"
topDir=$(pwd)
genomeID="hg38"

${juiceDir}/CPU/common/juicer_tools arrowhead --threads 70   ${topDir}/NPC_D868D/mega/aligned/inter_30.hic ${topDir}/NPC_D868D/mega/arrowhead 

${juiceDir}/CPU/common/juicer_tools arrowhead --threads 70   ${topDir}/NPC_D868N/mega/aligned/inter_30.hic ${topDir}/NPC_D868N/mega/arrowhead 

${juiceDir}/CPU/common/juicer_tools arrowhead --threads 70   ${topDir}/Neu_D868D/mega/aligned/inter_30.hic ${topDir}/Neu_D868D/mega/arrowhead 

${juiceDir}/CPU/common/juicer_tools arrowhead --threads 70   ${topDir}/Neu_D868N/mega/aligned/inter_30.hic ${topDir}/Neu_D868N/mega/arrowhead

date >> myjob.log;
echo "all done!!" >> myjob.log
