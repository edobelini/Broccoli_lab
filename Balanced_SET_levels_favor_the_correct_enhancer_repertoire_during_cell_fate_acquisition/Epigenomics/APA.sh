#!/bin/bash # This script was used to perform APA analysis on NPCs D868D loops
#SBATCH --job-name=APA
#SBATCH --account ric.broccoli
#SBATCH --mem=100GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=50  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=ALL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=zaghi.mattia@hsr.it
#SBATCH --error="APA.err"
#SBATCH --output="APA.out"

echo "my job strart now" > myjob.log;
date >> myjob.log; 

juiceDir="/beegfs/scratch/ric.broccoli/ric.broccoli/juicer1.5"
topDir=$(pwd)
genomeID="hg38"

${juiceDir}/CPU/common/juicer_tools apa -r 5000,10000 /beegfs/scratch/ric.broccoli/ric.broccoli/HiC/NPC_D868D/mega/aligned/inter_30.hic /beegfs/scratch/ric.broccoli/ric.broccoli/HiC/NPC_D868D/hiccups/merged_loops.bedpe /beegfs/scratch/ric.broccoli/ric.broccoli/HiC/NPC_D868D/APA

${juiceDir}/CPU/common/juicer_tools apa -r 5000,10000 /beegfs/scratch/ric.broccoli/ric.broccoli/HiC/NPC_D868N/mega/aligned/inter_30.hic /beegfs/scratch/ric.broccoli/ric.broccoli/HiC/NPC_D868D/hiccups/merged_loops.bedpe /beegfs/scratch/ric.broccoli/ric.broccoli/HiC/NPC_D868N/APA

${juiceDir}/CPU/common/juicer_tools apa -r 5000,10000 /beegfs/scratch/ric.broccoli/ric.broccoli/HiC/NPC_D868D/mega/aligned/inter_30.hic /beegfs/scratch/ric.broccoli/ric.broccoli/HiC/Loop_NPC_D868D_promoter_enhancers_contact.bedpe /beegfs/scratch/ric.broccoli/ric.broccoli/HiC/NPC_D868D/APA/decreased_loop

${juiceDir}/CPU/common/juicer_tools apa -r 5000,10000 /beegfs/scratch/ric.broccoli/ric.broccoli/HiC/NPC_D868N/mega/aligned/inter_30.hic /beegfs/scratch/ric.broccoli/ric.broccoli/HiC/Loop_NPC_D868D_promoter_enhancers_contact.bedpe /beegfs/scratch/ric.broccoli/ric.broccoli/HiC/NPC_D868N/APA/decreased_loop
