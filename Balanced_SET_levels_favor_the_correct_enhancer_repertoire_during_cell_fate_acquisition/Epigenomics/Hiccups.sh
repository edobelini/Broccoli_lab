#!/bin/bash #This bash script was used to calculated Loops on mega maps of Neurons and NPCs D868
#SBATCH --job-name=hiccups
#SBATCH --account ric.broccoli
#SBATCH --mem=200GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=70  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=ALL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=zaghi.mattia@hsr.it
#SBATCH --error="hiccups.err"
#SBATCH --output="hiccups.out"

echo "my job strart now" > myjob.log;
date >> myjob.log; 

juiceDir="/beegfs/scratch/ric.broccoli/ric.broccoli/juicer1.5"
topDir=$(pwd)
genomeID="hg38"

${juiceDir}/CPU/common/juicer_tools hiccups --cpu --threads 70 -r 5000,10000 -k KR  -f .1,.1 -p 4,2 -i 7,5 -t 0.02,1.5,1.75,2 -d 20000,20000  ${topDir}/NPC_D868D/mega/aligned/inter_30.hic ${topDir}/NPC_D868D/mega/hiccups/up_in_dev_D868D ${topDir}/Loop_NPC_D868D_promoter_enhancers_contact.bedpe

${juiceDir}/CPU/common/juicer_tools hiccups --cpu --threads 70 -r 5000,10000 -k KR  -f .1,.1 -p 4,2 -i 7,5 -t 0.02,1.5,1.75,2 -d 20000,20000  ${topDir}/NPC_D868N/mega/aligned/inter_30.hic ${topDir}/NPC_D868N/mega/hiccups/up_in_dev_D868D ${topDir}/Loop_NPC_D868D_promoter_enhancers_contact.bedpe

${juiceDir}/CPU/common/juicer_tools hiccups --cpu --threads 70 -r 5000,10000 -k KR  -f .1,.1 -p 4,2 -i 7,5 -t 0.02,1.5,1.75,2 -d 20000,20000  ${topDir}/Neu_D868D/mega/aligned/inter_30.hic ${topDir}/Neu_D868D/mega/hiccups/up_in_dev_D868D ${topDir}/Loop_NPC_D868D_promoter_enhancers_contact.bedpe

${juiceDir}/CPU/common/juicer_tools hiccups --cpu --threads 70 -r 5000,10000 -k KR  -f .1,.1 -p 4,2 -i 7,5 -t 0.02,1.5,1.75,2 -d 20000,20000  ${topDir}/Neu_D868N/mega/aligned/inter_30.hic ${topDir}/Neu_D868N/mega/hiccups/up_in_dev_D868D ${topDir}/Loop_NPC_D868D_promoter_enhancers_contact.bedpe 

date >> myjob.log;
echo "all done!!" >> myjob.log