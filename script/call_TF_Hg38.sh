for i in New_targets_cluster1 New_targets_cluster2 New_targets
do
	cmd="/home/edoardo/homer/bin/findMotifsGenome.pl $i.bed hg38 "${i}_hg38/" -size 200 -mask -p 20"
	echo $cmd
	eval $cmd
done
