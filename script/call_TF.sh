for i in SES SOX2
do
	cmd="/home/edoardo/homer/bin/findMotifsGenome.pl $i.bed hg19 $i/ -size 200 -mask -p 20"
	echo $cmd
	eval $cmd
done
