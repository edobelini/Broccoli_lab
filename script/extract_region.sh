for i  in sox10 sox15 sox17 sox21 sox2 sox3 sox4 sox6 sox7 sox9 oct4-sox17 oct4-sox2 Sox_fam
do
	CMD="/home/edoardo/homer/bin/findMotifsGenome.pl SES.bed hg19 /home/edoardo/sox2/TF_enrichment/Cut_Tag/SES -find /home/edoardo/homer/motifs/$i.motif > /home/edoardo/sox2/TF_enrichment/Cut_Tag/SES/$i.txt"
	echo $CMD
	eval $CMD
	CMD2="/home/edoardo/homer/bin/findMotifsGenome.pl SOX2.bed hg19 /home/edoardo/sox2/TF_enrichment/Cut_Tag/SOX2 -find /home/edoardo/homer/motifs/$i.motif > /home/edoardo/sox2/TF_enrichment/Cut_Tag/SOX2/$i.txt"
	echo $CMD2
	eval $CMD2

 done 
