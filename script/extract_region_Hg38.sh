mkdir  /home/edoardo/sox2/TF_enrichment/Cut_Tag/New_targets_hg38/sox_motif
mkdir  /home/edoardo/sox2/TF_enrichment/Cut_Tag/New_targets_cluster1_hg38/sox_motif
mkdir  /home/edoardo/sox2/TF_enrichment/Cut_Tag/New_targets_cluster1_hg38/sox_motif
for i  in sox10 sox15 sox17 sox21 sox2 sox3 sox4 sox6 sox7 sox9 oct4-sox17 oct4-sox2 Sox_fam
do
	mkdir  /home/edoardo/sox2/TF_enrichment/Cut_Tag/New_targets_hg38/sox_motif
	CMD="/home/edoardo/homer/bin/findMotifsGenome.pl New_targets.bed hg38 /home/edoardo/sox2/TF_enrichment/Cut_Tag/New_targets_hg38 -find /home/edoardo/homer/motifs/$i.motif > /home/edoardo/sox2/TF_enrichment/Cut_Tag/New_targets_hg38/sox_motif/$i.txt"
	echo $CMD
	eval $CMD

	mkdir  /home/edoardo/sox2/TF_enrichment/Cut_Tag/New_targets_cluster1_hg38/sox_motif
	CMD2="/home/edoardo/homer/bin/findMotifsGenome.pl New_targets_cluster1.bed hg38 /home/edoardo/sox2/TF_enrichment/Cut_Tag/New_targets_cluster1_hg38 -find /home/edoardo/homer/motifs/$i.motif > /home/edoardo/sox2/TF_enrichment/Cut_Tag/New_targets_cluster1_hg38/sox_motif/$i.txt"
	echo $CMD2
	eval $CMD2

	mkdir  /home/edoardo/sox2/TF_enrichment/Cut_Tag/New_targets_cluster2_hg38/sox_motif
	CMD3="/home/edoardo/homer/bin/findMotifsGenome.pl New_targets_cluster2.bed hg38 /home/edoardo/sox2/TF_enrichment/Cut_Tag/New_targets_cluster2_hg38 -find /home/edoardo/homer/motifs/$i.motif > /home/edoardo/sox2/TF_enrichment/Cut_Tag/New_targets_cluster2_hg38/sox_motif/$i.txt"
	echo $CMD3
	eval $CMD3

done
