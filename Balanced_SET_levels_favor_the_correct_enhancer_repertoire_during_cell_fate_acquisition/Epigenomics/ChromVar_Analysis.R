library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)


#Calculate ATAC fragments counts inside regions of interest 

Neu_up_dev <- getPeaks("SETBP1_epigenomics/pipeline/Deseq2/Neu_up_dev.bed", sort_peaks = TRUE) #Bed files containing genomic regions to analyze 

bamfiles_NPC_Neu <-c("/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/NPCD868D1.bam", #Bam files to counts fragments 
                     "/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/NPCD868D2.bam",
                     "/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/NPCD868D3.bam",
                     "/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/NPCD868N1.bam",
                     "/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/NPCD868N2.bam",
                     "/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/NPCD868N3.bam",
                     "/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/NeuD868D1.bam",
                     "/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/NeuD868D2.bam",
                     "/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/NeuD868D3.bam",
                     "/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/NeuD868N1.bam",
                     "/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/NeuD868N2.bam",
                     "/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/NeuD868N3.bam")

fragment_counts_Neu_gain_all_peaks <- getCounts(bamfiles_NPC_Neu,Neu_up_dev, 
                                                      paired =  TRUE, 
                                                      by_rg = FALSE, 
                                                      format = "bam", 
                                                      colData = DataFrame(celltype =c("NPC_D868D",
                                                                                      "NPC_D868D",
                                                                                      "NPC_D868D",
                                                                                      "NPC_D868N",
                                                                                      "NPC_D868N",
                                                                                      "NPC_D868N",
                                                                                      "Neu_D868D",
                                                                                      "Neu_D868D",
                                                                                      "Neu_D868D",
                                                                                      "Neu_D868N",
                                                                                      "Neu_D868N",
                                                                                      "Neu_D868N")))


#analysis Neurons-----------------------------------------------------------------------


fragment_counts_Neu <- addGCBias(fragment_counts_Neu_gain_all_peaks,genome = BSgenome.Hsapiens.UCSC.hg38) %>% 
  filterSamples(min_depth = 1500, min_in_peaks = 0.15) %>% 
  filterPeaks()


motifs <- getJasparMotifs()

motif_ix <- matchMotifs(getJasparMotifs(), 
                        fragment_counts_Neu,
                        genome = BSgenome.Hsapiens.UCSC.hg38)

dev_Neu <- computeDeviations(object =fragment_counts_Neu, 
                         annotations = motif_ix, 
                         background_peaks = getBackgroundPeaks(fragment_counts_Neu ),
                         expectation = computeExpectations(fragment_counts_Neu ))

dev_Neu_data_frame <- deviationScores(dev_Neu)

head(dev_Neu)

#Differential Deviations calculations inside TFBS---------------------------------------------------------------------------------------

diff_acc_motifs_Neu <- differentialDeviations(dev_Neu, "celltype") %>% 
  as.data.frame() %>% 
  write_delim("SETBP1_epigenomics/pipeline/ChromVar/diff_acc_motifs_Neu_up_dev",
              delim="\t",col_names = T) 

diff_acc_motifs_Neu <- diff_acc_motifs_Neu %>% 
  rownames_to_column()

Neu_all_peaks_scores <- deviationScores(dev_Neu) %>% 
  as.data.frame() %>% 
  rownames_to_column()

Neu_all_peaks_statistical_significance <- inner_join(diff_acc_motifs_Neu,Neu_all_peaks_scores) %>% 
  write_delim("SETBP1_epigenomics/pipeline/ChromVar/Neu_all_peaks_statistical_significance_up_dev",
              delim="\t",col_names = T)

Neu_all_peaks_statistical_significance <- Neu_all_peaks_statistical_significance %>% 
  dplyr::rename(Motif=1,
                NPC_D868D_1=4,
                NPC_D868D_2=5,
                NPC_D868D_3=6,
                NPC_D868N_1=7,
                NPC_D868N_2=8,
                NPC_D868N_3=9,
                Neu_D868D_1=10,
                Neu_D868D_2=11,
                Neu_D868D_3=12,
                Neu_D868N_1=13,
                Neu_D868N_2=14,
                Neu_D868N_3=15) %>%
  rowwise() %>% 
  dplyr::mutate(NPC_D868D=mean(NPC_D868D_1+NPC_D868D_2+NPC_D868D_3)/3,
                NPC_D868N=mean(NPC_D868N_1+NPC_D868N_2+NPC_D868N_3)/3,
                Neu_D868D=mean(Neu_D868D_1+Neu_D868D_2+Neu_D868D_3)/3,
                Neu_D868N=mean(Neu_D868N_1+Neu_D868N_2+Neu_D868N_3)/3) %>% 
  dplyr::select(Motif,p_value,p_value_adjusted,NPC_D868D,NPC_D868N,NPC_D868D_1,NPC_D868D_2,NPC_D868D_3,NPC_D868N_1,NPC_D868N_2,NPC_D868N_3,Neu_D868D,Neu_D868N,Neu_D868D_1,Neu_D868D_2,Neu_D868D_3,Neu_D868N_1,Neu_D868N_2,Neu_D868N_3) %>% 
  write_delim("SETBP1_epigenomics/pipeline/ChromVar/Neu_all_peaks_statistical_significance_table",
              delim="\t",col_names = T)