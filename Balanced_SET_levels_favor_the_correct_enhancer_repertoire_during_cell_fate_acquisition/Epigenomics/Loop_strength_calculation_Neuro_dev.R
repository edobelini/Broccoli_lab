library(tidyverse)
library(dplyr)
library(diffloop)
library(LSD)

#Loop strength calculation NPCs D868D vs Neu D868D (Neu D868D loops lists)

Loop_Neu_D868D_5kb <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/Neu_D868D/hiccups/enriched_pixels_5000.bedpe",
                                 delim="\t",col_names = T)

Loop_Neu_D868D_10kb <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/Neu_D868D/hiccups/enriched_pixels_10000.bedpe",
                                  delim="\t",col_names = T)


Loop_NPC_D868D_inNeuD868D_5kb <- read_delim("SETBP1_epigenomics/HiC/Neu_D868D_vs_NPC_D868D/hiccupsdiff/file2/requested_list_5000.bedpe",
                                            delim="\t",col_names = T)

Loop_NPC_D868D_inNeuD868D_10kb <- read_delim("SETBP1_epigenomics/HiC/Neu_D868D_vs_NPC_D868D/hiccupsdiff/file2/requested_list_10000.bedpe",
                                             delim="\t",col_names = T)


Loop_Neu_D868D_merged <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/Neu_D868D/hiccups/merged_loops.bedpe",
                                    delim="\t",col_names = T)

#join different resolution:


Loop_Neu_D868D_merged <- Loop_Neu_D868D_merged[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6)%>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2)


Loop_Neu_D868D_5kb <- Loop_Neu_D868D_5kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_Neu_D868D=12,
                expectedBL_Neu_D868D=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_Neu_D868D,
                expectedBL_Neu_D868D)

Loop_NPC_D868D_inNeuD868D_5kb <- Loop_NPC_D868D_inNeuD868D_5kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_NPC_D868D=12,
                expectedBL_NPC_D868D=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_NPC_D868D,
                expectedBL_NPC_D868D)

Loop_Neu_D868D_vs_NPC_D868D_5kb <- inner_join(Loop_Neu_D868D_5kb,Loop_NPC_D868D_inNeuD868D_5kb)


Loop_Neu_D868D_10kb <- Loop_Neu_D868D_10kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_Neu_D868D=12,
                expectedBL_Neu_D868D=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_Neu_D868D,
                expectedBL_Neu_D868D)

Loop_NPC_D868D_inNeuD868D_10kb <- Loop_NPC_D868D_inNeuD868D_10kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_NPC_D868D=12,
                expectedBL_NPC_D868D=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_NPC_D868D,
                expectedBL_NPC_D868D)

Loop_Neu_D868D_vs_NPC_D868D_10kb <- inner_join(Loop_Neu_D868D_10kb,Loop_NPC_D868D_inNeuD868D_10kb)

Loop_Neu_D868D_vs_NPCD868D_5kb_merged <- inner_join(Loop_Neu_D868D_merged,Loop_Neu_D868D_vs_NPC_D868D_5kb)

Loop_Neu_D868D_vs_NPCD868D_10kb_merged <- inner_join(Loop_Neu_D868D_merged,Loop_Neu_D868D_vs_NPC_D868D_10kb)

Loop_Neu_D868D_vs_NPC_D868D_merged_Neu_Loops <- rbind.data.frame(Loop_Neu_D868D_vs_NPCD868D_5kb_merged,Loop_Neu_D868D_vs_NPCD868D_10kb_merged) %>% 
  dplyr::mutate(intensity_NPC_D868D=observed_NPC_D868D/expectedBL_NPC_D868D) %>% 
  dplyr::mutate(intensity_Neu_D868D=observed_Neu_D868D/expectedBL_Neu_D868D) %>% 
  dplyr::mutate(fold=foldchange(intensity_Neu_D868D,intensity_NPC_D868D))


#Loop strength calculation NPCs D868D vs Neu D868D (NPCs D868D loops lists)


Loop_NPC_D868D_5kb <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/NPC_D868D/hiccups/enriched_pixels_5000.bedpe",
                                 delim="\t",col_names = T)

Loop_NPC_D868D_10kb <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/NPC_D868D/hiccups/enriched_pixels_10000.bedpe",
                                  delim="\t",col_names = T)

Loop_NPC_D868D_inNeuD868D_5kb <- read_delim("SETBP1_epigenomics/HiC/Neu_D868D_vs_NPC_D868D/hiccupsdiff/file1/requested_list_5000.bedpe",
                                            delim="\t",col_names = T)

Loop_NPC_D868D_inNeuD868D_10kb <- read_delim("SETBP1_epigenomics/HiC/Neu_D868D_vs_NPC_D868D/hiccupsdiff/file1/requested_list_10000.bedpe",
                                             delim="\t",col_names = T)

Loop_NPC_D868D_merged <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/NPC_D868D/hiccups/merged_loops.bedpe",
                                    delim="\t",col_names = T)


#join different resolution:


Loop_NPC_D868D_merged <- Loop_NPC_D868D_merged[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6)%>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2)


Loop_NPC_D868D_5kb <- Loop_NPC_D868D_5kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_NPC_D868D=12,
                expectedBL_NPC_D868D=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_NPC_D868D,
                expectedBL_NPC_D868D)

Loop_NPC_D868D_inNeuD868D_5kb <- Loop_NPC_D868D_inNeuD868D_5kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_Neu_D868D=12,
                expectedBL_Neu_D868D=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_Neu_D868D,
                expectedBL_Neu_D868D)

Loop_Neu_D868D_vs_NPC_D868D_5kb <- inner_join(Loop_NPC_D868D_5kb,Loop_NPC_D868D_inNeuD868D_5kb)


Loop_NPC_D868D_10kb <- Loop_NPC_D868D_10kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_NPC_D868D=12,
                expectedBL_NPC_D868D=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_NPC_D868D,
                expectedBL_NPC_D868D)

Loop_NPC_D868D_inNeuD868D_10kb <- Loop_NPC_D868D_inNeuD868D_10kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_Neu_D868D=12,
                expectedBL_Neu_D868D=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_Neu_D868D,
                expectedBL_Neu_D868D)

Loop_Neu_D868D_vs_NPC_D868D_10kb <- inner_join(Loop_NPC_D868D_10kb,Loop_NPC_D868D_inNeuD868D_10kb)


Loop_Neu_D868D_vs_NPCD868D_5kb_merged <- inner_join(Loop_NPC_D868D_merged,Loop_Neu_D868D_vs_NPC_D868D_5kb)

Loop_Neu_D868D_vs_NPCD868D_10kb_merged <- inner_join(Loop_NPC_D868D_merged,Loop_Neu_D868D_vs_NPC_D868D_10kb)

#Calculate a list with all loops asssociated to NPCs or Neu D868D

Loop_Neu_D868D_vs_NPC_D868D_merged_NPC_Loops <- rbind.data.frame(Loop_Neu_D868D_vs_NPCD868D_5kb_merged,Loop_Neu_D868D_vs_NPCD868D_10kb_merged) %>% 
  dplyr::mutate(intensity_NPC_D868D=observed_NPC_D868D/expectedBL_NPC_D868D) %>% 
  dplyr::mutate(intensity_Neu_D868D=observed_Neu_D868D/expectedBL_Neu_D868D) %>% 
  dplyr::mutate(fold=foldchange(intensity_Neu_D868D,intensity_NPC_D868D))

Loop_NPC_D868D_vs_Neu_D868D_all_loops <-  read_delim("/media/zaghi/SETBP1_Epigenomics/Hi-C/Loop_NPC_D868D_vs_Neu_D868D_all_loops",
                                                      delim="\t",col_names = T)

Loop_NPC_D868D_vs_D868D_all_loops_fc_minus2 <- Loop_NPC_D868D_vs_Neu_D868D_all_loops %>% 
  dplyr::filter(fold<=-2)

Loop_NPC_D868D_vs_D868D_all_loops_fc_minus1.5 <- Loop_NPC_D868D_vs_Neu_D868D_all_loops %>% 
  dplyr::filter(fold<=-1.5) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_Neu_D868D_all_loops_fc_minus1.5.bedpe")

Loop_NPC_D868D_vs_D868D_all_loops__fc_plus2 <- Loop_NPC_D868D_vs_Neu_D868D_all_loops %>% 
  dplyr::filter(fold>=2)

Loop_NPC_D868D_vs_D868D_all_loops_fc_plus1.5 <- Loop_NPC_D868D_vs_Neu_D868D_all_loops  %>% 
  dplyr::filter(fold>=1.5)%>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_Neu_D868D_all_loops_fc_plus1.5.bedpe")

Loop_NPC_D868D_vs_Neu_D868D_all_loops_unchanged <-  Loop_NPC_D868D_vs_Neu_D868D_all_loops  %>% 
  dplyr::filter(fold<=1.5) %>% 
  dplyr::filter(fold>=-1.5) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_Neu_D868D_all_loops_unchanged.bedpe")


#Loop strength calculation NPCs D868N vs Neu D868N (Neu D868N loops lists)

Loop_Neu_D868N_5kb <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/Neu_D868N/hiccups/enriched_pixels_5000.bedpe",
                                 delim="\t",col_names = T)

Loop_Neu_D868N_10kb <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/Neu_D868N/hiccups/enriched_pixels_10000.bedpe",
                                  delim="\t",col_names = T)


Loop_NPC_D868N_inNeuD868N_5kb <- read_delim("SETBP1_epigenomics/HiC/Neu_D868N_vs_NPC_D868N/hiccupsdiff/file2/requested_list_5000.bedpe",
                                            delim="\t",col_names = T)

Loop_NPC_D868N_inNeuD868N_10kb <- read_delim("SETBP1_epigenomics/HiC/Neu_D868N_vs_NPC_D868N/hiccupsdiff/file2/requested_list_10000.bedpe",
                                             delim="\t",col_names = T)


Loop_Neu_D868N_merged <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/Neu_D868N/hiccups/merged_loops.bedpe",
                                    delim="\t",col_names = T)

#join different resolution:


Loop_Neu_D868N_merged <- Loop_Neu_D868N_merged[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6)%>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2)


Loop_Neu_D868N_5kb <- Loop_Neu_D868N_5kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_Neu_D868N=12,
                expectedBL_Neu_D868N=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_Neu_D868N,
                expectedBL_Neu_D868N)

Loop_NPC_D868N_inNeuD868N_5kb <- Loop_NPC_D868N_inNeuD868N_5kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_NPC_D868N=12,
                expectedBL_NPC_D868N=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_NPC_D868N,
                expectedBL_NPC_D868N)

Loop_Neu_D868N_vs_NPC_D868N <- inner_join(Loop_Neu_D868N_5kb,Loop_NPC_D868N_inNeuD868N_5kb)


Loop_Neu_D868N_10kb <- Loop_Neu_D868N_10kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_Neu_D868N=12,
                expectedBL_Neu_D868N=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_Neu_D868N,
                expectedBL_Neu_D868N)

Loop_NPC_D868N_inNeuD868N_10kb <- Loop_NPC_D868N_inNeuD868N_10kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_NPC_D868N=12,
                expectedBL_NPC_D868N=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_NPC_D868N,
                expectedBL_NPC_D868N)

Loop_Neu_D868N_vs_NPC_D868N_10kb <- inner_join(Loop_Neu_D868N_10kb,Loop_NPC_D868N_inNeuD868N_10kb)

Loop_Neu_D868N_vs_NPCD868N_5kb_merged <- inner_join(Loop_Neu_D868N_merged,Loop_Neu_D868N_vs_NPC_D868N)

Loop_Neu_D868N_vs_NPCD868N_10kb_merged <- inner_join(Loop_Neu_D868N_merged,Loop_Neu_D868N_vs_NPC_D868N_10kb)

Loop_Neu_D868N_vs_NPC_D868N_merged_Neu_Loops <- rbind.data.frame(Loop_Neu_D868N_vs_NPCD868N_5kb_merged,Loop_Neu_D868N_vs_NPCD868N_10kb_merged) %>% 
  dplyr::mutate(intensity_NPC_D868N=observed_NPC_D868N/expectedBL_NPC_D868N) %>% 
  dplyr::mutate(intensity_Neu_D868N=observed_Neu_D868N/expectedBL_Neu_D868N) %>% 
  dplyr::mutate(fold=foldchange(intensity_Neu_D868N,intensity_NPC_D868N))


#Loop strength calculation NPCs D868N vs Neu D868N (NPCs D868N loops lists)

Loop_NPC_D868N_5kb <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/NPC_D868N/hiccups/enriched_pixels_5000.bedpe",
                                 delim="\t",col_names = T)

Loop_NPC_D868N_10kb <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/NPC_D868N/hiccups/enriched_pixels_10000.bedpe",
                                  delim="\t",col_names = T)

Loop_NPC_D868N_inNeuD868N_5kb <- read_delim("SETBP1_epigenomics/HiC/Neu_D868N_vs_NPC_D868N/hiccupsdiff/file1/requested_list_5000.bedpe",
                                            delim="\t",col_names = T)

Loop_NPC_D868N_inNeuD868N_10kb <- read_delim("SETBP1_epigenomics/HiC/Neu_D868N_vs_NPC_D868N/hiccupsdiff/file1/requested_list_10000.bedpe",
                                             delim="\t",col_names = T)

Loop_NPC_D868N_merged <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/NPC_D868N/hiccups/merged_loops.bedpe",
                                    delim="\t",col_names = T)


#join different resolution:


Loop_NPC_D868N_merged <- Loop_NPC_D868N_merged[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6)%>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2)


Loop_NPC_D868N_5kb <- Loop_NPC_D868N_5kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_NPC_D868N=12,
                expectedBL_NPC_D868N=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_NPC_D868N,
                expectedBL_NPC_D868N)

Loop_NPC_D868N_inNeuD868N_5kb <- Loop_NPC_D868N_inNeuD868N_5kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_Neu_D868N=12,
                expectedBL_Neu_D868N=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_Neu_D868N,
                expectedBL_Neu_D868N)

Loop_Neu_D868N_vs_NPC_D868N_5kb <- inner_join(Loop_NPC_D868N_5kb,Loop_NPC_D868N_inNeuD868N_5kb)


Loop_NPC_D868N_10kb <- Loop_NPC_D868N_10kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_NPC_D868N=12,
                expectedBL_NPC_D868N=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_NPC_D868N,
                expectedBL_NPC_D868N)

Loop_NPC_D868N_inNeuD868N_10kb <- Loop_NPC_D868N_inNeuD868N_10kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_Neu_D868N=12,
                expectedBL_Neu_D868N=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_Neu_D868N,
                expectedBL_Neu_D868N)

Loop_Neu_D868N_vs_NPC_D868N_10kb <- inner_join(Loop_NPC_D868N_10kb,Loop_NPC_D868N_inNeuD868N_10kb)

Loop_Neu_D868N_vs_NPCD868N_5kb_merged <- inner_join(Loop_NPC_D868N_merged,Loop_Neu_D868N_vs_NPC_D868N_5kb)

Loop_Neu_D868N_vs_NPCD868N_10kb_merged <- inner_join(Loop_NPC_D868N_merged,Loop_Neu_D868N_vs_NPC_D868N_10kb)

Loop_Neu_D868N_vs_NPC_D868N_merged_NPC_Loops <- rbind.data.frame(Loop_Neu_D868N_vs_NPCD868N_5kb_merged,Loop_Neu_D868N_vs_NPCD868N_10kb_merged) %>% 
  dplyr::mutate(intensity_NPC_D868N=observed_NPC_D868N/expectedBL_NPC_D868N) %>% 
  dplyr::mutate(intensity_Neu_D868N=observed_Neu_D868N/expectedBL_Neu_D868N) %>% 
  dplyr::mutate(fold=foldchange(intensity_Neu_D868N,intensity_NPC_D868N))


#Calculate a list with all loops asssociated to NPCs or Neu D868N

Loop_NPC_D868N_vs_Neu_D868N_all_loops <- rbind.data.frame(Loop_Neu_D868N_vs_NPC_D868N_merged_NPC_Loops,Loop_Neu_D868N_vs_NPC_D868N_merged_Neu_Loops) %>% 
  unique()%>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868N_vs_Neu_D868N_all_loops.bedpe")

Loop_NPC_D868N_vs_D868N_all_loops_fc_minus2 <- Loop_NPC_D868N_vs_Neu_D868N_all_loops %>% 
  dplyr::filter(fold<=-2)

Loop_NPC_D868N_vs_D868N_all_loops_fc_minus1.5 <- Loop_NPC_D868N_vs_Neu_D868N_all_loops  %>% 
  dplyr::filter(fold<=-1.5)%>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868N_vs_Neu_D868N_all_loops_fc_minus1.5.bedpe")

Loop_NPC_D868N_vs_D868N_all_loops__fc_plus2 <- Loop_NPC_D868N_vs_Neu_D868N_all_loops %>% 
  dplyr::filter(fold>=2)

Loop_NPC_D868N_vs_D868N_all_loops_fc_plus1.5 <- Loop_NPC_D868N_vs_Neu_D868N_all_loops  %>% 
  dplyr::filter(fold>=1.5)%>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868N_vs_Neu_D868N_all_loops_fc_plus1.5.bedpe")

Loop_NPC_D868N_vs_Neu_D868N_all_loops_unchanged <-  Loop_NPC_D868N_vs_Neu_D868N_all_loops  %>% 
  dplyr::filter(fold<=1.5) %>% 
  dplyr::filter(fold>=-1.5) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868N_vs_Neu_D868N_all_loops_unchanged.bedpe")
