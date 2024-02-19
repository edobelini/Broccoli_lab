library(tidyverse)
library(dplyr)
library(plyranges)


#NPC Loop strength:

Loop_NPC_D868D_5kb <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/NPC_D868D/hiccups/enriched_pixels_5000.bedpe",
                                 delim="\t",col_names = T)

Loop_NPC_D868D_10kb <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/NPC_D868D/hiccups/enriched_pixels_10000.bedpe",
                                 delim="\t",col_names = T)

Loop_NPC_D868D_inD868N_5kb <- read_delim("SETBP1_epigenomics/HiC/NPC_D868D_vs_NPC_D868N//hiccupsdiff/file2/requested_list_5000.bedpe",
                                         delim="\t",col_names = T)

Loop_NPC_D868D_inD868N_10kb <- read_delim("SETBP1_epigenomics/HiC/NPC_D868D_vs_NPC_D868N/hiccupsdiff/file2/requested_list_10000.bedpe",
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
                end2=6) %>% 
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
                observed_D868D=12,
                expectedBL_D868D=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_D868D,
                expectedBL_D868D)

Loop_NPC_D868D_inD868N_5kb <- Loop_NPC_D868D_inD868N_5kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_D868N=12,
                expectedBL_D868N=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_D868N,
                expectedBL_D868N)


Loop_NPC_D868D_vs_D868N_5kb <- inner_join(Loop_NPC_D868D_5kb,Loop_NPC_D868D_inD868N_5kb)


Loop_NPC_D868D_10kb <- Loop_NPC_D868D_10kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_D868D=12,
                expectedBL_D868D=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_D868D,
                expectedBL_D868D)

Loop_NPC_D868D_inD868N_10kb <- Loop_NPC_D868D_inD868N_10kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_D868N=12,
                expectedBL_D868N=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_D868N,
                expectedBL_D868N)

Loop_NPC_D868D_vs_D868N_10kb <- inner_join(Loop_NPC_D868D_10kb,Loop_NPC_D868D_inD868N_10kb)

Loop_NPC_D868D_vs_D868N_5kb_merged <- inner_join(Loop_NPC_D868D_merged,Loop_NPC_D868D_vs_D868N_5kb)

Loop_NPC_D868D_vs_D868N_10kb_merged <- inner_join(Loop_NPC_D868D_merged,Loop_NPC_D868D_vs_D868N_10kb)

Loop_NPC_D868D_vs_D868N_merged_Ctrl_loops <- rbind.data.frame(Loop_NPC_D868D_vs_D868N_5kb_merged,Loop_NPC_D868D_vs_D868N_10kb_merged) %>% 
  dplyr::mutate(intensity_D868D=observed_D868D/expectedBL_D868D) %>% 
  dplyr::mutate(intensity_D868N=observed_D868N/expectedBL_D868N) %>% 
  dplyr::mutate(fold=foldchange(intensity_D868N,intensity_D868D)) %>% 
  dplyr::mutate(chr2=chr1)

Loop_NPC_D868D_vs_D868N_fc_minus2 <- Loop_NPC_D868D_vs_D868N_merged_Ctrl_loops %>% 
  dplyr::filter(fold<=-2)%>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_D868N_fc_plus2.bedpe")

Loop_NPC_D868D_vs_D868N_fc_minus1.5 <- Loop_NPC_D868D_vs_D868N_merged_Ctrl_loops %>% 
  dplyr::filter(fold<=-1.5) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_D868N_fc_minus1.5.bedpe")


Loop_NPC_D868D_vs_D868N_fc_plus2 <- Loop_NPC_D868D_vs_D868N_merged_Ctrl_loops %>% 
  dplyr::filter(fold>=2)%>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_D868N_fc_plus2.bedpe")

Loop_NPC_D868D_vs_D868N_fc_plus1.5 <- Loop_NPC_D868D_vs_D868N_merged_Ctrl_loops %>% 
  dplyr::filter(fold>=1.5)%>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_D868N_fc_plus1.5.bedpe")


Loop_NPC_D868D_vs_D868N_fc_unchanged <- Loop_NPC_D868D_vs_D868N_merged_Ctrl_loops %>% 
  dplyr::filter(fold<=1.5)%>% 
  dplyr::filter(fold>=-1.5) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_D868N_fc_unchanged.bedpe")