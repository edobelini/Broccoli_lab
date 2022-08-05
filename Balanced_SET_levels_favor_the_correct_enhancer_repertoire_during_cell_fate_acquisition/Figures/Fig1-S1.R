library(dplyr)
library(tidyverse)
library(ggpubr)
library(clipr)
library(GenomicRanges)
library(regioneR)
library(rtracklayer)
library(IRanges)
library(dplyr)
library(plyranges)
library(ChIPseeker)
library(org.Hs.eg.db)


#Load multiBigWigSummary containig all annotation of SGS ATAC peaks patients
Open_chromatin_all_ATAC_SGS <- read_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate", col_names = T) %>% 
  dplyr::rename(chr=1,
                start=2,
                end=3,
                NPC_D868D=4,
                NPC_D868N=5,
                NPC_I871I=6,
                NPC_I871T=7,
                Neu_D868D=8,
                Neu_D868N=9,
                dir_Neu_D868D=10,
                dir_Neu_D868N=11,
                NPC_D868D_peaks=12,
                NPC_D868N_peaks=13,
                NPC_I871I_peaks=14,
                NPC_I871T_peaks=15,
                Neu_D868D_peaks=16,
                Neu_D868N_peaks=17,
                dir_Neu_D868D_peaks=18,
                dir_Neu_D868N_peaks=19,
                Neu_D868D_up_dev=20,
                Neu_D868D_down_dev=21,
                Neu_D868N_up_dev=22,
                Neu_D868N_down_dev=23) %>% 
  dplyr::mutate( NPC_D868D_peaks=NPC_D868D_peaks/NPC_D868D_peaks,
                 NPC_D868N_peaks=NPC_D868N_peaks/NPC_D868N_peaks,
                 NPC_I871I_peaks=NPC_I871I_peaks/NPC_I871I_peaks,
                 NPC_I871T_peaks=NPC_I871T_peaks/NPC_I871T_peaks,
                 Neu_D868D_peaks=Neu_D868D_peaks/Neu_D868D_peaks,
                 Neu_D868N_peaks=Neu_D868N_peaks/Neu_D868N_peaks,
                 dir_Neu_D868D_peaks=dir_Neu_D868D_peaks/dir_Neu_D868D_peaks,
                 dir_Neu_D868N_peaks=dir_Neu_D868N_peaks/dir_Neu_D868N_peaks) %>% 
  mutate_each(funs(replace(., is.na(.), 0))) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate")



