library(tidyverse)


# The goal of this script is to combine all tsv files obtained from script.R into a single data frame. 


# Combine all *filelist.tsv files
##Adapted from https://stackoverflow.com/questions/69341214/merging-thousands-of-csv-files-into-a-single-dataframe-in-r & https://dcl-prog.stanford.edu/purrr-parallel.html



## Pasens
filenames <- list.files("/Users/anaelizondo-ramos/Desktop/test/Pasens_script.R/", pattern = "*.tsv")

pasens_combined <- purrr::map_dfr(paste0("/Users/anaelizondo-ramos/Desktop/test/Pasens_script.R/", filenames), read_tsv, .id = "index", col_types = list(
  element_symbol = col_character(),
  contig_acc = col_character(),
  contig_start = col_number(),
  contig_stop = col_number(),
  orientation = col_character(),
  lesion = col_character(),
  lesion_type = col_character(),
  aa_identity = col_number(),
  nt_identity = col_number(),
  ref_accession = col_character(),
  ref_desc = col_number(),
  ref_start = col_number(),
  ref_stop = col_number(),
  coverage = col_character()
)) %>% 
  mutate(asm_acc = filenames[as.numeric(index)]) %>% 
  mutate(asm_acc = str_extract(asm_acc, regex('GCA_[0-9]*.[0-9]'))) %>% 
  select(element_symbol, 
         asm_acc, 
         contig_acc, 
         contig_start, 
         contig_stop,
         orientation,
         lesion,
         lesion_type,
         aa_identity,
         nt_identity,
         ref_accession,
         ref_desc,
         ref_start,
         ref_stop,
         coverage
         )

write_tsv(pasens_combined, "/Users/anaelizondo-ramos/Desktop/test/Pasens_combined_scriptR.tsv")


## Paresist
filenames2 <- list.files("/Users/anaelizondo-ramos/Desktop/test/Paresis_script.R/", pattern = "*.tsv")

paresis_combined <- purrr::map_dfr(paste0("/Users/anaelizondo-ramos/Desktop/test/Paresis_script.R/", filenames2), read_tsv, .id = "index", col_types = list(
  element_symbol = col_character(),
  contig_acc = col_character(),
  contig_start = col_number(),
  contig_stop = col_number(),
  orientation = col_character(),
  lesion = col_character(),
  lesion_type = col_character(),
  aa_identity = col_number(),
  nt_identity = col_number(),
  ref_accession = col_character(),
  ref_desc = col_number(),
  ref_start = col_number(),
  ref_stop = col_number(),
  coverage = col_character()
))%>% 
  mutate(asm_acc = filenames2[as.numeric(index)]) %>% 
  mutate(asm_acc = str_extract(asm_acc, regex('GCA_[0-9]*.[0-9]'))) %>% 
  select(element_symbol, 
         asm_acc, 
         contig_acc, 
         contig_start, 
         contig_stop,
         orientation,
         lesion,
         lesion_type,
         aa_identity,
         nt_identity,
         ref_accession,
         ref_desc,
         ref_start,
         ref_stop,
         coverage
  )
  
write_tsv(paresis_combined, "/Users/anaelizondo-ramos/Desktop/test/Paresis_combined_scriptR.tsv")


## Kpsens
filenames3 <- list.files("/Users/anaelizondo-ramos/Desktop/test/Kpsens_script.R/", pattern = "*.tsv")

kpsens_combined <- purrr::map_dfr(paste0("/Users/anaelizondo-ramos/Desktop/test/Kpsens_script.R/", filenames3), read_tsv, .id = "index", col_types = list(
  element_symbol = col_character(),
  contig_acc = col_character(),
  contig_start = col_number(),
  contig_stop = col_number(),
  orientation = col_character(),
  lesion = col_character(),
  lesion_type = col_character(),
  aa_identity = col_number(),
  nt_identity = col_number(),
  ref_accession = col_character(),
  ref_desc = col_number(),
  ref_start = col_number(),
  ref_stop = col_number(),
  coverage = col_character()
)) %>% 
  mutate(asm_acc = filenames3[as.numeric(index)]) %>% 
  mutate(asm_acc = str_extract(asm_acc, regex('GCA_[0-9]*.[0-9]'))) %>% 
  select(element_symbol, 
         asm_acc, 
         contig_acc, 
         contig_start, 
         contig_stop,
         orientation,
         lesion,
         lesion_type,
         aa_identity,
         nt_identity,
         ref_accession,
         ref_desc,
         ref_start,
         ref_stop,
         coverage
  )

write_tsv(kpsens_combined, "/Users/anaelizondo-ramos/Desktop/test/Kpsens_combined_scriptR.tsv")


## Paresist
filenames4 <- list.files("/Users/anaelizondo-ramos/Desktop/test/Kpresis_script.R/", pattern = "*.tsv")

kpresis_combined <- purrr::map_dfr(paste0("/Users/anaelizondo-ramos/Desktop/test/Kpresis_script.R/", filenames4), read_tsv, .id = "index", col_types = list(
  element_symbol = col_character(),
  contig_acc = col_character(),
  contig_start = col_number(),
  contig_stop = col_number(),
  orientation = col_character(),
  lesion = col_character(),
  lesion_type = col_character(),
  aa_identity = col_number(),
  nt_identity = col_number(),
  ref_accession = col_character(),
  ref_desc = col_number(),
  ref_start = col_number(),
  ref_stop = col_number(),
  coverage = col_character()
)) %>% 
  mutate(asm_acc = filenames4[as.numeric(index)]) %>% 
  mutate(asm_acc = str_extract(asm_acc, regex('GCA_[0-9]*.[0-9]'))) %>% 
  select(element_symbol, 
         asm_acc, 
         contig_acc, 
         contig_start, 
         contig_stop,
         orientation,
         lesion,
         lesion_type,
         aa_identity,
         nt_identity,
         ref_accession,
         ref_desc,
         ref_start,
         ref_stop,
         coverage
  )

write_tsv(kpresis_combined, "/Users/anaelizondo-ramos/Desktop/test/Kpresis_combined_scriptR.tsv")

