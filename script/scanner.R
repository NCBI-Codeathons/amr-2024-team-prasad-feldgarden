#!/usr/bin/env Rscript

## ---------------------------
##
## Script name: Scanner.R (to be adjusted)
##
## Project: "Codeathon to combat antimicrobial resistance at NCBI"
##
## Author: Dr. Adrien Assié
## Date Created: 2024-09-23
##
## ---------------------------

#Libraries

suppressMessages(library(tidyverse))
library(XML)
suppressMessages(library(Biostrings))
library(optparse)
suppressMessages(library(parallel))
suppressMessages(library(doParallel))

## Functions 

ConscMatch <- function(ref_seq, tgt_seq) {
  t1=gsub("-+\\*?$","*",as.character(tgt_seq))
  ref_split <- str_split(as.character(ref_seq), "", simplify = TRUE)
  tgt_split <- str_split(as.character(t1), "", simplify = TRUE)
  
  seq_length <- min(ncol(ref_split), ncol(tgt_split))
  if(seq_length<max(ncol(ref_split), ncol(tgt_split))){
    ref_split = ref_split[1:length(tgt_split)]
  }
  
  comparison_tbl <- tibble(
    Pos = 1:length(tgt_split),  
    ref_aa = as.character(ref_split),         
    tgt_aa = as.character(tgt_split)
  ) %>%
    mutate(
      Altype = case_when(           
        tgt_aa == "*" ~ "STOP",
        ref_aa == tgt_aa ~ "Match",          
        tgt_aa == "-" ~ "gap",                          
        tgt_aa != ref_aa ~ "Mismatch"
      ),
      cumulative_match = cumsum(Altype != "Mismatch" & Altype != "")
    ) 
  
  # Count consecutive matches
  consecutive_matches <- comparison_tbl %>%
    filter(Altype == "Match") %>%
    nrow()
  
  return(list(consecutive_matches, comparison_tbl))
}

detect_frameshift_positions <- function(ref, target, verbose) {
  tryCatch(
    {
  ref_chars <- as.character(str_split(as.character(ref), "", simplify = TRUE))
  target_chars <- as.character(str_split(as.character(target), "", simplify = TRUE))
  
  frame_offset <- 0  # tracks cumulative changes in frame
  
  i=1
  a="A"
  while (i < length(ref_chars) & a=="A") {
    if (target_chars[i] == "-") {
      #print("FOUND")
      frame_offset <- frame_offset + 1
      a="tplus"
      frame=i
      i=i+1
      } else if (ref_chars[i] == "-") {
      #print("FOUND")
      frame_offset <- frame_offset + 1
      a="rplus"
      frame=i
      i=i+1
      } else {
      i=i+1
    }
  }
  
  if(a=="tplus"){
    while(i < length(ref_chars) & a=="tplus") {
      if (target_chars[i] == "-") {
        # A gap in the target sequence
        frame_offset <- frame_offset + 1
        i=i+1
        } else {
          a="C"
        }
    }
  } else if(a=="rplus"){
    while(i < length(ref_chars) & a=="rplus") {
      if (ref_chars[i] == "-") {
        # A gap in the target sequence
        frame_offset <- frame_offset + 1
        i=i+1
      } else {
        a="C"
      }
    }
      
    }
  if (i == length(ref_chars)){
    ltype="sbs"
  } else if (frame_offset %% 3 != 0) {
    ltype=paste0(floor((frame-1)/3+1),"fs")
  } else {
    ltype="nfs"
  }
  return(ltype)
    },
  error = function(e) {
    # What to do if an error occurs
    if(verbose==TRUE){
      message("Error in frameshift encountered")
    }
    return("ERROR")
  }
  )
}

#Looking for arguments
arg_list <- list(
  make_option(c("-i", "--in"), type = "character", dest="input", default = NULL,
              help = "input file path - REQUIRED", metavar = "FILE"),
  make_option(c("-d", "--db"), type = "character", default = NULL,dest="db",
              help = "input reference file path - REQUIRED", metavar = "FILE"),
  make_option(c("-o", "--out"), type = "character", default = NULL,dest="out",
              help = "output file path - REQUIRED", metavar = "FILE"),
  make_option(c("-v", "--verbose"), action="store_true",default = FALSE, 
              dest="verbose", help="Print debug text - Default:FALSE"),
  make_option(c("-p", "--prefix"), default = format(Sys.time(), "%Y%m%d_%H%M%S"), type="character",metavar = "STRING",
              dest="pref", help="Prefix name for outputfile - Default:YYYYMMDD_HHMMSS"),
  make_option(c("-c", "--cpu"), default = 1, type="numeric",metavar = "NUMERIC",
              dest="core", help="Number of cpus to use - Default:1"),
  make_option(c("-m", "--more"), action="store_true", default = FALSE,
              dest="more", help="Create an additional table with sequence information and more - Default:FALSE")
)

# Create a parser object and parse arguments
opt_parser <- OptionParser(option_list = arg_list)
opt <- parse_args(opt_parser)

#Sanity check
if (is.null(opt$input) || is.null(opt$out) || is.null(opt$db)) {
  print_help(opt_parser)
  stop("All -i, -d and -o arguments must be supplied.", call. = FALSE)
}

query_fa <- opt$input
ref_fa <- opt$db
finalout <- opt$out

tmpDIR=tempdir()
output_db <- file.path(tempdir(), "blast_tmp")

#Preparing files
if(opt$verbose==TRUE){
  cat("Creating BLAST database \n")
}
system2(command = "makeblastdb",
        args = c("-dbtype nucl", 
                 "-in", ref_fa,
                 "-out",output_db,
                 "-title blasttmp"
        ),
        stdout = FALSE)

##Debug
if(opt$verbose==TRUE){
  if (file.exists(paste0(output_db, ".ndb"))) {
    cat("BLAST database successfully created at:", output_db, "\n")
  } else {
    cat("Failed to create BLAST database.\n")
  }
  
  cat("BLASTING \n")
}

output_res <- file.path(tmpDIR, "blast_res.xml")

system2(command = "blastn",
        args = c("-query", query_fa, 
                 "-db", output_db, 
                 "-out", output_res,
                 "-num_threads", opt$core,
                 "-outfmt", "5"),  # Output format 6: tabular
        stdout = TRUE)

#Code to clean xml when there is multithreading

if(opt$core>1){
  # Remove specific patterns
  system2(command = "gsed",
          args = c("-i", "'/BlastOutput/d;/xml/d;/DOCTYPE/d;/Paramet/d'", output_res),
          stdout = TRUE)
  
  #Rebuild xml
  
  system2(command = "gsed",
          args = c("-i", "'1i<BlastOutput_iterations>'", output_res),
          stdout = TRUE)
  
  system2(command = "gsed",
          args = c("-i", "'1i<BlastOutput>'", output_res),
          stdout = TRUE)
  
  system2(command = "gsed",
          args = c("-i", "'1i<!DOCTYPE BlastOutput>'", output_res),
          stdout = TRUE)
  
  system2(command = "gsed",
          args = c("-i", "'1i<?xml version=\"1.0\"?>'", output_res),
          stdout = TRUE)
  
  system2(command = "echo",
          args = c("'</BlastOutput_iterations>'", ">>", output_res),
          stdout = TRUE)
  
  system2(command = "echo",
          args = c("'</BlastOutput>'", ">>", output_res),
          stdout = TRUE)
}

##Debug
if(opt$verbose==TRUE){
  if (file.exists(paste0(output_res))) {
    cat("BLAST results successfully created at:", output_res, "\n")
  } else {
    cat("Failed to find BLAST results\n")
  }
  cat("Processing Results\n")
}

#Loading blast results
xml_data <- xmlParse(output_res)

#Getting sequences into the right format
ref=readDNAStringSet(ref_fa)
oref=ref
ref=as.data.frame(ref)

# Extract queries (Iterations)
queries <- getNodeSet(xml_data, "//Iteration")

# Detect the number of available cores
num_cores <- opt$core

# Parallel 
results <- mclapply(queries, function(query) {
  query_def <- xmlValue(getNodeSet(query, "Iteration_query-def")[[1]])
  hits <- getNodeSet(query, ".//Hit")
  
  query_df <- data.frame(Query = character(),
                         Hit_ID = character(),
                         Hit_Description = character(),
                         E_value = numeric(),
                         Score = numeric(),
                         start = numeric(),
                         stop = numeric(),
                         Hit_Sequence = character(),
                         Amino_Acid_Sequence = character(),
                         ntid = numeric(),
                         aaid = numeric(),
                         alilen = numeric(),
                         lesion_type = character(),
                         orientation = numeric(),
                         ref_start = numeric(),
                         ref_stop = numeric(),
                         stringsAsFactors = FALSE)
  
  hits <- getNodeSet(query, ".//Hit")
  for (hit in hits) {
    hit_id <- xmlValue(getNodeSet(hit, "Hit_id")[[1]])
    hit_def <- xmlValue(getNodeSet(hit, "Hit_def")[[1]])
    hsps <- getNodeSet(hit, ".//Hsp")
    hlen <- xmlValue(getNodeSet(hit, "Hit_len")[[1]])
    
    for (hsp in hsps) {
      e_value <- as.numeric(xmlValue(getNodeSet(hsp, "Hsp_evalue")[[1]]))
      score <- as.numeric(xmlValue(getNodeSet(hsp, "Hsp_score")[[1]]))
      hit_sequence <- xmlValue(getNodeSet(hsp, "Hsp_qseq")[[1]])
      hfrom <- xmlValue(getNodeSet(hsp, "Hsp_hit-from")[[1]])
      hto <- xmlValue(getNodeSet(hsp, "Hsp_hit-to")[[1]])
      qfrom <- xmlValue(getNodeSet(hsp, "Hsp_query-from")[[1]])
      qto <- xmlValue(getNodeSet(hsp, "Hsp_query-to")[[1]])
      alilen <- xmlValue(getNodeSet(hsp, "Hsp_align-len")[[1]])
      
      if (hto != 1 & hfrom != 1) {
        next
      }
      
      # Sequence extraction and translation
      dna_seq <- DNAString(gsub("-", "", hit_sequence))
      
      if ((as.numeric(hto) - as.numeric(hfrom)) < 0) {
        dna_seq=reverseComplement(dna_seq)
        alref <- translate((oref[hit_def]), genetic.code = getGeneticCode("11"))
        altgt <- translate(dna_seq, genetic.code = getGeneticCode("11"))
        orientation <- 2
      } else if ((as.numeric(qto) - as.numeric(qfrom)) < 0) {
        dna_seq=reverseComplement(dna_seq)
        alref <- translate(oref[hit_def], genetic.code = getGeneticCode("11"))
        altgt <- translate(dna_seq, genetic.code = getGeneticCode("11"))
        orientation <- 2
      } else {
        alref <- translate(oref[hit_def], genetic.code = getGeneticCode("11"))
        altgt <- translate((dna_seq), genetic.code = getGeneticCode("11"))
        orientation <- 2
      }
      
      altgt <- AAString(gsub("[*].*$", "*", as.character(altgt)))
      
      local_align <- pwalign::pairwiseAlignment(altgt, alref, type = "global-local", substitutionMatrix = "BLOSUM62", gapOpening = 20, gapExtension = 0.5)
      
      aamatch <- nmatch(local_align)
      aalength <- nchar(alref)
      aa_id <- aamatch / aalength * 100
      
      # Calculate percentage of reference covered
      ref_length <- width(subject(local_align))
      align_start <- start(subject(local_align))
      align_end <- end(subject(local_align))
      aligned_length <- align_end - align_start + 1
      percentage_covered <- (aligned_length / ref_length) * 100
      
      # Identity calculation
      hsp_identity <- as.numeric(xmlValue(getNodeSet(hsp, "Hsp_identity")[[1]]))
      hsp_align_len <- as.numeric(xmlValue(getNodeSet(hsp, "Hsp_align-len")[[1]]))
      ntid_per <- (hsp_identity / hsp_align_len) * 100
      
      CM <- ConscMatch(subject(local_align), pattern(local_align))
      
      if (length(altgt) < length(alref[[1]]) & as.character(altgt[length(altgt)]) == "*") {
        if (opt$verbose == TRUE) {
          cat("Maybe STOP mutation\n")
        }
        ltype <- CM[[2]] %>% filter(Pos >= (max(Pos) - 5)) %>% 
          group_by(Altype) %>%
          summarise(n = n()) %>% 
          spread(Altype, n)
        
        if (!"Mismatch" %in% colnames(ltype)) {
          if (pull(CM[[2]][nrow(CM[[2]]), 2]) == pull(CM[[2]][nrow(CM[[2]]), 3]) & pull(CM[[2]][nrow(CM[[2]]), 4]) == "STOP") {
            ltype <- "Ok?"
          } else if (pull(CM[[2]][nrow(CM[[2]]), 2]) != "*" & pull(CM[[2]][nrow(CM[[2]]), 3]) =="*") {
            ltype <- paste0(pull(CM[[2]][nrow(CM[[2]]), 2]),pull(CM[[2]][nrow(CM[[2]]), 1]),"STOP")
          } else {
            ltype <- paste0(pull((CM[[2]] %>% filter(Altype == "STOP") %>% select(Pos, ref_aa) %>% unite("combined", ref_aa, Pos, sep = ""))[1, ]), "STOP")
          }
        } else {
          ntalm=pwalign::pairwiseAlignment(oref[hit_def],dna_seq)
          ltmp <- detect_frameshift_positions(pattern(ntalm), subject(ntalm), opt$verbose)
          if (ltmp == "sbs") {
            ltype <-"sbs"
          } else if (!ltmp %in% c("STOP", "nfs")) {
            ltmp2 <- as.numeric(gsub("fs", "", ltmp))
            AAp <- alref[[1]][CM[[2]] %>% filter(Pos == floor(ltmp2-1/3+1)) %>% select(Pos) %>% pull]
            ltype <- paste0(AAp, ltmp)
          } else {
            ltype <- "Error"
          }
        }
      } else {
        ltype <- "Intact?"
      }
      
      # Add row to the temporary data frame
      query_df <- query_df %>% add_row(Query = query_def,
                                       Hit_ID = hit_id,
                                       Hit_Description = hit_def,
                                       start = as.numeric(qfrom),
                                       stop = as.numeric(qto),
                                       E_value = e_value,
                                       Score = score,
                                       Hit_Sequence = hit_sequence,
                                       Amino_Acid_Sequence = as.character(altgt),
                                       ntid = ntid_per,
                                       aaid = aa_id,
                                       alilen = percentage_covered,
                                       lesion_type = ltype,
                                       orientation = orientation,
                                       ref_start = as.numeric(hfrom),
                                       ref_stop = as.numeric(hto))
    }
  }
  
  return(query_df)
}, mc.cores = num_cores)

df <- bind_rows(results)

df2=df %>%
  separate(Hit_Description,into=c("Gene"), sep="_",extra="drop") %>% 
  select(Gene,ntid,aaid,alilen,lesion_type)

cat("Done\n")

print(df2)

fdf=df %>% 
  filter(!lesion_type %in% c("Intact?","sbs","Error","Ok?")) %>% 
  separate(Hit_Description,into=c("ref_accession","ref_desc"), sep=" ",extra="merge")%>%
  separate(ref_accession,into=c("element_symbol","ref_accession"), sep="_",extra="merge") %>%
  mutate(ref_accession=gsub(":.*","",ref_accession)) %>% 
  dplyr::rename(aa_identity=aaid,
                nt_identity=ntid,
                coverage=alilen,
                contig_start=start,
                contig_stop=stop,
                lesion=lesion_type) %>% 
  mutate(
    contig_acc=gsub(" .*","",Query),
    lesion_type = case_when(           
      str_detect(lesion, "fs") ~ "FrameShift",
      str_detect(lesion, "STOP") ~ "PointMutation"
    )) %>% 
  select(element_symbol,
         contig_acc,
         contig_start,
         contig_stop, 
         orientation,
         lesion,
         lesion_type,
         aa_identity,
         nt_identity,
         ref_accession,
         ref_desc,ref_start,
         ref_stop,coverage)

write_tsv(fdf,file.path(finalout,paste0(opt$pref,"_results.tsv")))

if(opt$more==TRUE){
  cat("Creating extra Table \n")
  write_tsv(df,file.path(finalout,paste0(opt$pref,"_extended_results.tsv")))
}

