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

#Basic AA colors - Shapely
aa_colors <- c(
  "M" = "#E6E600", "A" = "#C8C8C8", "L" = "#0F820F", "V" = "#0F820F", "I" = "#0F820F",
  "F" = "#3232AA", "W" = "#B45AB4", "Y" = "#3232AA", "T" = "#FA9600", "S" = "#FA9600",
  "R" = "#145AFF", "N" = "#00DCDC", "D" = "#E60A0A", "C" = "#E6E600", "E" = "#E60A0A",
  "Q" = "#00DCDC", "G" = "#EBEBEB", "H" = "#8282D2", "K" = "#145AFF", "P" = "#3232AA",
  "Match" = "dodgerblue",  
  "Mismatch" = "darkorchid",  
  "STOP" = "black",  
  "-" = "white" ,
  "*"="black"
)

## Functions 

ConscMatch <- function(ref_seq, tgt_seq) {
  ref_split <- str_split(as.character(ref_seq), "", simplify = TRUE)
  tgt_split <- str_split(as.character(tgt_seq), "", simplify = TRUE)
  
  seq_length <- min(ncol(ref_split), ncol(tgt_split))
  if(seq_length<max(ncol(ref_split), ncol(tgt_split))){
    tgt_split = c(tgt_split,rep("-",ncol(ref_split)-ncol(tgt_split)))
  }
  
  comparison_tbl <- tibble(
    Pos = 1:length(ref_split),  
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

#Looking for arguments
arg_list <- list(
  make_option(c("-i", "--in"), type = "character", dest="input", default = NULL,
              help = "input file path", metavar = "FILE"),
  make_option(c("-d", "--db"), type = "character", default = NULL,dest="db",
              help = "input reference file path", metavar = "FILE"),
  make_option(c("-o", "--out"), type = "character", default = NULL,dest="out",
              help = "output file path", metavar = "FILE"),
  make_option(c("-v", "--verbose"), action="store_true",default = FALSE, 
              dest="verbose", help="Print debug text"),
  make_option(c("-p", "--prefix"), default = format(Sys.time(), "%Y%m%d_%H%M%S"), type="character",metavar = "STRING",
              dest="pref", help="Prefix name for outputfile")
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
                 "-outfmt", "5"),  # Output format 6: tabular
        stdout = TRUE)
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

# Initialize a data frame to store the results
df <- data.frame(Query = character(),
                 Hit_ID = character(),
                 Hit_Description = character(),
                 E_value = numeric(),
                 Score = numeric(),
                 start=numeric(),
                 stop=numeric(),
                 Hit_Sequence = character(),
                 Amino_Acid_Sequence = character(),
                 ntid=numeric(),
                 aaid=numeric(),
                 alilen=numeric(),
                 lesion_type=character(),
                 stringsAsFactors = FALSE)

for (query in queries) {
  query_def <- xmlValue(getNodeSet(query, "Iteration_query-def")[[1]])
  
  hits <- getNodeSet(query, ".//Hit")
  
  for (hit in hits) {
    hit_id <- xmlValue(getNodeSet(hit, "Hit_id")[[1]])
    hit_def <- xmlValue(getNodeSet(hit, "Hit_def")[[1]])
    hsps <- getNodeSet(hit, ".//Hsp")
    hlen <-  xmlValue(getNodeSet(hit, "Hit_len")[[1]])
    
    for (hsp in hsps) {
      e_value <- as.numeric(xmlValue(getNodeSet(hsp, "Hsp_evalue")[[1]]))
      score <- as.numeric(xmlValue(getNodeSet(hsp, "Hsp_score")[[1]]))
      hit_sequence <- xmlValue(getNodeSet(hsp, "Hsp_qseq")[[1]])
      hfrom <- xmlValue(getNodeSet(hsp,"Hsp_hit-from")[[1]])
      hto <- xmlValue(getNodeSet(hsp,"Hsp_hit-to")[[1]])
      qfrom <- xmlValue(getNodeSet(hsp,"Hsp_query-from")[[1]])
      qto <- xmlValue(getNodeSet(hsp,"Hsp_query-to")[[1]])
      alilen <- xmlValue(getNodeSet(hsp,"Hsp_align-len")[[1]])
      
      #Sequence extraction and translation
      dna_seq <- DNAString(gsub("-","",hit_sequence))
      
      #Alignment
      if((as.numeric(hto)-as.numeric(hfrom))<0){
        alref=translate(reverseComplement(oref[hit_def]), genetic.code = getGeneticCode("11"))
      } else {
        alref=translate(oref[hit_def], genetic.code = getGeneticCode("11"))
      }
      if((as.numeric(qto)-as.numeric(qfrom))<0){
        altgt <- translate(reverseComplement(dna_seq), genetic.code = getGeneticCode("11"))
      } else {
        altgt <- translate(dna_seq, genetic.code = getGeneticCode("11"))
      }
      altgt <-AAString(gsub("[*].*$","*",as.character(altgt)))
      
      local_align <- pwalign::pairwiseAlignment(altgt, alref, type = "global-local", substitutionMatrix = "BLOSUM62", gapOpening = 20, gapExtension = 0.5)
      
      aamatch=nmatch(local_align)
      aalength=nchar(alref)
      aa_id=aamatch/aalength*100
      
      #Id calc
      hsp_identity <- as.numeric(xmlValue(getNodeSet(hsp, "Hsp_identity")[[1]]))
      hsp_align_len <- as.numeric(xmlValue(getNodeSet(hsp, "Hsp_align-len")[[1]]))
      ntid_per <- (hsp_identity / hsp_align_len) * 100
      
      CM=ConscMatch(subject(local_align),pattern(local_align))
      
      if(length(altgt)<length(alref[[1]]) & as.character(altgt[length(altgt)])=="*"){
        if(opt$verbose==TRUE){
          cat("Maybe STOP mutation\n")
        }
        ltype=CM[[2]] %>% filter(Pos >= (max(Pos) - 20)) %>% 
          group_by(Altype) %>%
          summarise(n=n()) %>% 
          spread(Altype,n)
        if(!"Mismatch" %in% colnames(ltype)){
          ltype=cbind(ltype,data.frame("Mismatch"=0))
        }
        ltype=ltype %>% 
          mutate(lesion_type=ifelse(Mismatch>Match,"STOP frame shift","STOP point")) %>% 
          select(lesion_type) %>% 
          pull()
      }else{
        ltype="Intact?"
      }
      
      df <- df %>% add_row(Query = query_def,
                           Hit_ID = hit_id,
                           Hit_Description = hit_def,
                           start=as.numeric(qfrom),
                           stop=as.numeric(qto),
                           E_value = e_value,
                           Score = score,
                           Hit_Sequence = hit_sequence,
                           Amino_Acid_Sequence = as.character(altgt),
                           ntid=ntid_per,
                           aaid=aa_id,
                           alilen=as.numeric(alilen)/as.numeric(hlen)*100,
                           lesion_type=ltype)
    }
  }
}

#filtering - On hold need to check for overlapping position
#df <- df %>% 
# group_by(Query) %>% 
# filter(ntid==max(ntid)) %>% 
# unique()

df2=df %>%
  separate(Hit_Description,into=c("Gene"), sep="_",extra="drop") %>% 
  select(Gene,ntid,aaid,alilen,lesion_type)

cat(" ")
cat(query_def)
print(df2)

write_tsv(df,file.path(finalout,paste0(opt$pref,"_results.tsv")))
