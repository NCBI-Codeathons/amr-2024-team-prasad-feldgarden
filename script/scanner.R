library(tidyverse)
library(XML)
library(Biostrings)

xml_data <- xmlParse("data/result.xml")

ref=readDNAStringSet("data/ompk36/ref.fna")
aaref=translate(ref, genetic.code = getGeneticCode("11"))
ref=as.data.frame(ref)


meta=read_tsv("data/ompk36/meta.tsv", name_repair = "universal")

# Extract queries (Iterations)
queries <- getNodeSet(xml_data, "//Iteration")

# Initialize a data frame to store the results
df <- data.frame(Query = character(),
                        Hit_ID = character(),
                        Hit_Description = character(),
                        E_value = numeric(),
                        Score = numeric(),
                        Hit_Sequence = character(),
                        Amino_Acid_Sequence = character(),
                        ntid=numeric(),
                        aaid=numeric(),
                        stringsAsFactors = FALSE)
for (query in queries) {
  query_def <- xmlValue(getNodeSet(query, "Iteration_query-def")[[1]])
  
  hits <- getNodeSet(query, ".//Hit")
  
  for (hit in hits) {
    hit_id <- xmlValue(getNodeSet(hit, "Hit_id")[[1]])
    hit_def <- xmlValue(getNodeSet(hit, "Hit_def")[[1]])
    hsps <- getNodeSet(hit, ".//Hsp")
    
    for (hsp in hsps) {
      e_value <- as.numeric(xmlValue(getNodeSet(hsp, "Hsp_evalue")[[1]]))
      score <- as.numeric(xmlValue(getNodeSet(hsp, "Hsp_score")[[1]]))
      hit_sequence <- xmlValue(getNodeSet(hsp, "Hsp_hseq")[[1]])
      
      #Sequence extraction and translation
      dna_seq <- DNAString(gsub("-","",hit_sequence))
      amino_acid_seq <- translate(dna_seq, genetic.code = getGeneticCode("11"))
      
      #Alignment
      local_align <- pwalign::pairwiseAlignment(amino_acid_seq, aaref, type = "local", substitutionMatrix = "BLOSUM62", gapOpening = 10, gapExtension = 1)
      aamatch=nmatch(local_align)
      aalength=nchar(aaref)
      aa_id=aamatch/aalength*100
      
      #Id calc
      hsp_identity <- as.numeric(xmlValue(getNodeSet(hsp, "Hsp_identity")[[1]]))
      hsp_align_len <- as.numeric(xmlValue(getNodeSet(hsp, "Hsp_align-len")[[1]]))
      ntid_per <- (hsp_identity / hsp_align_len) * 100
      
      
      df <- df %>% add_row(Query = query_def,
                                         Hit_ID = hit_id,
                                         Hit_Description = hit_def,
                                         E_value = e_value,
                                         Score = score,
                                         Hit_Sequence = hit_sequence,
                                         Amino_Acid_Sequence = as.character(amino_acid_seq),
                                         ntid=ntid_per,
                                         aaid=aa_id)
    }
  }
}

return(df)
}