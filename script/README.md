# Script folder

## Require blast and R

```
sudo apt-get install r-base
sudo apt-get install ncbi-blast+
```

## Blastn and translation comparison approach


```
sudo apt-get install libxml2 libxml2-dev
sudo apt-get install libcurl4-openssl-dev
R
install.packages('tidyverse')
install.packages('XML')
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages('devtools')
devtools:::install_github("Bioconductor/GenomeInfoDb")

```

### Scanner.R script

This script is called as follows:

```bash
Rscript DGW.R query.fa ref.fa output_folder
```
The steps are:
- Create nucleotide blast database
- blast query to database and export as xml
- extract the nucleotide hit sequences
- calculate the nucleotide identity to the reference in %,
- translate nucleotide hit sequences to Amino acid sequences
- align amino acid sequences to the reference one
- calculate the amino acid identity to the reference in %,

## Blastx and nonsense mutation 

Requires Bio
```
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
micromamba create -n codeathon
eval "$(micromamba shell hook --shell bash)"
micromamba activate codeathon
micromamba install python=3.12
micromamba install -c conda-forge biopython
```

Starts with a blastx scan of the query sequences against the reference protein database.

```js
blastx -query query.fasta -db reference_db -out results.xml -outfmt 5
```

Run the python script to check nonsense mutation 

```js
check_nonsense_mutations.py -i test/input.xml 
```
