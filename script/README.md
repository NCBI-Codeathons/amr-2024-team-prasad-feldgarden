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
sudo apt-get install -y libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
sudo apt-get install -y libfontconfig1-dev
sudo apt-get install libharfbuzz-dev libfribidi-dev
R
install.packages('tidyverse')
install.packages('XML')
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages('devtools')
devtools:::install_github("Bioconductor/GenomeInfoDb")
install.packages('optparse')

```

### Scanner.R script

This script is called as follows:

```bash
DGW.R -i query.fa -d ref.fa -o output_folder
```
Options are:
```bash
-i input file
-d reference file
-o output folder
-v verbose
-p prefix for output files, default is '"YYYYMMDD_HHMMSS"
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
micromamba install -y conda-forge::biopython bioconda::blast bioconda::diamond
```

Starts with a blastx scan of the query sequences against the reference protein database.

```js
blastx -query query.fna -db reference_db -out results.xml -outfmt 5
```

```js
## for frameshift detection, diamond (v2.1.8) blastx has option 
## --frameshift to allow frameshift gap aligned and 
## show '/' or '\' in the alignement result

diamond blastx -q query.fna -d database/references.faa -f 5 --min-orf 1 --frameshift 15  -o results.xml
```

Run the python script to check nonsense mutation 

```js
DGW.py -i test/input.xml 
```

Including blast step.
```js
## still need pre-build the database either from blast or diamond.

 ./DGW_blast.py -i query.fna -d database/references.faa --diamond -o -
```

```js
usage: DGW_blast.py [-h] -i INPUT_FASTA -d DATABASE [-o OUTPUT] [--cov COV] [--id ID] [--indels] [--diamond] [--threads THREADS] [--verbose] [--version]

    ____                 _    ____                  __        __    _ _    _             
   |  _ \  ___  __ _  __| |  / ___| ___ _ __   ___  \ \      / __ _| | | _(_)_ __   __ _ 
   | | | |/ _ \/ _` |/ _` | | |  _ / _ | '_ \ / _ \  \ \ /\ / / _` | | |/ | | '_ \ / _` |
   | |_| |  __| (_| | (_| | | |_| |  __| | | |  __/   \ V  V | (_| | |   <| | | | | (_| |
   |____/ \___|\__,_|\__,_|  \____|\___|_| |_|\___|    \_/\_/ \__,_|_|_|\_|_|_| |_|\__, |
                                                                                   |___/                                                                                                                                                                                                                                                                                                                                                                                          
  Detect nonsensus/frameshift mutations by running blastx on contig fasta against target protein database

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTA, --input_fasta INPUT_FASTA
                        contig fasta file
  -d DATABASE, --database DATABASE
                        protein (diamond) blast DB
  -o OUTPUT, --output OUTPUT
                        output tsv file
  --cov COV             target coverage
  --id ID               hit identity
  --indels              check in-frame insertion and deletion
  --diamond             use diamond blastx
  --threads THREADS     threads
  --verbose             Show more information in log
  --version             show program's version number and exit
```

