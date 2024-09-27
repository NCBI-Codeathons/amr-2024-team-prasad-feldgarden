# Script folder

## Require blast and R

```
sudo apt-get install r-base
sudo apt-get install ncbi-blast+
```

## BLASTn and translation comparison approach

### Requirements

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
install.packages('doParallel')
```

### Scanner.R script

This script is called as follows:

```bash
Rscript scanner.R -i query.fa -d ref.fa -o output_folder
```
Options are:
```bash
Options:
	-i FILE, --in=FILE
		input file path - REQUIRED

	-d FILE, --db=FILE
		input reference file path - REQUIRED

	-o FILE, --out=FILE
		output file path - REQUIRED

	-v, --verbose
		Print debug text - Default:FALSE

	-p STRING, --prefix=STRING
		Prefix name for outputfile - Default:YYYYMMDD_HHMMSS

	-c NUMERIC, --cpu=NUMERIC
		Number of cpus to use - Default:1

	-m, --more
		Create an additional table with sequence information and more - Default:FALSE

	-h, --help
		Show this help message and exit
```

The `scanner.R` script follows the following steps:
- Create nucleotide BLAST reference database
- BLAST query to database and export as XML
- Extract the nucleotide hit sequences
- Calculate the nucleotide identity to the reference in %
- Translate nucleotide hit sequences to amino acid sequences
- Align amino acid sequences to the reference one
- Calculate the amino acid identity to the reference in %

*NOTE: Other R packages used during testing are listed in `dependencies/R_packages.md`*

## BLASTx and nonsense mutation 

### Requirements

```
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
micromamba create -n codeathon
eval "$(micromamba shell hook --shell bash)"
micromamba activate codeathon
micromamba install python=3.12
micromamba install -y conda-forge::biopython bioconda::blast bioconda::diamond
```

### run_dgw.sh script

The `run_dgw.sh` script follows the following steps:
- Builds a protein BLAST reference database. 
- Executes a BLASTx scan of the query sequences against the protein reference database.

```bash
## for frameshift detection, diamond (v2.1.8) blastx has option 
## --frameshift to allow frameshift gap aligned and 
## show '/' or '\' in the alignement result

diamond blastx -q query.fna -d database/references.faa -f 5 --min-orf 1 --frameshift 15  -o results.xml
```

- Runs the python script to check nonsense mutation 

```bash
DGW.py -i test/input.xml 
```

```bash
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

