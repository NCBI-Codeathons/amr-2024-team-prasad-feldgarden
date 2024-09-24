# Script folder

## Require blast and R

```
sudo apt-get install r-base
sudo apt-get install ncbi-blast+
```

## Blastn and translation comparison approach


```
R
install.packages('tidyverse')
install.packages('XML')
install.packages('Biostrings')
```

Starts with a blastn scan of the query sequences against the reference database.

```
blastn -query query.fasta -db reference_db -out results.xml -outfmt 5
```

The `scanner.R` script is used to import the data from the BLAST scan and parse it.
The script:
- extract the nucleotide hit sequences
- calculate the nucleotide identity to the reference in %,
- translate nucleotide hit sequences to Amino acid sequences
- align amino acid sequences to the reference one
- calculate the amino acid identity to the reference in %,
- merge with metadata
- plot the data

## Blastx and nonsense mutation 

Requires Bio
```
micromamba create -n codeathon
eval "$(micromamba shell hook --shell bash)"
micromamba activate codeathon
micromamba install python=3.12
micromamba install Bio
```


Starts with a blastx scan of the query sequences against the reference protein database.

```js
blastx -query query.fasta -db reference_db -out results.xml -outfmt 5
```

Run the python script to check nonsense mutation 

```js
check_nonsense_mutations.py -i test/input.xml 
```
