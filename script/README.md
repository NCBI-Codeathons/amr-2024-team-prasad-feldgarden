# Script folder

## Blastn and translation comparison approach

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

Starts with a blastn scan of the query sequences against the reference database.

```js
blastx -query query.fasta -db reference_db -out results.xml -outfmt 5
```

```js
check_nonsense_mutations.py -i test/input.xml 
```
