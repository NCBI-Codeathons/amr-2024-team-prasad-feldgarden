![ncbilogo](https://github.com/user-attachments/assets/4b2da250-0b29-4298-8a04-dcc8e1b61e92)
# Dead gene walking

List of participants and affiliations:
- Arjun Prasad, Team co-Lead
- Michael Feldgarden, Team co-Lead
- Adrien Assie, Tech Lead
- EB Dickinson
- Chienchi Lo
- Ana Ramos, Writer

## Project Goals

This project aims to develop and test methods to identify the presence of loss-of-function mutations in bacterial genes that cause antibiotic resistance. The first part of the project will be to develop a tool to identify loss-of-function mutations, such as stop codons and frameshifts, in genes where these mutations have been demonstrated to affect function (e.g., OmpK35/K36 in _Klebsiella pneumoniae_).

## Project Output

Two scripts have been developed to identify frameshift mutations and nonsense mutations: [scanner.R](https://github.com/NCBI-Codeathons/amr-2024-team-prasad-feldgarden/blob/main/script/scanner.R) and [DGW.py](https://github.com/NCBI-Codeathons/amr-2024-team-prasad-feldgarden/blob/main/script/DGW.py).

A wrapper shell script, [run_DGW.sh](https://github.com/NCBI-Codeathons/amr-2024-team-prasad-feldgarden/blob/main/script/run_DGW.sh), was developed to run the software across a large set of genomes.

A full description of the software, their use and installation instructions can be found in the [script README](script/README.md).

The slides for the [codeathon final presentations](https://www.nlm.nih.gov/ncbi/workshops/2024-09_AMR-Codeathon/final-presentations.html) are [here](/Codeathon%20presentation/DeadGeneWalking.pptx.pdf).

## Dev log:

### Scope of this project:

1. Identifying and collecting sequences for validated examples of positives (broken genes that lead to resistance), and negatives (non-broken genes and/or genes broken that are not resistant). The identified sequences can be found in the `data` directory, which includes sequences previously identified by the team co-leads. Some examples include:
    - *K. pneumoniae* genes:  
      - OmpK35 stop in SAMN32518855
      - circA frameshift in SAMN31181384
    - *P. aeruginosa* genes: 
      - ampD stop in SAMN11110537
      - mexR frameshift in SAMN11110706
      - nalD frameshift in SAMN11110430

Specifically to start we will look for lesions in:

- _Acinetobacter baumanii_
  - AdeS
- _Klebsiella pneumoniae_
  - OmpK35
  - OmpK36
  - CirA
- _Pseudomonas aeruginosa_
  - OprD
  - AmpD
  - NalD
  - MexR

2. Developing software to identify broken genes caused by frameshift mutations and nonsense mutations.

## Structure of this repository

```
|-- Codeathon\ presentation
|-- data - Data used for analysis and some analysis results
|   |-- dev-examples - test data we have manually inspected
|   |   |-- resis_and_sens_Kp_and_Pa - Set 
|   |   `-- testset - A compilation of the manually inspected test data to use for software development
|   |-- randomset1 - A set of 1000 random isolates from each of the three taxonomic groups investigated here
|   |-- randomset2 - A set of 10,000 random isolates from each of the three taxonomic groups
|   `-- resis_and_sens_Kp_and_Pa - A set of isolates with AST data segregated into resistant and sensitive categories
|-- results - Results of our analyses
`-- script - Software and database
    |-- database - The reference database we use
    |-- dependencies - Dependencies for installing the software
    `-- test - Test data for software
```

## Approach

The solution developed in this project uses multiple BLAST alignments to identify frameshift mutations and nonsense mutations. 

### Assembly sequence to identify broken genes

- Frameshift detection
  - Nucleotide alignment with BLASTN and post-processing to identify frame shifts
      - The difference in identity (%) between nucleotide and translated alignment is used to identify frame-shifted proteins
- Stop codon detection
  - Translated alignment with BLASTX to identify stop codons

### Testing     
- Use test sequences in [`data/dev-examples/`](https://github.com/NCBI-Codeathons/amr-2024-team-prasad-feldgarden/tree/main/data/dev-examples) to make sure the software is working
    - So far we have collected > 90 example mutations for six genes from three taxa
    - The majority of the reference sequences for AMR genes are found in NCBI's [Reference Gene Catalog](https://www.ncbi.nlm.nih.gov/pathogens/refgene/#)
- Run the software across a large selection of genomes to characterize gene disruption and AMR
- In cases where BLASTN alignments needed to be confirmed, [Vadr](https://github.com/ncbi/vadr?tab=readme-ov-file#vadr---viral-annotation-definer-) was used for QA comparisons. 
  - Vadr has proven effective when running against specific bacterial genes. For now, we will use Vadr as a "control" to help us compare results from the tools we are developing

### Output format

The output is provided in a tabular form. An example of the output schema is provided [here](https://github.com/NCBI-Codeathons/amr-2024-team-prasad-feldgarden/blob/main/mock_output.csv). 

### Updates (September 24, 2024)

- Broken genes for testing have been identified. Reference sequences and assembly sequences have been obtained for those
  - We will now build a BLAST database with these (test database)
- The assembly sequence will rely on BLAST alignments to obtain both the [frameshift](https://github.com/NCBI-Codeathons/amr-2024-team-prasad-feldgarden/blob/main/script/scanner.R) (BLASTN) and identify [stop codons](https://github.com/NCBI-Codeathons/amr-2024-team-prasad-feldgarden/blob/main/script/check_nonsense_mutations.py) (BLASTX) lessions. 
  - The [output format](https://github.com/NCBI-Codeathons/amr-2024-team-prasad-feldgarden/blob/main/output_format.md) was defined to ensure compatibility of our tools
  - Both scripts will be tested against the test database to confirm the same BLASTN parameters can be used to identify both lessions
  - Once confirmed, they will be added to a bash script
  - Vadr has proven effective when running against specific bacterial genes. For now, we will use Vadr as a "control" to help us compare results from the tools we are developing
- A virtual machine was built to assist in tool testing across environments

### Updates (September 26, 2024)

- Software has been improved and we have run tests on a large set of random isolates
- Put together a manually curated test set with sample output manually created
- Downloaded 30,000 assemblies to run our software on
- Need to refine output, standardize, and run on random isolates

## Results

### AST testset

A test set derived from AST consisting of 928 isolates was analyzed by DGW. The AST dataset included:
  - 359 *Pseudomonas aeruginosa* isolates
    - 292 with resistant phenotype
    - 67 with sensitive phenotype
  - 569 *Klebsiella pneumoniae* isolates
    - 297 with resistant phenotype
    - 272 with sensitive phenotype
    
Of the 589 resistant isolates DGW was able to identify 493 mutations, of which 130 are present in MicroBIGG-E. 

DWG was able to identify 363 new mutations in resistant isolates. 

A more comprehensive analysis of the results can be found in [`results/Analisys-of-DGW-results--AST-set-.html`](https://github.com/NCBI-Codeathons/amr-2024-team-prasad-feldgarden/blob/main/results/Analysis-of-DGW-results.html)


### Random testset

A random set of 3000 genome accession IDs for *Acinetobacter baumannii*, *Pseudomonas aeruginosa* and *Klebsiella pneumoniae* (1000 isolates for each) were obtained using NCBI’s datasets as described [here](https://github.com/NCBI-Codeathons/amr-2024-team-prasad-feldgarden/tree/main/results).

DGW identified 877 mutations out of the 3000 isolates analyzed
 - 150 mutations were present in MicroBIGG-E too
 - 727 new mutations in 353 isolates

A more comprehensive analysis of the results can be found in [`results/Analisys-of-DGW-results--randomset-.html`](https://github.com/NCBI-Codeathons/amr-2024-team-prasad-feldgarden/tree/main/results).

### Random testset2

A random set of 30,000 genome accession IDs for *Acinetobacter baumannii*, *Pseudomonas aeruginosa* and *Klebsiella pneumoniae* (10,000 isolates for each) were obtained using NCBI’s datasets as described [here](https://github.com/NCBI-Codeathons/amr-2024-team-prasad-feldgarden/tree/main/results).

- DGW identified 8377 mutations out of the 30,000 isolates analyzed
  - 1361 mutations were present in MicroBIGG-E too
  - 7016 new mutations in 353 isolates

A more comprehensive analysis of the results can be found in [`results/Analisys-of-DGW-results--randomset-.html`](https://github.com/NCBI-Codeathons/amr-2024-team-prasad-feldgarden/tree/main/results).

  
## Future Work

With the promising results from the DGW tool produced during the codeathon, future work should focus on expanding the reference database of genes reported in the to lead to AMR after a loss-of-function mutation. 

## NCBI Codeathon Disclaimer
This software was created as part of an NCBI codeathon, a hackathon-style event focused on rapid innovation. While we encourage you to explore and adapt this code, please be aware that NCBI does not provide ongoing support for it.

For general questions about NCBI software and tools, please visit: [NCBI Contact Page](https://www.ncbi.nlm.nih.gov/home/about/contact/)

