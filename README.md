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


## Approach

The solution developed in this project uses BLAST alignments to identify frameshift mutations and nonsense mutations. 

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

### Updates (September 25, 2024)

- Software has been improved and we have run tests on a large set of random isolates
- Put together a manually curated test set with sample output manually created
- Downloaded 30,000 assemblies to run our software on
- Need to refine output, standardize, and run on random isolates
- 

## Results

## Future Work

A second, more challenging goal (the 'reach project') will be to identify IS element insertions that affect AMR phenotypes, using ISAba3 and OXA-58 family carbapenemases in Acinetobacter baumannii as a test system.

The reach project will focus on the role of ISAba's in causing carbapenem resistance in *Acinetobacter baumannii*, specifically ISAba3 and blaOXA-58 family carbapenemases.

There are studies in the literature showing that when the ISAba3 is in the opposite orientation from the blaOXA-58, the carbapenemase is expressed (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3993105/). 

In addition, there have been reports ifd ISAba and ISAba-like elements disrupting AMR genes. 

### ISAba Papers of interest:

1. For A.b., this paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9771954/ has a lot of good examples for the ISAba project.
2. Some additional background papers:
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC538857/
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10150277/
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4517069/
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6635527/
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2916316/
  - https://pubmed.ncbi.nlm.nih.gov/16630258/
  - https://pubmed.ncbi.nlm.nih.gov/16441449/
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3498117/ -this is an interrupted AdeS gene in A. baumannii
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10655397/
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10253819/ also IS-interrupted genes
  - https://pubmed.ncbi.nlm.nih.gov/16441449/
  - https://pubmed.ncbi.nlm.nih.gov/32712382/


## NCBI Codeathon Disclaimer
This software was created as part of an NCBI codeathon, a hackathon-style event focused on rapid innovation. While we encourage you to explore and adapt this code, please be aware that NCBI does not provide ongoing support for it.

For general questions about NCBI software and tools, please visit: [NCBI Contact Page](https://www.ncbi.nlm.nih.gov/home/about/contact/)

