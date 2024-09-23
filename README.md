# Team Project Name TBD

List of participants and affiliations:
- Arjun Prasad, Team co-Lead
- Michael Feldgarden, Team co-Lead
- Adrien Assie, Tech Lead
- Emma Dickinson
- Chienchi Lo
- Ana Ramos, Writer
- Erin Young

## Project Goals

This project aims to develop and test methods to identify the presence of loss-of-function mutations in bacterial genes that cause antibiotic resistance. The first part of the project will be to develop a tool to identify loss-of-function mutations, such as stop codons and frameshifts, in genes where these mutations have been demonstrated to affect function (e.g., OmpK35/K36 in _Klebsiella pneumoniae_). A second, more challenging goal (the 'reach project') will be to identify IS element insertions that affect AMR phenotypes, using ISAba3 and OXA-58 family carbapenemases in _Acinetobacter baumannii_ as a test system.

### There will be three aims for this project:

1. Identifying and collecting sequences for validated examples of positives (broken genes that lead to resistance), and negatives (non-broken genes and/or genes broken that are not resistant) for the examples above. This will require some literature search, and often some digging to find the actual sequences and make sure they’re really what they say they are. Mike has done a little work in this direction, but we should have more examples for validation.
2. Developing software to identify broken genes (frame-shift mutations) and nonsense mutations. Most likely by running BLAST and post-processing the results.
3. If people are interested, identifying and coming up with methods to identify the ISAba or other similar insertion element mutations that cause resistance (will require both identifying more examples and developing methods, but method development and testing will need to be a little more involved)

### Some of the genes we could focus on are:

- OmpK35/K36 lesions (stops, frameshifts etc.) in Klebsiella pneumoniae.
  - OmpK35/K36 refs: ADG27466.1/ADG56549.1 (note that many reference sequences for AMR genes will be found in NCBI's Reference Gene Catalog https://www.ncbi.nlm.nih.gov/pathogens/refgene/#)
- P. aeruginosa genes: 
  - ampD stops: example genome is SAMN11110537
  - mexR frameshifts: example genome is SAMN11110706
  - nalD frameshifts: example genome is SAMN11110430

### Reach project information:

The reach project will be the role of ISAba's in causing carbapenem resistance in Acinetobacter baumannii.  Specifically:

1. We will look at ISAba3 and blaOXA-58 family carbapenemases.
2. When the ISAba3 is in the opposite orientation from the blaOXA-58, the carbapenemase is expressed (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3993105/)
3. A second reach project could be to look for ISAba that disrupt AMR genes.

### Papers of interest:

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

## Approach

## Results

## Future Work

## NCBI Codeathon Disclaimer
This software was created as part of an NCBI codeathon, a hackathon-style event focused on rapid innovation. While we encourage you to explore and adapt this code, please be aware that NCBI does not provide ongoing support for it.

For general questions about NCBI software and tools, please visit: [NCBI Contact Page](https://www.ncbi.nlm.nih.gov/home/about/contact/)

