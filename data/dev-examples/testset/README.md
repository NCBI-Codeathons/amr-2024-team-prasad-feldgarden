# Test sequences

- `SAMN11110706.nuc.fa` is a Pseudomonas genome that has a MexR frameshift after pos. 107 (aa); in nucleotide space, in the query, pos. 116-126 are missing, and in the query this corresponds to pos. 245917-245925 in DAHMTF010000002.1
    - It looks to me (Arjun) like the lesion appears at position 125 in the reference and position 335615 in DAHMTF010000002.1 (it's on the negative strand)
    - `SAMN11110706.expected` is a mockup of output based on the format we discussed. It may not be perfect

- `SAMN11110537.nuc.fa` is a Pseudomonas that has an ampD stop at pos. 108 E_to stop; pos. 37418 in DAHNQV010000008.1
    - `SAMN11110537.expected` is a mockup of output baesd on the format we discussed.

- `synthetic_fs.fna` - Manually mutated ompK35 sequence with an insertion after site 300 in the coding sequence (100 in the amino-acid sequence)
    - `synthetic_fs.expected` is a mockup of expected output

- `OmpK35_frameshift.nuc.fa` - OmpK35 frameshift sequence from test data for Kleborate. A frameshift at pos. 116; in nucleotide space, pos. 345 in the reference has an insertion (A), pos. 3084 in the test sequence 
    - `OmpK35_frameshift.expected` is a mockup of expected output

- `OmpK35_early_stop.nuc.fa` - OmpK35 early stop from test data for Kleborate. Nonsense mutation (early stop) at position 36
    - `OmpK35_early_stop.nuc.expected` mockup of expected output


----------

## The following are not ready yet

- `SAMN03105473.nuc.fa` is a Pseudomonas genome that has a nalD frameshift at pos. 131; in nucleotide space, there is a deletion at 394-397 in the ref. and pos. 59736-59739 in JTXK01000022.1
    - I couldn't find the nalD mutation ?

I could not find the following:
- The two ompK35 sequences are taken from the test data for Kleborate (https://github.com/klebgenomics/Kleborate/tree/main/test/test_res_omp); these are found in actualKP_testgenomes:
1. An early stop at pos. 36 (Y to *) ; in nucleotide space, pos. 1035 has C to T in the reference, 3220 in the test sequence
2. A frameshift at pos. 116; in nucleotide space, pos. 345 in the reference has an insertion (A), pos. 3084 in the test sequence
