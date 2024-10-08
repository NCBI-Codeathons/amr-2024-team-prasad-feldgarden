# Comparison of results from `DGW.py` and `scanner.R` against Microbigg-E subset

1. Combine output files

The output file for each genome in the AST testset was combined into a single `.tsv` file using the `results_combined.R` script. 

The `*_combined_scriptR.tsv` files were then used as input to analyze the results.

2. Combine outputs from `scanner.R` and `DGW.py` with `microbigge_subset`

Update the SQL model to import the correct output files (from step 1). Then run:

```
 duckdb < results/script_comparison/dgw_combined_outputs.sql
```

3. Analyze results in Rmd notebook (e.g. `ASTset_results.Rmd`)
