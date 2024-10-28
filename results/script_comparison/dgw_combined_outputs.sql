SET temp_directory='results/tmp.tmp';

.mode csv
.output results/script_comparison/dgw_combined_outputs_ast.csv

WITH 
-- Replace the DGW.py output file in each union block as required
dgw_combined_p AS (
    SELECT 
    *,
    md5(CONCAT(LOWER(element_symbol), asm_acc, contig_acc, lesion)) AS mutation_id,
    'kp_resistant' AS cohort
    FROM read_csv('results/Kpresis.dgw')
    UNION ALL 
    SELECT 
    *,
    md5(CONCAT(LOWER(element_symbol), asm_acc, contig_acc, lesion)) AS mutation_id,
    'pa_resistant' AS cohort
    FROM read_csv('results/Paresis.dgw')
    UNION ALL 
    SELECT 
    *,
    md5(CONCAT(LOWER(element_symbol), asm_acc, contig_acc, lesion)) AS mutation_id,
    'kp_sensitive' AS cohort
    FROM read_csv('results/Kpsens.dgw')
    UNION ALL 
    SELECT 
    *,
    md5(CONCAT(LOWER(element_symbol), asm_acc, contig_acc, lesion)) AS mutation_id,
    'pa_sensitive' AS cohort
    FROM read_csv('results/Pasens.dgw')
),

-- Replace the scanner.R output file in each union block as required
dgw_combined_r AS (
    SELECT 
    *,
    md5(CONCAT(LOWER(element_symbol), asm_acc, contig_acc, lesion)) AS mutation_id,
    'kp_resistant' AS cohort
    FROM 'results/scanner.R_AST/Kpresis_combined_scriptR.tsv'
    UNION ALL 
    SELECT 
    *,
    md5(CONCAT(LOWER(element_symbol), asm_acc, contig_acc, lesion)) AS mutation_id,
    'pa_resistant' AS cohort
    FROM 'results/scanner.R_AST/Paresis_combined_scriptR.tsv' 
    UNION ALL 
    SELECT 
    *,
    md5(CONCAT(LOWER(element_symbol), asm_acc, contig_acc, lesion)) AS mutation_id,
    'kp_sensitive' AS cohort
    FROM 'results/scanner.R_AST/Kpsens_combined_scriptR.tsv'
    UNION ALL 
    SELECT 
    *,
    md5(CONCAT(LOWER(element_symbol), asm_acc, contig_acc, lesion)) AS mutation_id,
    'pa_sensitive' AS cohort
    FROM read_csv('results/scanner.R_AST/Pasens_combined_scriptR.tsv')
),

-- Import subset of Microbigg-E results
microbigge AS (
    SELECT * FROM 'results/microbigge_subset_ast.csv'
),

-- Get list of all the genome accession numbers in the testset with identified lesions
all_acc AS 
(
SELECT DISTINCT asm_acc FROM 
  (SELECT DISTINCT asm_acc FROM dgw_combined_p UNION ALL 
  SELECT DISTINCT asm_acc FROM dgw_combined_r UNION ALL 
  SELECT DISTINCT asm_acc FROM microbigge)
),

p_mutations AS (
  SELECT DISTINCT
    all_acc.asm_acc,
    dgw_combined_p.element_symbol AS p_element_symbol,
    dgw_combined_p.lesion AS p_lesion,
    dgw_combined_p.lesion_type AS p_lesion_type,
    dgw_combined_p.mutation_id AS p_mutation_id,
    CONCAT(all_acc.asm_acc, LOWER(dgw_combined_p.element_symbol)) AS p_mutations_concat,
    COUNT(DISTINCT dgw_combined_p.mutation_id) OVER (PARTITION BY all_acc.asm_acc) AS p_mut_count,
    dgw_combined_p.cohort
  FROM all_acc
  LEFT JOIN dgw_combined_p
  ON all_acc.asm_acc = dgw_combined_p.asm_acc
),

r_mutations AS (
  SELECT DISTINCT
    all_acc.asm_acc,
    dgw_combined_r.element_symbol AS r_element_symbol,
    dgw_combined_r.lesion AS r_lesion,
    dgw_combined_r.lesion_type AS r_lesion_type,
    dgw_combined_r.mutation_id AS r_mutation_id,
    CONCAT(all_acc.asm_acc, LOWER(dgw_combined_r.element_symbol)) AS r_mutations_concat,
    COUNT(DISTINCT dgw_combined_r.mutation_id) OVER (PARTITION BY all_acc.asm_acc) AS r_mut_count,
    dgw_combined_r.cohort
  FROM all_acc
  LEFT JOIN dgw_combined_r
  ON all_acc.asm_acc = dgw_combined_r.asm_acc
),

combined AS (
  SELECT DISTINCT
    all_acc.asm_acc,
    -- p_mutations.p_mutations_concat,
    -- r_mutations.r_mutations_concat,
    p_mutations.p_element_symbol,
    r_mutations.r_element_symbol,
    p_mutations.p_lesion,
    r_mutations.r_lesion,
    p_mutations.p_lesion_type,
    r_mutations.r_lesion_type,
    p_mutations.p_mutation_id,
    r_mutations.r_mutation_id,
    STRING_AGG(DISTINCT(p_mutations.cohort, r_mutations.cohort)) OVER (PARTITION BY all_acc.asm_acc) AS cohort,
    p_mutations.p_mut_count,
    r_mutations.r_mut_count
  FROM all_acc
  LEFT JOIN p_mutations
  ON all_acc.asm_acc = p_mutations.asm_acc
  LEFT JOIN r_mutations
  ON all_acc.asm_acc = r_mutations.asm_acc
  AND p_mutations.p_mutations_concat = r_mutations.r_mutations_concat
  WHERE p_mut_count > 0

  UNION BY NAME

  SELECT DISTINCT
    all_acc.asm_acc,
    r_mutations.r_element_symbol,
    r_mutations.r_lesion,
    r_mutations.r_lesion_type,
    r_mutations.r_mutation_id,
    STRING_AGG(DISTINCT(p_mutations.cohort, r_mutations.cohort)) OVER (PARTITION BY all_acc.asm_acc) AS cohort,
    r_mutations.r_mut_count
  FROM all_acc
  LEFT JOIN p_mutations
  ON all_acc.asm_acc = p_mutations.asm_acc
  LEFT JOIN r_mutations
  ON all_acc.asm_acc = r_mutations.asm_acc
  WHERE p_mutations.p_mutations_concat != r_mutations.r_mutations_concat
  AND p_mutations.p_mut_count > 0
  AND LOWER(r_element_symbol) NOT IN (SELECT DISTINCT LOWER(p_element_symbol) FROM p_mutations WHERE asm_acc = r_mutations.asm_acc)

  UNION BY NAME

  SELECT DISTINCT
    all_acc.asm_acc,
    r_mutations.r_element_symbol,
    r_mutations.r_lesion,
    r_mutations.r_lesion_type,
    r_mutations.r_mutation_id,
    STRING_AGG(DISTINCT(p_mutations.cohort, r_mutations.cohort)) OVER (PARTITION BY all_acc.asm_acc) AS cohort,
    r_mutations.r_mut_count
  FROM all_acc
  LEFT JOIN p_mutations
  ON all_acc.asm_acc = p_mutations.asm_acc
  LEFT JOIN r_mutations
  ON all_acc.asm_acc = r_mutations.asm_acc
  WHERE p_mut_count = 0
) , 

mbe AS (
  SELECT 
    asm_acc,
    SPLIT_PART(element_symbol, '_', 1) AS mb_element_symbol,
    SPLIT_PART(element_symbol, '_', 2) AS mb_lesion,
  FROM microbigge
  WHERE LOWER(microbigge.element_symbol) LIKE '%oprd%' OR 
        LOWER(microbigge.element_symbol) LIKE '%ampd%' OR
        LOWER(microbigge.element_symbol) LIKE '%ompk%' OR
        LOWER(microbigge.element_symbol) LIKE '%cira%' OR
        LOWER(microbigge.element_symbol) LIKE '%nald%' OR
        LOWER(microbigge.element_symbol) LIKE '%mexr%'
),

combined_2 AS (
  SELECT DISTINCT
    combined.asm_acc,
    combined.p_element_symbol,
    combined.r_element_symbol,
    mbe.mb_element_symbol,
    CASE 
      WHEN LOWER(mbe.mb_element_symbol) IN (LOWER(combined.p_element_symbol), LOWER(combined.r_element_symbol)) THEN 1
      WHEN COALESCE(LOWER(combined.p_element_symbol), LOWER(combined.r_element_symbol)) IS NULL AND mbe.mb_element_symbol IS NOT NULL THEN 1
      WHEN mbe.mb_element_symbol IS NOT NULL THEN 0
      ELSE 3
    END AS mb_element_symbol_check,
    combined.p_lesion,
    combined.r_lesion,
    mbe.mb_lesion,
    combined.p_lesion_type,
    combined.r_lesion_type,
    combined.p_mutation_id,
    combined.r_mutation_id,
    cohort
  FROM combined
  LEFT JOIN mbe
  ON combined.asm_acc = mbe.asm_acc
  AND LOWER(mbe.mb_element_symbol) IN (LOWER(combined.p_element_symbol), LOWER(combined.r_element_symbol))
  QUALIFY COUNT(COALESCE(combined.p_element_symbol, combined.r_element_symbol)) OVER (PARTITION BY combined.asm_acc) > 0

  UNION BY NAME

  SELECT DISTINCT
    combined.asm_acc,
    mbe.mb_element_symbol,
    CASE 
      WHEN LOWER(mbe.mb_element_symbol) IN (LOWER(combined.p_element_symbol), LOWER(combined.r_element_symbol)) THEN 1
      WHEN COALESCE(LOWER(combined.p_element_symbol), LOWER(combined.r_element_symbol)) IS NULL AND mbe.mb_element_symbol IS NOT NULL THEN 1
      WHEN mbe.mb_element_symbol IS NOT NULL THEN 0
      ELSE 3
    END AS mb_element_symbol_check,
    mbe.mb_lesion,
    cohort
  FROM combined
  LEFT JOIN mbe
  ON combined.asm_acc = mbe.asm_acc
  QUALIFY COUNT(COALESCE(combined.p_element_symbol, combined.r_element_symbol)) OVER (PARTITION BY combined.asm_acc) = 0
    UNION BY NAME

  SELECT DISTINCT
    combined.asm_acc,
    mbe.mb_element_symbol,
    CASE 
      WHEN LOWER(mbe.mb_element_symbol) IN (LOWER(combined.p_element_symbol), LOWER(combined.r_element_symbol)) THEN 1
      WHEN COALESCE(LOWER(combined.p_element_symbol), LOWER(combined.r_element_symbol)) IS NULL AND mbe.mb_element_symbol IS NOT NULL THEN 1
      WHEN mbe.mb_element_symbol IS NOT NULL THEN 0
      ELSE 3
    END AS mb_element_symbol_check,
    mbe.mb_lesion,
    cohort
  FROM combined
  LEFT JOIN mbe
  ON combined.asm_acc = mbe.asm_acc
  WHERE mbe.mb_element_symbol IS NOT NULL
  AND LOWER(mbe.mb_element_symbol) NOT IN (SELECT DISTINCT LOWER(p_element_symbol) FROM combined WHERE asm_acc = mbe.asm_acc)
  AND LOWER(mbe.mb_element_symbol) NOT IN (SELECT DISTINCT LOWER(r_element_symbol) FROM combined WHERE asm_acc = mbe.asm_acc)
  QUALIFY COUNT(COALESCE(combined.p_element_symbol, combined.r_element_symbol)) OVER (PARTITION BY combined.asm_acc) > 0
),

final AS (
SELECT 
  *,
  md5(CONCAT(asm_acc, p_element_symbol, r_element_symbol, mb_element_symbol, p_lesion, r_lesion, mb_lesion, p_lesion_type, r_lesion_type, p_mutation_id, r_mutation_id)) AS record_id,  
FROM combined_2
ORDER BY asm_acc
) 

SELECT * FROM final;