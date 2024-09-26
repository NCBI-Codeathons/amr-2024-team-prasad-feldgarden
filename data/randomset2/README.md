```
gcloud auth login --no-launch-browser

bq query --use_legacy_sql=false --format=csv --max_rows 300000 '
SELECT Run, computed_types.serotype, computed_types.antigen_formula,
    target_acc, biosample_acc, asm_acc, scientific_name
FROM `ncbi-pathogen-detect.pdbrowser.isolates`
WHERE taxgroup_name = "Acinetobacter baumannii"
  AND asm_acc is not NULL
ORDER by rand() limit 10000
' > Ab_random.csv

bq query --use_legacy_sql=false --format=csv --max_rows 300000 '
SELECT Run, computed_types.serotype, computed_types.antigen_formula,
    target_acc, biosample_acc, asm_acc, scientific_name
FROM `ncbi-pathogen-detect.pdbrowser.isolates` tablesample system(50 percent)
WHERE taxgroup_name = "Klebsiella pneumoniae"
  AND asm_acc is not NULL
ORDER by rand() limit 10000
' > Kp_random.csv

bq query --use_legacy_sql=false --format=csv --max_rows 300000 '
SELECT Run, computed_types.serotype, computed_types.antigen_formula,
    target_acc, biosample_acc, asm_acc, scientific_name
FROM `ncbi-pathogen-detect.pdbrowser.isolates` tablesample system(50 percent)
WHERE taxgroup_name = "Pseudomonas aeruginosa"
  AND asm_acc is not NULL
ORDER by rand() limit 10000
' > Pa_random.csv

# download

cut -d',' -f6 Ab_random.csv | grep -v asm_acc > Ab_random.asm_acc
cut -d',' -f6 Kp_random.csv | grep -v asm_acc > Kp_random.asm_acc
cut -d',' -f6 Pa_random.csv | grep -v asm_acc > Pa_random.asm_acc

datasets download genome accession --dehydrated --include genome --inputfile Ab_random.asm_acc
unzip ncbi_dataset.zip -d Ab_genomes
datasets rehydrate --directory Ab_genomes

datasets download genome accession --dehydrated --include genome --inputfile Kp_random.asm_acc --filename Kp_datasets.zip
unzip Kp_datasets.zip -d Kp_genomes
datasets rehydrate --directory Kp_genomes

datasets download genome accession --dehydrated --include genome --inputfile Pa_random.asm_acc --filename Pa_datasets.zip
unzip Pa_datasets.zip -d Pa_genomes
datasets rehydrate --directory Pa_genomes

