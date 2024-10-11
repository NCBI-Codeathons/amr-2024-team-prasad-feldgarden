# test on a random set of assemblies

## commands used to generate a random set of assemblies for testing

```
gcloud auth login --no-launch-browser

bq query --use_legacy_sql=false --format=csv --max_rows 300000 '
SELECT Run, computed_types.serotype, computed_types.antigen_formula,
    target_acc, biosample_acc, asm_acc, scientific_name
FROM `ncbi-pathogen-detect.pdbrowser.isolates` 
WHERE taxgroup_name = "Acinetobacter baumannii"
  AND asm_acc is not NULL
ORDER by rand() limit 1000
' > Ab_random.csv


bq query --use_legacy_sql=false --format=csv --max_rows 300000 '
SELECT Run, computed_types.serotype, computed_types.antigen_formula,
    target_acc, biosample_acc, asm_acc, scientific_name
FROM `ncbi-pathogen-detect.pdbrowser.isolates` tablesample system(50 percent)
WHERE taxgroup_name = "Klebsiella pneumoniae"
  AND asm_acc is not NULL
ORDER by rand() limit 1000
' > Kp_random.csv


bq query --use_legacy_sql=false --format=csv --max_rows 300000 '
SELECT Run, computed_types.serotype, computed_types.antigen_formula,
    target_acc, biosample_acc, asm_acc, scientific_name
FROM `ncbi-pathogen-detect.pdbrowser.isolates` tablesample system(50 percent)
WHERE taxgroup_name = "Pseudomonas aeruginosa"
  AND asm_acc is not NULL
ORDER by rand() limit 1000
' > Pa_random.csv
```

# Now download the data
```
cut -d',' -f6 Ab_random.csv | grep -v asm_acc > Ab_random.asm_acc
cut -d',' -f6 Kp_random.csv | grep -v asm_acc > Kp_random.asm_acc
cut -d',' -f6 Pa_random.csv | grep -v asm_acc > Pa_random.asm_acc

datasets download genome accession --dehydrated --include genome --inputfile Ab_random.asm_acc
unzip ncbi_dataset.zip -d Ab_genomes
datasets rehydrate --directory Ab_genomes


datasets download genome accession --dehydrated --include genome --inputfile Kp_random.asm_acc
unzip ncbi_dataset.zip -d Kp_genomes
datasets rehydrate --directory Kp_genomes


datasets download genome accession --dehydrated --include genome --inputfile Pa_random.asm_acc
unzip ncbi_dataset.zip -d Pa_genomes
datasets rehydrate --directory Pa_genomes

```

Now copy the files to get them in a place that's easier to work with

```
cd Ab_genomes/
mv ncbi_dataset/data/*/*_genomic.fna .

cd ../Kp_genomes
mv ncbi_dataset/data/*/*_genomic.fna .

cd ../Pa_genomes
mv ncbi_dataset/data/*/*_genomic.fna .
```

# Run DGW.py

```
pushd ../script
export PATH="$PATH:$PWD"
popd

mkdir Ab_dgw.py
for file in Ab_genomes/*
do
    base=`basename $file .fna`
    echo $file
    run_dgw.sh $file > Ab_dgw.py/$base.dgw.py
done


mkdir Kp_dgw.py
for file in Kp_genomes/*
do
    base=`basename $file .fna`
    echo $file
    run_dgw.sh $file > Kp_dgw.py/$base.dgw.py
done


mkdir Pa_dgw.py
for file in Pa_genomes/*
do
    base=`basename $file .fna`
    echo $file
    run_dgw.sh $file > Pa_dgw.py/$base.dgw.py
done


##### Run scanner.R

run_scanner.sh Kp_genomes.fa > Kp_genomes.scanner
run_scanner.sh Ab_genomes.fa > Ab_genomes.scanner
run_scanner.sh Pa_genomes.fa > Pa_genomes.scanner

