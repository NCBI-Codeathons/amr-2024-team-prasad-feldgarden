# Run on known sensitive and resistant isolates

## Download assemblies

```
datasets download genome accession --dehydrated --include genome --inputfile Pasens_GCAacc.txt
unzip ncbi_dataset.zip -d Pasens
datasets rehydrate --directory Pasens
cd Pasens
mv ncbi_dataset/data/*/*_genomic.fna .
cd ..

datasets download genome accession --include genome --inputfile Paresis_GCAacc.txt
unzip ncbi_dataset.zip -d Paresis
cd Paresis
mv ncbi_dataset/data/*/*_genomic.fna .
cd ..


datasets download genome accession --include genome --inputfile Kpsens_GCAacc.txt
unzip ncbi_dataset.zip -d Kpsens
cd Kpsens
mv ncbi_dataset/data/*/*_genomic.fna .
cd ..

datasets download genome accession --include genome --inputfile Kpresis_GCAacc.txt
unzip ncbi_dataset.zip -d Kpresis
cd Kpresis
mv ncbi_dataset/data/*/*_genomic.fna .
cd ..

```

