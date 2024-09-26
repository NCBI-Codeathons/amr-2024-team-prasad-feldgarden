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

# run DGW on test set

```Shell
mkdir Pasens.dgw
for file in Pasens/*.fna
do
    base=`basename $file .fna`
    echo $file
    run_dgw.sh $file > Pasens.dgw/$base.dgw
done

mkdir Paresis.dgw
for file in Paresis/*.fna
do
    base=`basename $file .fna`
    echo $file
    run_dgw.sh $file > Paresis.dgw/$base.dgw
done

mkdir Kpsens.dgw
for file in Kpsens/*.fna
do
    base=`basename $file .fna`
    echo $file
    run_dgw.sh $file > Kpsens.dgw/$base.dgw
done

mkdir Kpresis.dgw
for file in Kpresis/*.fna
do
    base=`basename $file .fna`
    echo $file
    run_dgw.sh $file > Kpresis.dgw/$base.dgw
done

