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
```

Redone as:
```
rename_contig.pl Kpsens/*.fna > Kpsens.fna
run_dgw.sh Kpsens.fna > Kpsens.dgw1
R --no-save <<END
library(dplyr)
library(tidyr)
tabin('Kpsens.dgw1') %>% separate(contig_acc, into=c('asm_acc','contig_acc'), sep='-', remove=F) %>% tabout('Kpsens.dgw2')
END

rename_contig.pl Kpresis/*.fna > Kpresis.fna
run_dgw.sh Kpresis.fna > Kpresis.dgw1
R --no-save <<END
library(dplyr)
library(tidyr)
tabin('Kpresis.dgw1') %>% separate(contig_acc, into=c('asm_acc','contig_acc'), sep='-', remove=F) %>% tabout('Kpresis.dgw2')
END

rename_contig.pl Pasens/*.fna > Pasens.fna
run_dgw.sh Pasens.fna > Pasens.dgw1
R --no-save <<END
library(dplyr)
library(tidyr)
tabin('Pasens.dgw1') %>% separate(contig_acc, into=c('asm_acc','contig_acc'), sep='-', remove=F) %>% tabout('Pasens.dgw2')
END

rename_contig.pl Paresis/*.fna > Paresis.fna
run_dgw.sh Paresis.fna > Paresis.dgw1
R --no-save <<END
library(dplyr)
library(tidyr)
tabin('Paresis.dgw1') %>% separate(contig_acc, into=c('asm_acc','contig_acc'), sep='-', remove=F) %>% tabout('Paresis.dgw2')
END

```

cp Kpsens.dgw2 ../../../results/Kpsens.dgw
cp Kpresis.dgw2 ../../../results/Kpresis.dgw
cp Pasens.dgw2 ../../../results/Pasens.dgw
cp Paresis.dgw2 ../../../results/Paresis.dgw

```

