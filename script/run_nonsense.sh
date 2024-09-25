#!/bin/bash

if [ "$1" == "" ]
then
    echo "run_nonsense.sh - identify nonsense mutations"
    echo "Usage: run_nonsense.sh <NUCLEOTIDE_FASTA>"
    exit 1
fi

script_dir=`dirname $0`
DATABASE="$script_dir/database/references.faa"
tmpdir=`mktemp -d`


input_fasta="$1"

blastx -query $input_fasta -subject $DATABASE -out $tmpdir/blast.xml -outfmt 5

$script_dir/check_nonsense_mutations.py -i $tmpdir/blast.xml

rm -r $tmpdir

