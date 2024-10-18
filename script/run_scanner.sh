#!/bin/bash

if [ "$1" == "" ]
then
    echo "run_scanner.sh - identify nonsense mutations and frame shifts"
    echo "Usage: run_scanner.sh <NUCLEOTIDE_FASTA>"
    echo "Uses scanner.R and you will need to have the prerequisites installed"
    exit 1
fi

script_dir=`dirname $0`
DATABASE="$script_dir/database/references.fna"
db_base=`basename $DATABASE .faa`
db_dir=`dirname $DATABASE`
tmpdir=`mktemp -d`
input_fasta="$1"
input_base="`basename $input_fasta`.fasta"

scanner.R -i "$input_fasta" -o "$tmpdir" -d "$DATABASE" -p "$input_base" -c 16 > $tmpdir/output 2> $tmpdir/error
if [ "$?" -gt 0 ]
then
    cat $tmpdir/output
    >&2 cat $tmpdir/error
    exit 1
fi
cat $tmpdir/*.fasta_results.tsv

rm -r $tmpdir
