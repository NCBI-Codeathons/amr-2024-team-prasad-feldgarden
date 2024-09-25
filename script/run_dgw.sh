#!/bin/bash

if [ "$1" == "" ]
then
    echo "run_dgw.sh - identify nonsense mutations and frame shifts"
    echo "Usage: run_dgw.sh <NUCLEOTIDE_FASTA>"
    echo "Runs diamond then DGW.py on the output"
    exit 1
fi

script_dir=`dirname $0`
DATABASE="$script_dir/database/references.faa"
db_base=`basename $DATABASE .faa`
db_dir=`dirname $DATABASE`
if [ ! -e "$db_dir/$db_base.dmnd" ]
then
    diamond makedb --in $DATABASE --db "$db_dir/$db_base"
fi
tmpdir=`mktemp -d`

input_fasta="$1"

diamond blastx --query $input_fasta --db $DATABASE --out $tmpdir/blast.xml --frameshift 15 --min-orf 1 -f 5

micromamba run -n codeathon $script_dir/DGW.py -i $tmpdir/blast.xml


rm -r $tmpdir

