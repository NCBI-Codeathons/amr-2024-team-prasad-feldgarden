#!/bin/sh

script_dir=`dirname $0`
DATABASE="$script_dir/database/references.faa"
tmpdir=`mktemp`


my $input_fasta=$1

blastx -query $input_fasta -subject $DATABASE -out $tmpdir/blast.xml -outfmt 5

$script_dir/check_nonsense_mutations.py -i $tmpdir/blast.xml

rm -r $tmpdir

