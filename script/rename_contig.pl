#!/bin/env perl
use strict;
use warnings;
use File::Basename;

if (@ARGV == 0) {
    print "rename_contig.pl - rename contigs to add the assembly accession to the\n";
    print "start of the contig identifier\n";
    print 'Uses s/(\.\d+)_.*/$1/; to identify the assembly accession from the filename';
    print "\nUsage: rename_contig.pl <FASTA_FILE> [<FASTA_FILE> ...] > <combined.fa>\n";
    exit 1;
}

my @files = @ARGV;

foreach my $file (@files) {
    open(my $fh, "<", $file)
        or die "Couldn't open $file: $!";
    my $base = basename($file, '_genomic.dgw');
    $base =~ s/(\.\d+)_.*/$1/;

    while (<$fh>) {
        if (/^>/) {
            $_ =~ s/^>/>${base}-/;
        }
        print $_;
    }
}

