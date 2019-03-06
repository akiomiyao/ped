#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     kmer.pl - pipeline of kmer method for standalone machine. 

     run  perl kmer.pl accession  reference
     e.g. perl kmer.pl SRR8181712 TAIR10

     This script is for stand alone machine.

     Before run the script, downloading target reads and reference genome data
     are required.

     To make the reference data,
     run  perl mkref.pl reference
     e.g. perl mkref.pl TAIR10

     Run without reference name, references in config file will be listed.
     You can select reference from the list.

     To download from SRA in NCBI,
     run  perl download.pl accession
     e.g. perl download.pl SRR8181712

     To analyze your data of short reads,
     mkdir mydata
     mkdir mydata/read
     and the copy your fastq data to mydata/read
     and then run perl kmer.pl mydata reference

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[1] eq ""){
    print $usage;
    exit;
}

$target = $ARGV[0];
$ref    = $ARGV[1];

report("kmer.pl start.");

if (! -e "$ref/$ref.lbc.AAA.gz"){
    report("Making kmer count of $ref.");
    system("perl count.pl $ref");
    system("gzip $ref/$ref.lbc.*");
}

report("Making kmer count of $target.");
system("perl count.pl $target");
report("Making snp");
system("perl snp.pl $target $ref");
report("Making map");
system("perl map.pl $target $ref");
report("Verifying");
system("perl verify_snp.pl $target default $ref AAA kmer");
report("Convert to vcf");
system("perl snp2vcf.pl $target kmer");
report("kmer.pl done.");


sub report{
    my $message = shift;
    my $now = `date`;
    chomp($now);
    print "$now $message\n";
}
