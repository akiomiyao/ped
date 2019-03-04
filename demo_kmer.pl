#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     demo_kmer.pl - demonstration of kmer method. 

e.g. perl demo_kmer.pl accession
     perl demo_kmer.pl SRR8181712
     
     This demonstration script is for sequence of Arabidopsis thaliana.

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] eq ""){
    print $usage;
    exit;
}

$target = $ARGV[0];

report("demo_kmer.pl start.");

if (! -e "TAIR10/TAIR10.sort_uniq"){
    chdir "TAIR10";
    report("Making reference data.");
    system("/usr/bin/sh setup.sh");
    chdir "..";
}

if (! -e "TAIR10/TAIR10.lbc.AAA.gz"){
    report("Making kmer count of TAIR10.");
    system("perl count.pl TAIR10");
    system("gzip TAIR10/TAIR10.lbc.*");
}

if (! -e "$target/$target.sort_uniq"){
    report("Downloading fastq files of $target.");
    system("perl download.pl $target");
    report("Making sort_uniq sequence of $target.");
    system("perl sort_uniq.pl $target");
}

report("Making kmer count of $target.");
system("perl count.pl $target");
report("Making snp");
system("perl snp.pl $target TAIR10");
report("Making map");
system("perl map.pl $target TAIR10");
report("Verifying");
system("perl verify_snp.pl $target default TAIR10 AAA kmer");
report("Convert to vcf");
system("perl snp2vcf.pl $target kmer");
report("demo_kmer.pl done.");


sub report{
    my $message = shift;
    my $now = `date`;
    chomp($now);
    print "$now $message\n";
}
