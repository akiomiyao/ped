#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     demo_bidirectional.pl - demonstration of bidirectional alignment method. 

e.g. perl demo_bidirectional.pl

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] ne ""){
    print $usage;
    exit;
}

chdir "TAIR10";
report("Making reference data.");
system("/usr/bin/sh setup.sh");

chdir "..";

report("Downloading fastq files of SRR8181712.");
system("perl download.pl SRR8181712");
report("Making sort_uniq sequence of SRR8181712.");
system("perl sort_uniq.pl SRR8181712");
report("Aligning of SRR8181712 sequence to reference genome.");
system("perl align.pl SRR8181712 TAIR10 5");
report("Splitting alignment of SNP.");
system("perl split_snp.pl SRR8181712");
report("Splitting alignment of Indel.");
system("perl split_indel.pl SRR8181712");
report("Verifying SNPs.");
system("perl verify_snp.pl SRR8181712 default TAIR10 01 bi");
report("Verifying Indels.");
system("perl verify_indel.pl SRR8181712 default TAIR10 01 bi");
report("Making vcf file of SNP");
system("perl snp2vcf.pl SRR8181712");
system("mv SRR8181712.indel.verify.01 SRR8181712.indel");
report("Complete.");

sub report{
    my $message = shift;
    my $now = `date`;
    chomp($now);
    print "$now $message\n";
}
