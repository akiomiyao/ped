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

e.g. perl demo_bidirectional.pl accession
     perl demo_bidirectional.pl SRR8181712

     This demonstration script is for sequence of Arabidopsis thaliana.

     Downloading time by fastq-dump from the Sequence Read Archive
     depends on the traffic from/to NCBI. Sometimes, it requires
     half day or more.

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] eq ""){
    print $usage;
    exit;
}

$target = $ARGV[0];

if (! -e "TAIR10/TAIR10.sort_uniq"){
    chdir "TAIR10";
    report("Making reference data.");
    system("/usr/bin/sh setup.sh");
    chdir "..";
}

if (! -e "$target/$target.sort_uniq"){
    report("Downloading fastq files of $target.");
    system("perl download.pl $target");
    report("Making sort_uniq sequence of $target.");
    system("perl sort_uniq.pl $target");
}

report("Aligning of $target sequence to reference genome.");
system("perl align.pl $target TAIR10 5");
report("Splitting alignment of SNP.");
system("perl split_snp.pl $target");
report("Splitting alignment of Indel.");
system("perl split_indel.pl $target");
opendir(DIR, $target);
while(readdir(DIR)){
    if (/$target.snp.[0-9][0-9]$/){
	$number = (split('\.', $_))[2];
	report("Verifying SNPs. $number");
	system("perl verify_snp.pl $target default TAIR10 $number bi");
    }
}
report("Verifying Indels.");
system("perl verify_indel.pl $target default TAIR10 01 bi");
report("Making vcf file of SNP");
system("perl snp2vcf.pl $target");
system("mv $target.indel.verify.01 $target.indel");
report("Complete.");

sub report{
    my $message = shift;
    my $now = `date`;
    chomp($now);
    print "$now $message\n";
}
