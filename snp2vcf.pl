#!/usr/local/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$target = $ARGV[0];

if ($target eq ""){
    print "e.g. perl snp2vcf.pl SRR5886855
SRR5886855.snp.verify files in target directoy will be converted to SRR5886855.snp.vcf.
" ;
    exit;
}

open(IN, "cat $target/$target.snp.verify.*|");
open(OUT, "> $target/$target.snp.vcf");

while(<IN>){
    chomp;
    @row = split;
    if (/M/){
	print OUT "chr$row[0]\t$row[1]\t.\t$row[2]\t$row[3]\t0\tPASS\tDP=$row[8];GT=1/1\n";
    }elsif(/H/){
	print OUT "chr$row[0]\t$row[1]\t.\t$row[2]\t$row[3]\t0\tPASS\tDP=$row[8];GT=0/1\n";
    }
}
