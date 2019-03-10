#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$target = $ARGV[0];
$type   = $ARGV[1];

if ($target eq ""){
    print "
e.g. perl snp2vcf.pl SRR5886855
     perl snp2vcf.pl SRR5886855 kmer

SRR5886855.snp.verify files in target directoy will be converted to SRR5886855.snp.vcf.

Default input file is result of bidirectional method.
For converting result of kmer method, option 'kmer' at 2nd argument is required.
" ;
    exit;
}

if ($type eq "kmer"){
    open(IN, "cat $target/$target.kmer.snp|");
    open(OUT, "> $target/$target.kmer.vcf");
}else{
    open(IN, "cat $target/$target.snp.verify.*|");
    open(OUT, "> $target/$target.snp.vcf");
}
print OUT "##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth)\">
##INFO=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$target\n";
while(<IN>){
    chomp;
    @row = split('\t', $_);
    $dp = $row[$#row] -1;
    if ($type eq "kmer"){
 	if (/M/){
	    print OUT "chr$row[0]\t$row[1]\t.\t$row[4]\t$row[5]\t1000\tPASS\tDP=$row[8];GT=1/1\tDP:GT\t$row[8]:1/1\n";
	}elsif(/H/){
	    print OUT "chr$row[0]\t$row[1]\t.\t$row[4]\t$row[5]\t1000\tPASS\tDP=$row[8];GT=0/1\tDP:GT\t$row[8]:1/0\n";
	}
   }else{
	if (/M/){
	    print OUT "chr$row[0]\t$row[1]\t.\t$row[2]\t$row[3]\t1000\tPASS\tDP=$row[8];GT=1/1\tDP:GT\t$row[8]:1/1\n";
	}elsif(/H/){
	    print OUT "chr$row[0]\t$row[1]\t.\t$row[2]\t$row[3]\t1000\tPASS\tDP=$row[8];GT=0/1\tDP:GT\t$row[8]:1/0\n";
	}
    }
}
