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
    open(IN, "cat $target/$target.bi.snp|");
    open(OUT, "> $target/$target.bi.vcf");
}
print OUT "##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth)\">
##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$target\n";
while(<IN>){
    chomp;
    @row = split('\t', $_);
    $dp = $row[$#row] -1;
    if (length($row[0]) <= 3 or $row[0] =~ /^[0-9]+$/){
	$row[0] = "chr" . $row[0];
    }
    if ($type eq "kmer"){
	$dp = $row[17] + $row[18];
	next if $dp == 0;
	$af = int(1000 * $row[18]/$dp)/1000;
 	if (/M/){
	    print OUT "$row[0]\t$row[1]\t.\t$row[4]\t$row[5]\t1000\tPASS\tGT=1/1;AF=$af;DP=$row[18]\tGT:AD:DP\t1/1:$row[17],$row[18]:$dp\n";
	}elsif(/H/){
	    print OUT "$row[0]\t$row[1]\t.\t$row[4]\t$row[5]\t1000\tPASS\tGT=0/1;AF=$af;DP=$row[18]\tGT:AD:DP\t0/1:$row[17],$row[18]:$dp\n";
	}
    }else{
	$dp = $row[7] + $row[8];
	next if $dp == 0;
	$af = int(1000 * $row[8]/$dp)/1000;
	if (/M/){
	    print OUT "$row[0]\t$row[1]\t.\t$row[2]\t$row[3]\t1000\tPASS\tGT=1/1;AF=$af;DP=$dp\tGT:AD:DP\t1/1:$row[7],$row[8]:$dp\n";
	}elsif(/H/){
	    print OUT "$row[0]\t$row[1]\t.\t$row[2]\t$row[3]\t1000\tPASS\tGT=0/1;AF=$af;DP=$dp\tGT:AD:DP\t0/1:$row[7],$row[8]:$dp\n";
	}
    }
}
