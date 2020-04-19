#!/usr/bin/perl
#
# This file is a script for Polymorphic Edge Detection.
#
# Multithreaded version of all in one.
#
# Copyright (C) 2019 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#
# Author: MIYAO Akio <miyao@affrc.go.jp>
#

$usage = '
 check_length.pl - output read counts for each sequence length.

 e.g. perl check_length.pl SRR11542243
 Fastq files should be saved in ./target/read directory. 

 Output data format: read lenth<TAB>number of reads

 If fastq data have multiple length reads, clipping sequence is required.
 For analisys sequence for RT-PCR product, check length using this script.
 
 For example, highest read count of fastq files is 101, 
 perl ped.pl target=SRR11542243,ref=COVID19,clipping=100 (or 99).

 Author: MIYAO Akio <miyao@affrc.go.jp>

';

$target = $ARGV[0];

die $usage if $target eq "";

opendir(DIR, "$target/read/");
foreach (readdir(DIR)){
    next if /^\./;
    if (/gz$/){
	open(IN, "zcat $target/read/$_ |");
    }elsif(/bz2$/){
	open(IN, "bzcat $target/read/$_ |");
    }elsif(/xz$/){
	open(IN, "xzcat $target/read/$_ |");
    }else{
	open(IN, "$target/read/$_");
    }
    while(<IN>){
	$line ++;
	if ($line % 4 == 2 and ! /N/){
	    chomp;
	    $len = length($_);
	    $length{$len}++; 
	}
    }
    close(IN);
}
close(DIR);


foreach (sort bynumber keys %length){
    print "$_\t$length{$_}\n";
}

sub bynumber{
    $a <=> $b;
}
