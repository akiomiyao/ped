#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2019 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

if ($ENV{PBS_O_WORKDIR} ne ""){
    $cwd = $ENV{PBS_O_WORKDIR};
    chdir $cwd;
    require "$cwd/common.pl";
}elsif($ENV{SGE_O_WORKDIR} ne ""){
    $cwd = $ENV{SGE_O_WORKDIR};
    chdir $cwd;
    require "$cwd/common.pl";
}else{
    require './common.pl';
}

$usage = '
     index.pl - indexing program for bidirectional alignments. 

     For example,
     perl index.pl target

     qsub -v target=ERR194147 index.pl

     Aftger the indexing, alignments can be found by search.pl.

     perl search.pl target chromosome position

     perl search,pl ERR194146 1 2915759

     Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] ne ""){
    $target = $ARGV[0];
}elsif($ENV{target} ne ""){
    $target    = $ENV{target};
}else{
    print $usage;
    exit;
}

opendir(DIR, $target);
foreach (readdir(DIR)){
    if (/$target.aln/){
	chomp;
	@row = split('\.', $_);
	$margin{$row[2]} = 1;
    }
}

open(OUT, "|sort $sort_opt -T . > $target/$target.index");
foreach  $margin (sort keys %margin){
    $pos  = 0;
    open(IN, "$target/$target.aln.$margin");
    while(<IN>){
	$length = length($_);
	if(/^#/){
	    @row = split;
	    if ($row[1] =~ /^[0-9]*$/){
		$row[1] = "000$row[1]";
		$row[1] = substr($row[1], length($row[1]) -3, 3);
	    }
	    $row[2] = "0000000000$row[2]";
	    $row[2] = substr($row[2], length($row[2]) - 10, 10);
	    print OUT "$row[1] $row[2] $margin $pos\n";
	}
	$pos += $length;
    }
    close(IN);
}
close(OUT);
