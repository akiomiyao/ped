#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     split_sort_uniq.pl - split sort_uniq file for kmer method. 

e.g. perl split_sort_uniq.pl target

     qsub -v target=ERR194147 split_sort_uniq.pl

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] ne ""){
    $target = $ARGV[0];
    $cwd = `pwd`;
    chomp($cwd);
    $workdir = "$cwd/$target";
}elsif($ENV{target} ne ""){
    $target    = $ENV{target};
    $cwd       = $ENV{PBS_O_WORKDIR};
    $cwd       = $ENV{SGE_O_WORKDIR} if $ENV{SGE_O_WORKDIR} ne "";
    $workdir = "$cwd/$target";
}else{
    print $usage;
    exit;
}

$file_count = "0001";
chdir $workdir;
system("mkdir tmp") if ! -d "tmp";

if (-e "$target.sort_uniq"){
    open(IN, "$target.sort_uniq");
}elsif (-e "$target.sort_uniq.gz"){
    open(IN, "zcat $target.sort_uniq.gz|");
}
open(OUT, "|gzip > tmp/$target.sort_uniq.$file_count.gz");
while(<IN>){
    if ($count == 10000000){
	$file_count ++;
	open(OUT, "| gzip > tmp/$target.sort_uniq.$file_count.gz");
	$count = 0;
    }
    $count++;
    print OUT;
}
close(IN);
close(OUT);
