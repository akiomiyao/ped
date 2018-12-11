#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     count.pl - create count of kmers from sort_uniq data. 

e.g. qsub -v target=ERR194147,number=0001,tmpdir=/mnt/ssd count.pl

     target is name of target.
     number is number of split target file (e.g. 0001).
     tmpdir can be ommited.

Author: Akio Miyao <miyao@affrc.go.jp>

';

if($ENV{target} ne ""){
    $target    = $ENV{target};
    $number    = $ENV{number};
    $cwd       = $ENV{PBS_O_WORKDIR};
    $tmpdir    = $ENV{tmpdir};
    $workdir = "$cwd/$target";
}else{
    print $usage;
    exit;
}

if ($tmpdir eq ""){
    $tmpdir = $workdir;
}else{
    if (-e "$tmpdir/$target"){
	system("rm -r $tmpdir/$target");
    }
    system("mkdir $tmpdir/$target");
    $tmpdir = "$tmpdir/$target";
}

chdir $tmpdir;

$input_file = "$target.sort_uniq.$number";

if ($tmpdir ne $workdir){
    system("cp $workdir/tmp/$input_file $tmpdir");
}
open(IN, "$input_file");
open(OUT, "|sort -S 2G -T $tmpdir | uniq -c | awk '{print \$2 \"\t\" \$1}' | /usr/bin/perl $cwd/split_count.pl");
while(<IN>){
    chomp;
    $length = length($_);
    for ($i = 0; $i <= $length - 20; $i++){
	$seq = substr($_, $i, 20);
	if ($seq !~ /N/){
	    print OUT "$seq\n";
	}
    }
}
