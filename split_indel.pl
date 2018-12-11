#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     split_indel.pl - split data of align.pl output. 

e.g. perl split_indel.pl target
     perl split_indel ERR194147

     qsub -v target=ERR194147,tmpdir=/mnt/ssd split_indel.pl

     tmpdir can be ommitted.

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] ne ""){
    $target = $ARGV[0];
    $workdir = `pwd`;
    chomp($workdir);
    $workdir .= "/$target";
}elsif($ENV{target} ne ""){
    $target    = $ENV{target};
    $tmpdir    = $ENV{tmpdir};
    $cwd       = $ENV{PBS_O_WORKDIR};
    if ($tmpdir ne ""){
	if (-d "$tmpdir/$target"){
	    system("rm -r $tmpdir/$target");
	}
	system("mkdir $tmpdir/$target");
	$workdir = "$tmpdir/$target";
    }else{
	$workdir = "$cwd/$target";
    }
}else{
    print $usage;
    exit;
}

chdir $workdir;

if ($cwd eq ""){
    open(IN, "cat $workdir/$target.aln.*|");
}else{
    open(IN, "cat $cwd/$target/$target.aln.*|");
}
open(OUT, "|sort -S 100M -T $workdir |uniq > $workdir/$target.indel.sort");
while(<IN>){
    chomp;
    if (/ion/){
	$flag = 1;
	$out .= "$_\t";
	$count = 0;
    }elsif($flag and $count != 1){
	if ($_ eq ""){
	    print OUT "$out\n";
	    $flag = 0;
	    $out = "";
	}else{
	    $out .= "$_\t";
	}
    }
    $count++;
}
close(IN);
close(OUT);

$fcount = "01";
open(IN, "$workdir/$target.indel.sort");
open(OUT, "> $workdir/$target.indel.$fcount");
while(<IN>){
    if ($count == 10000000){
        $fcount ++;
	print "$fcount\n";
        open(OUT, "> $workdir/$target.indel.$fcount");
        $count = 0;
    }
    $count++;
    print OUT;
}
close(IN);
close(OUT);

if ($tmpdir ne ""){
    system("cp $target.indel.sort $target.indel.?? $cwd/$target && rm -r $workdir");
}
