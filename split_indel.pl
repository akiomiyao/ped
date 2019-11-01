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

     For example,
     perl split_indel.pl target reference
     perl split_indel ERR194147 hg38

     qsub -v target=ERR194147,ref=hg38,tmpdir=/mnt/ssd split_indel.pl
 
     tmpdir is optional.
     tmpdir should be specified to a fast local disk (e.g. SSD) in each node.
     If tmpdir is not specified, target dir will be used for tmpdir.

     Author: Akio Miyao <miyao@affrc.go.jp>

';

$uname = `uname`;
chomp($uname);
if ($uname eq "FreeBSD"){
    $sort_opt = "-S 100M";
}

if ($ARGV[0] ne ""){
    $target  = $ARGV[0];
    $tmpdir  = $ARGV[1];
    $cwd = `pwd`;
    chomp($cwd);
}elsif($ENV{target} ne ""){
    $target    = $ENV{target};
    $tmpdir    = $ENV{tmpdir};
    $cwd       = $ENV{PBS_O_WORKDIR};
    $cwd       = $ENV{SGE_O_WORKDIR} if $ENV{SGE_O_WORKDIR} ne "";
}else{
    print $usage;
    exit;
}

$workdir = "$cwd/$target";

if ($tmpdir ne ""){
    $tmpdir = $tmpdir . "/pedtmp." . int(rand 1000000);
    system("mkdir $tmpdir");
}else{
    $tmpdir = "$cwd/$target";
}

open(IN, "cat $workdir/$target.aln.*|");
open(OUT, "|sort $sort_opt -T $tmpdir -k 1 -k 2 -n |uniq > $tmpdir/$target.indel.sort");
while(<IN>){
    chomp;
    if (/ion/){
	$flag = 1;
	@row = split;
	$chr = "00$row[1]";
	$chr = substr($chr, length($chr) - 3, 3);
	if ($row[7] eq ""){
	    $out .= "$row[0] $chr $row[2] $row[3] $row[4] $row[5] $row[6]\t";
	}else{
	    $out .= "$row[0] $chr $row[2] $row[3] $row[4] $row[5] $row[6] $row[7]\t";
	}
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

while(1){
    $mtime = (stat("$tmpdir/$target.indel.sort"))[9];
    if (time > $mtime + 10){
	last;
    }
    sleep 1;
}

$fcount = "01";
open(IN, "$tmpdir/$target.indel.sort");
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

if ($tmpdir ne $workdir){
    system("mv $tmpdir/$target.indel.sort $workdir");
    system("rm -r $tmpdir");
}
