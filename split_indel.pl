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
    $ref     = $ARGV[1];
    $tmpdir  = $ARGV[2];
    $cwd = `pwd`;
    chomp($cwd);
}elsif($ENV{target} ne ""){
    $target    = $ENV{target};
    $ref       = $ENV{ref};
    $tmpdir    = $ENV{tmpdir};
    $cwd       = $ENV{PBS_O_WORKDIR};
}else{
    print $usage;
    exit;
}

$workdir = "$cwd/$target";

if ($tmpdir ne ""){
    $tmpdir = $tmpdir . "/$target";
    system("mkdir $tmpdir");
}else{
    $tmpdir = $workdir;
}

open(IN, "$cwd/config");
while(<IN>){
    chomp;
    @row = split;
    if($row[0] eq $ref && $row[1] eq "chromosome"){
	if ($row[3] != 0){
	    for ($i = $row[2]; $i <= $row[3]; $i++){
		push(@chr, $i);
	    }
	}
	if ($row[4] ne ""){
	    foreach ($i = 4; $i <= $#row; $i++){
		push(@chr, $row[$i]);
	    }
	}
    }
}
close(IN);

if ($chr[0] eq ""){
    opendir(REF, "$cwd/$ref");
    foreach(sort readdir(REF)){
	chomp;
	if (/^chr/){
	    ($chr = $_) =~ s/^chr//; 
	    push(@chr, $chr);
	}
    }
    closedir(REF);
}

open(IN, "cat $workdir/$target.aln.*|");
open(OUT, "|sort $sort_opt -T $tmpdir |uniq > $tmpdir/$target.indel.tmp");
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

system("mkdir $tmpdir/indel.tmp");
open(IN, "$tmpdir/$target.indel.tmp");
while(<IN>){
    @row = split;
    if ($prev ne $row[1]){
	open(OUT, "> $tmpdir/indel.tmp/$row[1]");
    }
    print OUT;
    $prev = $row[1];
    $detected{$row[1]} =1
}
close(IN);
close(OUT);

system("rm $workdir/$target.indel.sort") if -e "$workdir/$target.indel.sort";
foreach $chr (sort keys %detected){
    system("sort $sort_opt -k 2 -n -T $tmpdir $tmpdir/indel.tmp/$chr >> $tmpdir/$target.indel.sort");
}
system("rm -r $tmpdir/indel.tmp $tmpdir/$target.indel.tmp");

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
