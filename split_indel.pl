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

e.g. perl split_indel.pl target reference
     perl split_indel ERR194147 hg38

     qsub -v target=ERR194147,ref=hg38 split_indel.pl


Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] ne ""){
    $target  = $ARGV[0];
    $ref     = $ARGV[1];
    $cwd = `pwd`;
    chomp($cwd);
}elsif($ENV{target} ne ""){
    $target    = $ENV{target};
    $ref       = $ENV{ref};
    $cwd       = $ENV{PBS_O_WORKDIR};
}else{
    print $usage;
    exit;
}

$workdir = "$cwd/$target";

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

open(IN, "cat $workdir/$target.aln.*|");
open(OUT, "|sort -T $workdir |uniq > $workdir/$target.indel.tmp");
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

system("rm $workdir/$target.indel.sort") if -e "$workdir/$target.indel.sort";
foreach $chr (@chr){
    print ERR "$chr\n";
    open(IN, "$workdir/$target.indel.tmp");
    open(OUT, "|sort -k 2 -n -T $workdir >> $workdir/$target.indel.sort");
    while(<IN>){
	@row = split;
	if ($row[1] eq $chr){
	    print OUT;
	}
    }
    close(IN);
    close(OUT);
}
system("rm $workdir/$target.indel.tmp");

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
