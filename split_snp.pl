#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     split_snp.pl - split data of align.pl output. 

e.g. perl split_snp.pl target ref
     perl split_snp ERR194147 hg38

     qsub -v target=ERR194147,ref=hg38,tmpdir=/mnt/ssd split_snp.pl

     tmpdir can be ommitted.

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] ne ""){
    $target = $ARGV[0];
    $ref = $ARGV[1];
    $cwd = `pwd`;
    chomp($cwd);
    $workdir = "$cwd/$target";
}elsif($ENV{target} ne ""){
    $target    = $ENV{target};
    $ref       = $ENV{ref};
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

open(IN, "$cwd/config");
while(<IN>){
    chomp;
    @row = split;
    if($row[0] eq $ref && $row[1] eq "chromosome"){
	if ($row[2] != 0){
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

if ($cwd eq ""){
    open(IN, "cat $target.aln.*|");
}else{
    open(IN, "cat $cwd/$target/$target.aln.*|");
}
open(OUT, "|sort -T $workdir |uniq > $workdir/$target.aln.sort");
while(<IN>){
    chomp;
    if (/snp/){
	$flag = 1;
	push(@dat, $_);
    }elsif($flag){
	if ($_ eq ""){
	    chomp($out);
	    foreach $dat (@dat){
		print OUT "$dat\t$out\n";
	    }
	    $flag = 0;
	    $out = "";
	    @dat = ();
	}else{
	    $out .= "$_\t";
	}
    }
}
close(IN);
close(OUT);

open(IN, "$workdir/$target.aln.sort");
open(OUT, "|sort -T $workdir |uniq -c |awk '{print \$2, \$3, \$4, \$5, \$1}' >$workdir/$target.snp.tmp");
while(<IN>){
    foreach $dat(split('\t', $_)){
	if ($dat =~/^#/){
	    @row = split(' ', $dat);
	    print OUT "$row[1] $row[2] $row[4] $row[5]\n";
	}
    }
}

system("rm $workdir/$target.snp") if -e "$workdir/$target.snp";
foreach $chr (@chr){
    open(IN, "$workdir/$target.snp.tmp");
    open(OUT, "|sort -k 2 -n -T $workdir >> $workdir/$target.snp");
    while(<IN>){
	@row = split;
	if ($row[0] eq $chr){
	    print OUT;
	}
    }
    close(IN);
    close(OUT);
}
system("rm $workdir/$target.snp.tmp");

$fcount = "01";
open(IN, "$workdir/$target.snp");
open(OUT, "> $workdir/$target.snp.$fcount");
while(<IN>){
    if ($count == 3000000){
        $fcount ++;
        open(OUT, "> $workdir/$target.snp.$fcount");
        $count = 0;
    }
    $count++;
    print OUT;
}
close(IN);
close(OUT);

if ($tmpdir ne ""){
    system("cp $target.aln.sort $target.snp $target.snp.?? $cwd/$target && rm -r $workdir");
}
