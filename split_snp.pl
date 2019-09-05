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

     For example,
     perl split_snp.pl target ref
     perl split_snp ERR194147 hg38

     qsub -v target=ERR194147,ref=hg38,tmpdir=/mnt/ssd split_snp.pl

     tmpdir is fast local disk, e.g. SSD, in nodes.
     If tmpdir is not specified, target dir will be used for tmpdir.

     Author: Akio Miyao <miyao@affrc.go.jp>

';

$uname = `uname`;
chomp($uname);
if ($uname eq "FreeBSD"){
    $sort_opt = "-S 100M";
}

if ($ARGV[0] ne ""){
    $target = $ARGV[0];
    $ref    = $ARGV[1];
    $tmpdir = $ARGV[2];
    $cwd    = `pwd`;
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
open(OUT, "|sort $sort_opt -T $tmpdir |uniq > $tmpdir/$target.aln.sort");
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

while(1){
    $mtime = (stat("$tmpdir/$target.aln.sort"))[9];
    if (time > $mtime + 10){
	last;
    }
    sleep 1;
}

open(IN, "$tmpdir/$target.aln.sort");
open(OUT, "|sort $sort_opt -T $tmpdir |uniq -c |awk '{print \$2, \$3, \$4, \$5, \$1}' >$tmpdir/$target.snp.tmp");
while(<IN>){
    foreach $dat(split('\t', $_)){
	if ($dat =~/^#/){
	    @row = split(' ', $dat);
	    print OUT "$row[1] $row[2] $row[4] $row[5]\n";
	}
    }
}
close(IN);
close(OUT);

system("mkdir $tmpdir/snp.tmp");
open(IN, "$tmpdir/$target.snp.tmp");
while(<IN>){
    @row = split;
    if ($prev ne $row[0]){
	open(OUT, "> $tmpdir/snp.tmp/$row[0]");
	$detected{$row[0]} = 1;
    }
    print OUT;
    $prev = $row[0];
}
close(IN);
close(OUT);

system("rm $workdir/$target.snp") if -e "$workdir/$target.snp";
foreach $chr (@chr){
    next if $detected{$chr} != 1 or $chr eq "NOP";
    system("sort -k 2 -n $sort_opt -T $tmpdir $tmpdir/snp.tmp/$chr>> $tmpdir/$target.snp");
}
system("rm -r $tmpdir/snp.tmp/ $tmpdir/$target.snp.tmp");

while(1){
    $mtime = (stat("$tmpdir/$target.snp"))[9];
    if (time > $mtime + 10){
	last;
    }
    sleep 1;
}

$fcount = "01";
open(IN, "$tmpdir/$target.snp");
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

if ($tmpdir ne $workdir){
    system("mv $tmpdir/$target.aln.sort $workdir");
    system("rm -r $tmpdir");
}
