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
     perl split_snp.pl target
     perl split_snp ERR194147

     qsub -v target=ERR194147,tmpdir=/mnt/ssd split_snp.pl

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
    $tmpdir = $ARGV[1];
    $cwd    = `pwd`;
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

if ($tmpdir ne ""){
    $tmpdir = $tmpdir . "/pedtmp." . int(rand 1000000);
    system("mkdir $tmpdir");
}else{
    $tmpdir = "$cwd/$target";
}

opendir(DIR, "$cwd/$target/");
foreach (readdir(DIR)){
    if (/$target.aln/){
	$aln_count ++;
    }
}
closedir(DIR);

open(IN, "cat $cwd/$target/$target.aln.*|");
open(OUT, "|sort $sort_opt -T $tmpdir -k 1 -k 2 -n |uniq -c > $tmpdir/$target.snp.tmp");
while(<IN>){
    if (/snp/){
	chomp;
	@row = split;
	$chr = "00$row[1]";
	$chr = substr($chr, length($chr) - 3, 3);
	print OUT "$chr\t$row[2]\t$row[4]\t$row[5]\n";
    }
}
close(IN);
close(OUT);

while(1){
    $mtime = (stat("$tmpdir/$target.snp.tmp"))[9];
    if (time > $mtime + 10){
	last;
    }
    sleep 1;
}

$fcount = "01";
open(IN, "$tmpdir/$target.snp.tmp");
open(OUT, "> $cwd/$target/$target.snp.$fcount");
while(<IN>){
    if ($count == 3000000){
        $fcount ++;
        open(OUT, "> $cwd/$target/$target.snp.$fcount");
        $count = 0;
    }
    $count++;
    chomp;
    @row = split;
    $row[0] = int($row[0] / $aln_count) + 1;
    print OUT "$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[0]\n";
}
close(IN);
close(OUT);

system("rm $tmpdir/$target.snp.tmp");

if ($tmpdir ne "$cwd/$target"){
    system("rm -r $tmpdir");
}
