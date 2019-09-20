#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     last_base_count.pl - convert to last base count from kmer count data. 

e.g. qsub -v target=ERR194147,tag=AAA,tmpdir=/mnt/ssd last_base_count.pl

     target is name of target.
     tmpdir can be ommited.

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] ne ""){
    print $usage;
    exit;
}

if($ENV{target} ne ""){
    $target    = $ENV{target};
    $tag       = $ENV{tag};
    $cwd       = $ENV{PBS_O_WORKDIR};
    $cwd       = $ENV{SGE_O_WORKDIR} if $ENV{SGE_O_WORKDIR} ne "";
    $tmpdir    = $ENV{tmpdir};
    $workdir = "$cwd/$target";
    &cluster;
}else{
    &standalone;
}

sub cluster{
    if ($tmpdir eq ""){
	$tmpdir = $workdir;
    }else{
	system("mkdir $tmpdir/$target");
	$tmpdir = "$tmpdir/$target";
	system("cp $workdir/$target.count.$tag.gz $tmpdir");
    }
    
    chdir $tmpdir;
    
    $output_file = "$target.lbc.$tag";
    
    open(IN, "/usr/bin/zcat $target.count.$tag.gz|");
    open(OUT, "> $output_file");
    while(<IN>){
	chomp;
	($seq, $count) = split;
	$tag = substr($seq, 0, 19);
	$nuc = substr($seq, 19, 1);
	
	if ($tag ne $prev and $prev ne ""){
	    print OUT "$prev\t$A\t$C\t$G\t$T\n";
	    $A = 0;
	    $C = 0;
	    $G = 0;
	    $T = 0;
	}
	$$nuc = $count;
	
	$prev = $tag;
    }
    print OUT "$prev\t$A\t$C\t$G\t$T\n";
    close(OUT);
    if ($tmpdir eq $workdir){
	system("gzip $output_file");
    }else{
	system("gzip $output_file && cp $output_file.gz $workdir && rm -r $tmpdir");
    }
}

sub standalone{
    while(<>){
	chomp;
	($seq, $count) = split;
	$tag = substr($seq, 0, 19);
	$nuc = substr($seq, 19, 1);
	
	if ($tag ne $prev and $prev ne ""){
	    print "$prev\t$A\t$C\t$G\t$T\n";
	    $A = 0;
	    $C = 0;
	    $G = 0;
	    $T = 0;
	}
	$$nuc = $count;
	
	$prev = $tag;
    }
    print "$prev\t$A\t$C\t$G\t$T\n";
}
