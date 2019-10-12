#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

if ($ENV{PBS_O_WORKDIR} ne ""){
    $cwd = $ENV{PBS_O_WORKDIR};
    chdir $cwd;
    require "$cwd/common.pl";
}elsif($ENV{SGE_O_WORKDIR} ne ""){
    $cwd = $ENV{SGE_O_WORKDIR};
    chdir $cwd;
    require "$cwd/common.pl";
}else{
    $cwd = `pwd`;
    chomp($cwd);
    require './common.pl';
}

$usage = '
     qsub_sort_uniq.pl - make sorted unique read data. 

e.g. perl sort_uniq.pl target
     perl sort_uniq.pl DRR054198
     perl sort_uniq.pl SRR8181712

     qsub -v target=ERR194147,tmpdir=/mnt/ssd sort_uniq.pl

     tmpdir can be ommitted.

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] ne ""){
    $target = $ARGV[0];
    $workdir = $target;
    $tmpdir = $workdir;
}elsif($ENV{target} ne ""){
    $target    = $ENV{target};
    $tmpdir    = $ENV{tmpdir};
    $cwd       = $ENV{PBS_O_WORKDIR};
    $cwd       = $ENV{SGE_O_WORKDIR} if $ENV{SGE_O_WORKDIR} ne "";
    $workdir = "$cwd/$target";
    if ($tmpdir eq ""){
	$tmpdir = $workdir;
    }else{
	system("mkdir $tmpdir/$target");
	$tmpdir = "$tmpdir/$target";
    }
}else{
    print $usage;
    exit;
}

&report("Job start making $target.sort_uniq");

opendir(DIR, "$workdir/read");
foreach(sort grep(! /^\.|download.sh/, readdir(DIR))){
    if (/fq$|fastq$/){
	if ($type ne "" and $type ne "fastq"){
	    print "Different compression format in read directory is not acceptable.";
	    exit;
	}
	$cmd = "cat $workdir/read/*|";
	$type = "fastq";
    }elsif(/gz$/){
	if ($type ne "" and $type ne "gz"){
	    print "Different compression format in read directory is not acceptable.";
	    exit;
	}
	$cmd = "zcat $workdir/read/* 2> /dev/null |";
	$type = "gz";
    }elsif(/bz2$/){
	if ($type ne "" and $type ne "bz2"){
	    print "Different compression format in read directory is not acceptable.";
	    exit;
	}
	$cmd = "bzcat $workdir/read/* 2> /dev/null |";
	$type = "bz2";
    }elsif(/xz$/){
	if ($type ne "" and $type ne "xz"){
	    print "Different compression format in read directory is not acceptable.";
	    exit;
	}
	$cmd = "xzcat $workdir/read/* 2> /dev/null |";
	$type = "xz";
    }
}
closedir(DIR);

foreach $nuca (@nuc){
    foreach $nucb (@nuc){
	foreach $nucc (@nuc){
	    $tag = $nuca . $nucb . $nucc;
	    open($tag, "> $workdir/sort_uniq/$tag.seq")
	}
    }
}

open(IN, $cmd);
while(<IN>){
    if ($count == 1 and !/N/){
	chomp;
	$tag = substr($_, 0, 3);
	print $tag "$_\n";
	$complement = complement($_);
	$tag = substr($complement, 0, 3);
	print $tag "$complement\n";
    }elsif($count == 4){
	$count = 0;
    }
    $count++;
}
close(OUT);

&closeTag;

chdir $cwd;

foreach $nuca (@nuc){
    foreach $nucb (@nuc){
	foreach $nucc (@nuc){
	    $tag = $nuca . $nucb . $nucc;
	    $qsub = "-v target=$target,tag=$tag sort_uniq_sub.pl";
	    &doQsub($qsub);
	}
    }
}
