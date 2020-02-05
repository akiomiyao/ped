#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

require "./common.pl";

$usage = '
     sort_uniq.pl - make sorted unique read data. 

e.g. perl sort_uniq.pl target
     perl sort_uniq.pl DRR054198
     perl sort_uniq.pl SRR8181712

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] ne ""){
    $target = $ARGV[0];
}else{
    print $usage;
    exit;
}

&report("sort_uniq.pl: execute. Making $target.sort_uniq");

opendir(DIR, "$cwd/$target/read");
foreach(sort grep(! /^\.|download.sh/, readdir(DIR))){
    if (/fq$|fastq$/){
	if ($type ne "" and $type ne "fastq"){
	    print "Different compression format in read directory is not acceptable.";
	    exit;
	}
	$cmd = "cat $target/read/*|";
	$type = "fastq";
    }elsif(/gz$/){
	if ($type ne "" and $type ne "gz"){
	    print "Different compression format in read directory is not acceptable.";
	    exit;
	}
	$cmd = "$zcat $target/read/* 2> /dev/null |";
	$type = "gz";
    }elsif(/bz2$/){
	if ($type ne "" and $type ne "bz2"){
	    print "Different compression format in read directory is not acceptable.";
	    exit;
	}
	$cmd = "bzcat $target/read/* 2> /dev/null |";
	$type = "bz2";
    }elsif(/xz$/){
	if ($type ne "" and $type ne "xz"){
	    print "Different compression format in read directory is not acceptable.";
	    exit;
ll
	}
	$cmd = "xzcat $target/read/* 2> /dev/null |";
	$type = "xz";
    }
}
closedir(DIR);

if (! -e "$cwd/$target/sort_uniq"){
    system("mkdir $cwd/$target/sort_uniq");
}

&report("Making $target.sort_uniq: Split to subfiles");

foreach $nuca (@nuc){
    foreach $nucb (@nuc){
	foreach $nucc (@nuc){
	    $tag = $nuca . $nucb . $nucc;
	    open($tag, "> $cwd/$target/sort_uniq/$tag.seq")
	}
    }
}

open(IN, $cmd);
while(<IN>){
    if ($count == 1 and !/N/){
	chomp;
	$tag = substr($_, 0, 3);
	print $tag "$_\n";
	$complement = &complement($_);
	$tag = substr($complement, 0, 3);
	print $tag "$complement\n";
	$total ++;
	if ($total % 1000000 == 0){
	    &report("Making $target.sort_uniq files: Split to subfiles. $total reads processed");
	}
    }elsif($count == 4){
	$count = 0;
    }
    $count++;
}
close(OUT);

&closeTag;
&report("Making $target.sort_uniq files: Split to subfiles. Done. $total reads processed");

foreach $nuca (@nuc){
    foreach $nucb (@nuc){
	foreach $nucc (@nuc){
	    $tag = $nuca . $nucb . $nucc;
	    &report("Making $target.sort_uniq files: Sorting for $target.sort_uniq.$tag.gz");
	    system("sort -T $cwd/$target/sort_uniq/ $sort_opt $cwd/$target/sort_uniq/$tag.seq | uniq | gzip > $cwd/$target/sort_uniq/$target.sort_uniq.$tag.gz && rm $cwd/$target/sort_uniq/$tag.seq");
	}
    }
}

&report("sort_uniq.pl for $target: complete");
