#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     sort_uniq.pl - make sorted unique read data. 

e.g. perl sort_uniq.pl target
     perl sort_uniq.pl DRR054198
     perl sort_uniq.pl SRR8181712

     qsub -v target=ERR194147,tmpdir=/mnt/ssd sort_uniq.pl

     tmpdir can be ommitted.

Author: Akio Miyao <miyao@affrc.go.jp>

';

$uname = `uname`;
chomp($uname);
if ($uname eq "FreeBSD"){
    $sort_opt = "-S 100M";
}

if ($ARGV[0] ne ""){
    $target = $ARGV[0];
    $workdir = $target;
    $tmpdir = $workdir;
}elsif($ENV{target} ne ""){
    $target    = $ENV{target};
    $tmpdir    = $ENV{tmpdir};
    $cwd       = $ENV{PBS_O_WORKDIR};
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

open(IN, $cmd);
open(OUT, "|sort $sort_opt -T $tmpdir |uniq > $workdir/$target.sort_uniq");
while(<IN>){
    if ($count == 1 and !/N/){
	chomp;
	print OUT "$_\n";
	$complement = complement($_);
	print OUT "$complement\n";
    }elsif($count == 4){
	$count = 0;
    }
    $count++;
}
close(OUT);

if ($tmpdir ne $workdir){
    system("rm -r $tmpdir");
}

sub complement{
    my $seq = shift;
    my @seq = split('', $seq);
    my $i = length($seq);
    my $out = "";
    while($i > 0){
        $i--;
        if ($seq[$i] eq "A"){
            $out .= "T";
        }elsif($seq[$i] eq "C"){
            $out .= "G";
        }elsif($seq[$i] eq "G"){
            $out .= "C";
        }elsif($seq[$i] eq "T"){
            $out .= "A";
        }else{
            $out .= "N";
        }
    }
    return $out;
}
