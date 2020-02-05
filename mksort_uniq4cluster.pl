#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2019 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$target    = $ENV{target};
$execution = $ENV{execution};
$number    = $ENV{number};
$tag       = $ENV{tag};
$tmpdir    = $ENV{tmpdir};

$usage = '
     mksort_uniq4cluster.pl - make sorted unique read data using computer cluster. 

e.g. perl mksort_uniq4cluster.pl target tmpdir
     perl mksort_uniq4cluster.pl DRR054198 /mnt/ssd

     Specify of tmpdir in each node is required.

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ENV{PBS_O_WORKDIR} ne ""){
    $cwd = $ENV{PBS_O_WORKDIR};
    chdir $cwd;
    require "$cwd/common.pl";
}elsif($ENV{SGE_O_WORKDIR} ne ""){
    $cwd = $ENV{SGE_O_WORKDIR};
    chdir $cwd;
    require "$cwd/common.pl";
}elsif($ARGV[1] ne ""){
    $target = $ARGV[0];
    $tmpdir = $ARGV[1];
    require "./common.pl";
}else{
    print $usage;
    exit;
}

if ($execution eq ""){
    &report("Job Start");
    &report("Split read");
    $qsub = "qsub -v target=$target,execution=splitRead,tmpdir=$tmpdir mksort_uniq4cluster.pl";
    system($qsub);
    while(1){
	last if -e "$cwd/$target/sort_uniq/splitRead.done";
	sleep 10;
    }
    &report("Split with tag");
    open(IN, "$cwd/$target/sort_uniq/splitRead.done");
    $fcount = <IN>;
    close(IN);
    for($i = 1; $i <= $fcount; $i++){
	$qsub = "qsub -v target=$target,execution=splitTag,number=$i,tmpdir=$tmpdir mksort_uniq4cluster.pl";
	print "$qsub\n";
	system($qsub);
	sleep 1;
    }
    while(1){
	$count = 0;
	opendir(DIR, "$cwd/$target/sort_uniq");
	foreach (sort readdir(DIR)){
	    if (/$target.read/){
		$count++;
	    }
	}
	sleep 10;
	last if ($count == 0);
    }
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		system("qsub -v target=$target,execution=sort_uniq,tag=$tag,tmpdir=$tmpdir mksort_uniq4cluster.pl");
		sleep 1;
	    }
	}
    }
}elsif($execution eq "splitRead"){
    &splitRead;
}elsif($execution eq "splitTag"){
    &splitTag;
}elsif($execution eq "sort_uniq"){
    &sort_uniq;
}

sub sort_uniq{
    if ($tmpdir ne ""){
	$tmpdir = "$tmpdir/$target";
	system("mkdir $tmpdir");
    }
    open(IN, "$cwd/$target/sort_uniq/splitRead.done");
    $fcount = <IN>;
    close(IN);
    return if $fcount == 0;
    while(1){
	$count = 0;
	opendir(DIR, "$cwd/$target/sort_uniq/");
	foreach(sort readdir(DIR)){
	    if (/$tag.seq/){
		$mtime = (stat($_))[9];
		if (time > $mtime + 5){
		    $count++;
		}
	    }
	}
	close(DIR);
	last if $count == $fcount;
    }
    system("cp $cwd/$target/sort_uniq/$tag.seq.* $tmpdir");
    system("cat $tmpdir/$tag.seq.* |sort $sort_opt -T $tmpdir | uniq |gzip > $tmpdir/$target.sort_uniq.$tag.gz");
    &waitFile("$tmpdir/$target.sort_uniq.$tag.gz");
    system("cp $tmpdir/$target.sort_uniq.$tag.gz $cwd/$target/sort_uniq/ && rm -r $tmpdir && rm $cwd/$target/sort_uniq/$tag.seq.*");
}

sub splitTag{
    if ($tmpdir ne ""){
	$tmpdir = "$tmpdir/$target";
	system("rm -r $tmpdir") if -e $tmpdir;
	system("mkdir $tmpdir");
    }
    system("cp $cwd/$target/sort_uniq/$target.read.$number $tmpdir");
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "> $tmpdir/$tag.seq.$number")
	    }
	}
    }
    open(IN, "$tmpdir/$target.read.$number");
    while(<IN>){
	chomp;
	$tag = substr($_, 0, 3);
	print $tag "$_\n";
	$complement = complement($_);
	$tag = substr($complement, 0, 3);
	print $tag "$complement\n";
    }
    close(IN);
    &closeTag;
    system("cp $tmpdir/*.seq.$number $cwd/$target/sort_uniq/ && rm $cwd/$target/sort_uniq/$target.read.$number && rm -r $tmpdir");
}

sub splitRead{
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
    system("mkdir $cwd/$target/sort_uniq");
    $fcount = 1;
    $fname = "$target.read.$fcount";
    open(OUT, "> $cwd/$target/sort_uniq/$fname");
    open(IN, $cmd);
    while(<IN>){
	if ($count == 1 and !/N/){
	    print OUT $_;
	    $total ++;
	    if ($total % 10000000 == 0){
		close(OUT);
		$fcount ++;
		$fname = "$target.read.$fcount";
		open(OUT, "> $cwd/$target/sort_uniq/$fname");
	    }
	}elsif($count == 4){
	    $count = 0;
	}
	$count++;
    }
    close(IN);
    close(OUT);
    open(OUT, "> $cwd/$target/sort_uniq/splitRead.done");
    print OUT $fcount;
    close(OUT);
}
