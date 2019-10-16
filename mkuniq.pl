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
    require './common.pl';
}

if ($ENV{target} ne ""){
    $target    = $ENV{target};
    $tag      = $ENV{tag};
}else{
    exit;
}

system("sort -T $cwd/$target $sort_opt $cwd/$target/ref20.$tag > $cwd/$target/ref20_sort.$tag");
&waitFile("$cwd/$target/ref20_sort.$tag");
open(OUT, "|gzip -f > $cwd/$target/ref20_uniq.$tag.gz");
open(IN, "$cwd/$target/ref20_sort.$tag");
while(<IN>){
    chomp;
    @row = split;
    if ($prev ne "" and $prev ne $row[0]){
	print OUT "$pline\n" if $count == 1;
	$count =0;
    }
    $prev = $row[0];
    $pline = $_;
    $count++;
}
close(IN);
close(OUT);
&waitFile("$cwd/$target/ref20_uniq.$tag.gz");
system("rm $cwd/$target/ref20.$tag $cwd/$target/ref20_sort.$tag");

