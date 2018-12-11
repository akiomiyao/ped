#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     split_count.pl - split and compress data from count.pl 

e.g. perl split_count.pl input_file workdir
     This script is called by count.pl.

Author: Akio Miyao <miyao@affrc.go.jp>
';

$target   = $ENV{target};
$number   = $ENV{number};
$tmpdir   = $ENV{tmpdir};
$cwd      = $ENV{PBS_O_WORKDIR};
$workdir = "$cwd/$target";

if ($target eq ""){
    print $usage;
    exit;
}

while(<>){
    $head = substr($_, 0, 3);
    if ($prev ne $head){
	$output_file = "$target.count.$head.$number.gz";
	open(OUT, "|gzip > $output_file");
    }
    print OUT "$_";
    $prev = $head;
}
close(OUT);

if ($tmpdir ne ""){
    system("cp $target.count* $workdir && rm -r $tmpdir/$target");
}
system("rm $workdir/tmp/$target.sort_uniq.$number");
