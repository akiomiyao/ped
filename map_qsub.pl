#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     map_qsub.pl - qsub launcher of map.pl. 

e.g. perl map_qsub.pl target ref tmpdir
     perl map_qsub.pl ERR194147 hg38 /mnt/ssd

     target is name of target.
     tmpdir can be ommited.
     cutoff (defult is 10) can be ommited.

Author: Akio Miyao <miyao@affrc.go.jp>

';

$target = $ARGV[0];
$ref = $ARGV[1];
$tmpdir = $ARGV[2];


if ($target eq ""){
    print $usage;
    exit;
}

@nuc = ('A', 'C', 'G', 'T');
foreach (@nuc){
    $tag[0] = $_;
    foreach (@nuc){
	$tag[1] = $_;
	foreach (@nuc){
	    $tag[2] = $_;
	    $tag = join('', @tag);
	    $cmd = "qsub -v target=$target,ref=$ref,tag=$tag,tmpdir=$tmpdir map.pl";
	    print "$cmd\n";
	    system($cmd);
	    sleep(1);
	}
    }
}
