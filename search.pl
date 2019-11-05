#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     search.pl - search alignment of structural variation. 

e.g. perl search.pl target chr pos
       perl search.pl ERR194147 1 100631715

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] eq ""){
    print $usage;
    exit;
}

$target = $ARGV[0];
$chr    = $ARGV[1];
$pos    = $ARGV[2];

open(IN, "$target/$target.aln.0");
seek(IN, 3062691,0);
while(<IN>){
    exit if $flag and /^#/;
    print;
    $flag = 1;
}

