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

$target = "$ARGV[0]/$ARGV[0].*.sort";

open(IN, "cat $target | grep $ARGV[2] |");
while(<IN>){
    chomp;
    @row = split;
    if ($row[1] eq $ARGV[1] and $row[2] == $ARGV[2]){
	s/\t/\n/g;
	print;
    }
}
