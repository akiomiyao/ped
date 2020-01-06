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

if ($chr =~ /^[0-9]*$/){
    $chr = "000$chr";
    $chr = substr($chr, length($chr) - 3, 3);
}

$pos = "00000000000" . $pos;
$pos = substr($pos, length($pos) - 11, 11);

$size = -s "$target/$target.index";
open(INDEX, "$target/$target.index");
binmode(INDEX);
$top = 0;
$bottom = $size;
$middle = int($size / 2);
while($bottom - $top > 1){
    seek(INDEX, $middle, 0);
    read(INDEX, $data, 1000);
    foreach (split('\n', $data)){
	@row = split;
	if ($row[1] =~ /^[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]$/){
	    $ichr = $row[0];
	    $ipos = $row[1];
	    if ($chr eq $ichr){
		if ($pos eq $ipos){
		    &printData;
		    $hit = 1;
		}
	    }
	}
    }
    exit if $hit;
    if ($ichr gt $chr){
	$bottom = $middle;
    }elsif($ichr eq $chr){
	if ($ipos gt $pos){
	    $bottom = $middle;
	}else{
	    $top = $middle;
	}
    }else{
	$top = $middle;
    }
    $middle = int(($top + $bottom) / 2);
 }

sub printData{
    if ($row[2] eq "S"){
	open(IN, "$target/$target.aln");
    }else{
	open(IN, "$target/$target.aln.$row[2]");
    }
    seek(IN, $row[3], 0);
    while(<IN>){
	if ($flag and /^$/){
	    $flag = 0;
	    print "\n";
	    return;
	}
	print;
	$flag = 1;
    }
}
