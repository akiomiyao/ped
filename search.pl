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
or   perl search.pl target=ERR194147,chr=1,pos=100631715
     perl search.pl target=ERR194147,chr=1,pos=100631715,wd=/work
     wd (working directory) is optional.

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] =~ /target/){
    my @arg = split(',', $ARGV[0]);
    foreach (sort @arg){
	next if $_ eq "";
	my ($name, $val) = split('=', $_);
	$$name = $val;
    }
}elsif ($ARGV[0] eq ""){
    print $usage;
    exit;
}else{
    $target = $ARGV[0];
    $chr    = $ARGV[1];
    $pos    = $ARGV[2];
}

if ($chr =~ /^[0-9]*$/){
    $chr = "000$chr";
    $chr = substr($chr, length($chr) - 3, 3);
}

$pos = "00000000000" . $pos;
$pos = substr($pos, length($pos) - 11, 11);

if ($wd eq ""){
    $wd = ".";
}

$size = -s "$wd/$target/$target.index";
open(INDEX, "$wd/$target/$target.index");
binmode(INDEX);
$top = 0;
$bottom = $size;
$middle = int($size / 2);
while($bottom - $top > 1){
    seek(INDEX, $middle - 100, 0);
    read(INDEX, $data, 1000);
    foreach (split('\n', $data)){
	@row = split;
	if ($row[1] =~ /^[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]$/){
	    $ichr = $row[0];
	    $ipos = $row[1];
	    if ($chr eq $ichr){
		if ($pos eq $ipos){
		    $dat = &printData;
		    $output{$dat} = 1;
		    $hit = 1;
		}
	    }
	}
    }
    foreach (sort keys %output){
	print $_;
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
    my (@tmp, $output);
    if ($row[2] eq "S"){
	open(IN, "$wd/$target/$target.aln");
    }else{
	open(IN, "$wd/$target/$target.aln.$row[2]");
    }
    $row[3] -= 1000;
    $row[3] = 0 if $row[3] < 0;
    seek(IN, $row[3], 0);
    $flag = 0;
    while(<IN>){
	if (/^$/){
	    if ($flag){
		return $output;
	    }
	    $output = "";
	}
	@tmp = split;
	$flag = 1 if $tmp[2] == $pos;
	$output .= $_;
    }
}
