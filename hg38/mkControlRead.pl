#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

open(OUT, "|sort -T ./ |uniq > hg38.sort_uniq");
for(1 .. 22, 'X', 'Y'){
    open(IN, "chr$_");
    $data = <IN>;
    close(IN);
    $i = 0;
    while(1){
	$read = substr($data, $i, 100);
	last if length($read) != 100;
	if ($read !~ /[MRWSYKVHDBN]/){
	    print OUT "$read\n";
	    print OUT &complement($read) . "\n";
	}
	$i += 2;
    }
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
