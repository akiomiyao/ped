#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#
open(IN, "zcat Gmax_275_v2.0.fa.gz|");
while(<IN>){
    chomp;
    if (/^>/){
	close(OUT);
	close(FA);
	($file = $_) =~ s/>//;
	$file =~ s/Chr/chr/;
	$file =~ s/chr0/chr/;
	$file =~ s/scaffold_/chr/;
	open(OUT, "> $file");
	open(FA, "> $file.fa");
	print FA ">$file\n";
    }else{
	y/a-z/A-Z/;
	print OUT;
	print FA "$_\n";
    }
}
close(IN);
close(OUT);
close(FA);
