#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

open(IN, "zcat dmel-all-chromosome-r6.26.fasta.gz|");
while(<IN>){
    chomp;
    if (/^>/){
	last if />2Cen/;
	close(OUT);
	close(FA);
	s/>//;
	@row = split;
	$file = "chr" . $row[0];
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
