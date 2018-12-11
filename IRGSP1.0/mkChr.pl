#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#
open(IN, "zcat IRGSP-1.0_genome.fasta.gz|");
while(<IN>){
    chomp;
    if (/^>/){
	close(OUT);
	close(FA);
	if (!/_|M/){
	    ($file = $_) =~ s/>//;
	    $file =~ s/chr0/chr/;
	    open(OUT, "> $file");
	    open(FA, "> $file.fa");
	    print FA ">$file\n";
	}
    }else{
	y/a-z/A-Z/;
	print OUT;
	print FA "$_\n";
    }
}
close(IN);
close(OUT);
close(FA);
