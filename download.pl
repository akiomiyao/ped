#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     download.pl - download fastq data from SRA of NCBI and setup directory. 

e.g. perl download.pl ACCESSION Working_directory
     perl download.pl ERR194147 /mnt/ssd/data

     or 

     perl download.pl accession=ACCESSION,wd=Working_directory
     perl download.pl accession=ERR194147,wd=/mnt/ssd/data

     Working_directory can be omitted.

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] =~ /accession/){
    my @arg = split(',', $ARGV[0]);
    foreach (sort @arg){
	next if $_ eq "";
	my ($name, $val) = split('=', $_);
	$$name = $val;
    }
}elsif ($ARGV[0] eq ""){
    print "$usage";
    exit;
}else{
    $accession = $ARGV[0];
    $wd = $ARGV[1];
}

if ($wd eq ""){
    $wd = `pwd`;
    chomp($wd);
}

system("mkdir $wd/$accession");
system("mkdir $wd/$accession/read");
chdir "$wd/$accession/read";
system("fastq-dump --split-files -A $accession");
