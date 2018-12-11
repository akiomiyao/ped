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

e.g. perl download.pl ACCESSION
     perl download.pl ERR194147

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] eq ""){
    print "$usage";
    exit;
}

system("mkdir $ARGV[0]");
system("mkdir $ARGV[0]/read");
chdir "$ARGV[0]/read";
system("fastq-dump --split-files -A $ARGV[0]");
