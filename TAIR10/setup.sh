#!/bin/sh
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

wget ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
perl mkChr.pl
perl mk20.pl
perl mkControlRead.pl
