#!/bin/sh
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

echo Making reference data set of Drosophila_melanogaster r6.26.
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.26_FB2019_01/fasta/dmel-all-chromosome-r6.26.fasta.gz
echo Making chromosome data.
perl mkChr.pl
echo Making 20mer reference.
perl mk20.pl
echo Making control sort_uniq file.
perl mkControlRead.pl
echo Done.
