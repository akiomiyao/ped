#!/bin/sh
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

if [ ! -f Gmax_275_v2.0.fa.gz ]; then
    echo
    echo Gmax_275_v2.0.fa.gz not found.
    echo Please download Gmax_275_v2.0.fa.gz in this directory.
    echo URL: https://phytozome.jgi.doe.gov/pz/portal.html
    echo
fi
perl mkChr.pl
perl mk20.pl
perl mkControlRead.pl
