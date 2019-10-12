#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2019 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$target    = $ENV{target};
$tag       = $ENV{tag};
$cwd       = $ENV{PBS_O_WORKDIR};

system("sort -S 100M -T . $cwd/$target/sort_uniq/$tag.seq | uniq | gzip > $cwd/$target/sort_uniq/$target.sort_uniq.$tag.gz && rm $cwd/$target/sort_uniq/$tag.seq");

