#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2019 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

if ($ENV{PBS_O_WORKDIR} ne ""){
    $cwd = $ENV{PBS_O_WORKDIR};
    chdir $cwd;
    require "$cwd/common.pl";
}elsif($ENV{SGE_O_WORKDIR} ne ""){
    $cwd = $ENV{SGE_O_WORKDIR};
    chdir $cwd;
    require "$cwd/common.pl";
}else{
    $cwd = `pwd`;
    chomp($cwd);
    require './common.pl';
}
$target    = $ENV{target};
$tag       = $ENV{tag};
$tmpdir    = $ENV{tmpdir};

if ($tmpdir eq ""){
    system("sort $sort_opt -T $cwd/$target $cwd/$target/sort_uniq/$tag.seq | uniq | gzip > $cwd/$target/sort_uniq/$target.sort_uniq.$tag.gz && rm $cwd/$target/sort_uniq/$tag.seq && touch $cwd/$target/done.$tag");
}else{
    $tmpdir = "$tmpdir/$target";
    system("mkdir $tmpdir");
    system("cp $cwd/$target/sort_uniq/$tag.seq $tmpdir && sort $sort_opt -T $tmpdir $tmpdir/$tag.seq | uniq | gzip > $tmpdir/$target.sort_uniq.$tag.gz && sleep 10 && cp $tmpdir/$target.sort_uniq.$tag.gz $cwd/$target/sort_uniq/ && rm $cwd/$target/sort_uniq/$tag.seq  && rm -r $tmpdir && touch $cwd/$target/done.$tag");
}


