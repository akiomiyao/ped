#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     count_qsub.pl - qsub launcher of count.pl. 

e.g. perl count_qsub.pl target tmpdir

     target is name of target.
     tmpdir can be ommited.

Author: Akio Miyao <miyao@affrc.go.jp>

';

$target = $ARGV[0];
$tmpdir = $ARGV[1];
$cwd = `pwd`;
chomp($cwd);

if ($target eq ""){
    print $usage;
    exit;
}

opendir(DIR, "$target/tmp");
@file = grep(/sort_uniq/, readdir(DIR));
closedir(DIR);

foreach (sort @file){
    $number = (split('\.', $_))[2];
    if ($tmpdir eq ""){
	$cmd = "qsub -v target=$target,number=$number count.pl";
    }else{
	$cmd = "qsub -v target=$target,number=$number,tmpdir=$tmpdir count.pl";
    }
    print "$cmd\n";
    system($cmd);
    sleep 1;
}
