#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     lbc_qsub.pl - launcher of last_base_count.pl. 

e.g. perl lbc_qsub.pl target tmpdir
     perl lbc_qsub.pl ERR194147 /mnt/ssd

     target is name of target.
     tmpdir can be ommited.

Author: Akio Miyao <miyao@affrc.go.jp>

';

$target = $ARGV[0];
$tmpdir = $ARGV[1];

if($target eq ""){
    print $usage;
    exit;
}

opendir(DIR, $target);
@file = grep(/count/, readdir(DIR));
closedir(DIR);

foreach $file (sort @file){
    ($target, $tag) = (split('\.', $file))[0, 2];
    $cmd = "qsub -v target=$target,tag=$tag,tmpdir=$tmpdir last_base_count.pl";
    print "$cmd\n";
    system($cmd);
    sleep(1);
}
