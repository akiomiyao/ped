#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     merge_qsub.pl - launcher of merge.pl. 

e.g. perl merge_qsub.pl target tmpdir
     perl merge_qsub.pl ERR194147 /mnt/ssd

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
foreach (grep(/count/, readdir(DIR))){
    @row = split('\.', $_);
    $tag{$row[$#row - 2]} = 1;
}
closedir(DIR);

foreach $tag (sort keys %tag){
    $cmd = "qsub -v target=$target,tag=$tag,tmpdir=$tmpdir merge.pl";
    print "$cmd\n";
    system($cmd);
    sleep(1);
}
