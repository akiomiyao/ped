#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     verify_snp_qsub.pl - qsub launcher for verify_snp.pl.

e.g.
     perl verify_snp_qsub.pl target control reference type tmpdir

     perl verify_snp_qsub.pl ERR194147 default hg38 bi /mnt/ssd

target  : target name e.g. ERR194147
control : e.g. ERR194146, or \'default\' if you want to use reference data for control
ref     : reference genome. e.g. hg38
number  : specify number of splited subfile
type    : specify type of input data (bi, kmer or vcf)
tmpdir  : specify temporary directofy on local disk (can be ommited)

Author  : Akio Miyao <miyao@affrc.go.jp>
';


$target  = $ARGV[0];
$control = $ARGV[1];
$ref     = $ARGV[2];
$type    = $ARGV[3];
$tmpdir  = $ARGV[4];

if ($target eq ""){
    print $usage;
    exit;
}

opendir(DIR, $target);
while(readdir(DIR)){
    if ($type eq "kmer"){
	if (/$target.map.[ACGT][ACGT][ACGT]$/){
	    $number = (split('\.', $_))[2];
	    push(@qsub, "qsub -v target=$target,control=$control,ref=$ref,number=$number,type=$type,tmpdir=$tmpdir verify_snp.pl");
	}
    }else{
	if (/$target.snp.[0-9][0-9]$/){
	    $number = (split('\.', $_))[2];
	    push(@qsub, "qsub -v target=$target,control=$control,ref=$ref,type=$type,number=$number,tmpdir=$tmpdir verify_snp.pl");
	}
    }
}

foreach (sort @qsub){
    print "$_\n";
    system($_);
    sleep 1;
}
