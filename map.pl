#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     map.pl - map SNPs from snp file. 

e.g. qsub -v target=ERR194147,ref=hg38,tag=AAA,tmpdir=/mnt/ssd map.pl

     For standalone machine
e.g. perl map.pl DRR054198 IRGSP1.0 
     perl map.pl SRR8181712 TAIR10

     target is name of target.
     ref is name of reference.
     tmpdir can be ommited.

Author: Akio Miyao <miyao@affrc.go.jp>

';

if($ENV{target} ne ""){
    $target      = $ENV{target};
    $ref         = $ENV{ref};
    $tag         = $ENV{tag};
    $cwd         = $ENV{PBS_O_WORKDIR};
    $cwd         = $ENV{SGE_O_WORKDIR} if $ENV{SGE_O_WORKDIR} ne "";
    $tmpdir      = $ENV{tmpdir};
    $workdir     = "$cwd/$target";
    $refdir      = "$cwd/$ref";
    $target_file = "$target.snp.$tag";
    $ref_file    = "zcat ref20_uniq.$tag.gz|";
    $output_file = "$target.map.$tag";
}elsif($ARGV[0] ne ""){
    $cwd       = `pwd`;
    chomp($cwd);
    $target      = $ARGV[0];
    $ref         = $ARGV[1];
    $workdir     = "$cwd/$target";
    $refdir      = "$cwd/$ref";
    $target_file = "cat $workdir/$target.snp.[ACGT][ACGT][ACGT]|";
    $output_file = "$target.map.AAA";
    $ref_file    = "zcat $refdir/ref20_uniq.*.gz|";
}else{
    print $usage;
    exit;
}

if ($tmpdir eq ""){
    $tmpdir = $workdir;
    $ref_file    = "zcat $refdir/ref20_uniq.*.gz|";
}else{
    $tmpdir = "$tmpdir/$target";
    if (-e $tmpdir){
	system("rm -r $tmpdir");
    }
    system("mkdir $tmpdir");
    system("cp $workdir/$target.snp.$tag $tmpdir && cp $refdir/ref20_uniq.$tag.gz $tmpdir");
}

chdir $tmpdir;

$s = {};
open(OUT, "> $output_file");
open(IN, $target_file);
while(<IN>){
    chomp;
    @row = split;
    $seq = shift(@row);
    if (length($row[0]) == 1){
	$dat{$seq} = "$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\t$row[9]";
    }
}

open(IN, $ref_file);
while(<IN>){
    chomp;
    @row = split;
    $seq = substr($row[0], 0, 19);
    $nuc = substr($row[0], 19, 1);
    if ($dat{$seq}){
	@dat = split('\t', $dat{$seq});
	if($row[3] eq "f"){
	    $pos = $row[2] + 19;
	    if ($nuc eq $dat[0]){
		print OUT "$row[1]\t$pos\t$seq\t$nuc\t$dat[0]\t$dat[1]\t$row[3]\t$dat[2]\t$dat[3]\t$dat[4]\t$dat[5]\t$dat[6]\t$dat[7]\t$dat[8]\t$dat[9]\n";
	    }
	}else{
	    $pos = $row[2] - 19;
	    $dat[0] = complement($dat[0]);
	    $dat[1] = complement($dat[1]);
	    $seq = complement($seq);
	    $nuc = complement($nuc);
	    if ($nuc eq $dat[0]){
		print OUT "$row[1]\t$pos\t$seq\t$nuc\t$dat[0]\t$dat[1]\t$row[3]\t$dat[2]\t$dat[3]\t$dat[4]\t$dat[5]\t$dat[6]\t$dat[7]\t$dat[8]\t$dat[9]\n";
	    }
	}
    }
}
close(OUT);

if ($tmpdir ne $workdir){
    system("cp $target.map.$tag $workdir && rm -r $tmpdir");
}

sub complement{
    my $seq = shift;
    my @seq = split('', $seq);
    my $i = 0;
    my $out = "";
    for ($i = length($seq) - 1 ; $i >= 0; $i--){
        if ($seq[$i] eq "A"){
            $out .= "T";
        }elsif($seq[$i] eq "C"){
            $out .= "G";
        }elsif($seq[$i] eq "G"){
            $out .= "C";
        }elsif($seq[$i] eq "T"){
            $out .= "A";
        }
    }
    return $out;
}

sub bynumber{
    $a <=> $b;
}
