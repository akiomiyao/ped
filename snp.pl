#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     snp.pl - detect snp candidate from lbc data. 

e.g. qsub -v target=ERR194147,control=hg38,tag=AAA,tmpdir=/mnt/ssd,cutoff=10 snp.pl

     target is name of target.
     control is reference or name of control data.
     tmpdir can be ommited.
     cutoff is cutoff value of noise.
      Default of cutoff is 10, can be ommited.

Author: Akio Miyao <miyao@affrc.go.jp>

';

if($ENV{target} ne ""){
    $target    = $ENV{target};
    $control   = $ENV{control};
    $tag       = $ENV{tag};
    $cutoff    = $ENV{cutoff};
    $cwd       = $ENV{PBS_O_WORKDIR};
    $tmpdir    = $ENV{tmpdir};
    $workdir = "$cwd/$target";
    $controldir = "$cwd/$control";
}else{
    print $usage;
    exit;
}

if ($tmpdir eq ""){
    $tmpdir = $workdir;
}else{
    $tmpdir = "$tmpdir/$target";
    if (-e $tmpdir){
	system("rm -r $tmpdir");
    }
    system("mkdir $tmpdir");
    system("cp $workdir/$target.lbc.$tag.gz $tmpdir && cp $controldir/$control.lbc.$tag.gz $tmpdir");
}

chdir $tmpdir;

$cutoff = 10 if $cutoff eq "";

$target_file = "$target.lbc.$tag";
$control_file = "$control.lbc.$tag";
$output_file = "$target.snp.$tag";

@nuc = ('A', 'C', 'G', 'T');

system("gzip -d *.lbc.$tag.gz");
open(OUT, "> $output_file");
open(IN, "join $control_file $target_file |");
while(<IN>){
    my $rep = 0;
    my $a = "";
    my $b = "";
    my $nohita = 0;
    my $nohitb = 0;
    my $pola = 0;
    my $polb = 0;
    @row = split;
    for ($i = 1; $i <= 4; $i++){
	if($row[$i] > 100){
	    $rep++;
	}elsif ($row[$i] >= $cutoff){
	    $a .= $nuc[$i -1];
	    if ($row[$i + 4] <= 1){
		$pola++;
	    }
	}elsif($row[$i] <= 1){
	    $nohita++;
	}	

	if($row[$i + 4] > 100){
	    $rep++;    
	}elsif ($row[$i + 4] >= $cutoff){
	    $b .= $nuc[$i -1];
	    if ($row[$i] <= 1){
		$polb++;
	    }
	}elsif($row[$i + 4] <= 1){
	    $nohitb++;
	}
    } 

    if ($a ne $b and $a ne "" and $b ne "" and $rep == 0 and $nohita >= 2 and $nohitb >= 2 and ($pola == 1 or  $polb == 1)){
	    print OUT "$row[0]\t$a\t$b\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\n";
    }
}
close(OUT);
close(IN);
if ($tmpdir ne $workdir){
    system("cp $output_file $workdir && rm -r $tmpdir");
}

