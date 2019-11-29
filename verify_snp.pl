#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
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
    require './common.pl';
}

$usage = '
     verify_snp.pl - verification program for SNPs.


e.g. qsub -v target=ERR194147,control=default,ref=hg38,number=01,type=bi,tmpdir=/mnt/ssd verify_snp.pl

     For standalone machine
e.g. perl verify_snp.pl target control reference number type tmpdir
     perl verify_snp.pl ERR194147 default hg38 01 bi /mnt/ssd
     perl verify_snp.pl DRR054198 default IRGSP1.0 01 bi
     perl verify_snp.pl DRR054198 default IRGSP1.0 AAA kmer

target  : target name e.g. ERR194147
control : e.g. ERR194146, or \'default\' if you want to use reference data for control
ref     : reference genome. e.g. hg38
number  : specify number of splited subfile
type    : specify type of input data (bi, kmer or vcf)
tmpdir  : specify temporary directofy on local disk (can be ommited)

Author  : Akio Miyao <miyao@affrc.go.jp>
';

if ($ARGV[0] ne ""){
    $target  = $ARGV[0];
    $control = $ARGV[1];
    $ref     = $ARGV[2];
    $number  = $ARGV[3];
    $type    = $ARGV[4];
    $tmpdir  = $ARGV[5];
}elsif($ENV{target} ne ""){
    $target    = $ENV{target};
    $control   = $ENV{control};
    $ref       = $ENV{ref};
    $number    = $ENV{number};
    $type      = $ENV{type};
    $tmpdir    = $ENV{tmpdir};
}else{
    print $usage;
    exit;
}

if($tmpdir eq ""){
    $tmpdir = ".";
    $refdir = "$cwd/$ref";
    $workdir = "$cwd/$target";
}else{
    if (! -d $tmpdir){
	print "$tmpdir is not directory.

$usage";
	exit;
    }
    $refdir = $tmpdir . "/pedtmp." . int(rand 1000000);
    system("mkdir $refdir");
    system("$rsync -a $cwd/$ref/sort_uniq $refdir");
    $workdir = $tmpdir . "/pedtmp." . int(rand 1000000);
    system("mkdir $workdir");
    system("$rsync -a $cwd/$target/sort_uniq $workdir");
    if (($control ne $ref) and ($control ne "default")){
	$controldir = $tmpdir . "/pedtmp." . int(rand 1000000);
	if(! -e "$tmpdir/$control"){
	    system("mkdir $controldir");
	}
	system("$rsync -a $cwd/$control/sort_uniq $controldir");
    }
}

$number = "01" if $number eq "";

chdir $workdir;

open(IN, "$cwd/config");
while(<IN>){
    chomp;
    @row = split;
    if($row[0] eq $ref && $row[1] eq "chromosome"){
	if ($row[3] != 0){
	    for ($i = $row[2]; $i <= $row[3]; $i++){
		push(@chr, $i);
	    }
	}
	if ($row[4] ne ""){
	    foreach ($i = 4; $i <= $#row; $i++){
		push(@chr, $row[$i]);
	    }
	}
    }
}
close(IN);

open(IN, "zcat sort_uniq/*.gz 2> /dev/null |");
while(<IN>){
    chomp;
    $length = length($_);
    last;
}
close(IN);

&openTag;

if($type eq "vcf"){
    open(IN,  "$cwd/$target/$target.vcf");
}elsif($type eq "bi"){
    open(IN,  "$cwd/$target/$target.snp.$number");
}elsif($type eq "kmer"){
    open(IN,  "sort $cwd/$target/$target.map.$number |");
}else{
    print $usage;
    exit;
}
while(<IN>){
    chomp;
    @row = split;
    if ($type eq "vcf"){
	if(length($row[3]) == 1 and length($row[4]) == 1){
	    $chr   = $row[0];
	    $pos   = $row[1];
	    $rf    = $row[3];
	    $alt   = $row[4];
	    $count = 20;
	}
    }elsif ($type eq "bi"){
	    $chr   = $row[0];
	    $pos   = $row[1];
	    $rf    = $row[2];
	    $alt   = $row[3];
	    $count = $row[4];
    }elsif ($type eq "kmer"){
	    $chr   = $row[0];
	    $pos   = $row[1];
	    $rf    = $row[4];
	    $alt   = $row[5];
	    $alt =~ s/$rf//;
	    $count = 5;
    }
    $chr_prev = $chr_name;
    ($chr_name = $chr) =~ s/^0+//;
    if ($chr_name ne $chr_prev){
	open(CHR, "$cwd/$ref/chr$chr_name");
	$chr_seq = <CHR>;
	close(CHR);
    }
    $ref_seq = substr($chr_seq, $pos - $length, $length * 2 -1);
    next if length($ref_seq) != $length * 2 -1;
    $head = substr($ref_seq, 0, $length-1);
    $tail = substr($ref_seq, $length, $length);
    $mut_seq = $head . $alt . $tail;
    for($h = 0; $h < 10; $h++){
	$k = $j + $h -5;
	if ($k >= 0 and $k <= $#dat){
	    ($ichr, $ipos, $irf, $ialt, $icount) = split(' ', $dat[$k]);
	    if ($icount > 2 and abs($pos - $ipos) < $length and $pos != $ipos){
		$i = $length - ($pos - $ipos) -1;
		$head = substr($mut_seq, 0, $i);
		$tail = substr($mut_seq, $i +1);
		$mut_seq = $head . $ialt . $tail;
	    }
	}
    }
    for($i = 0; $i < $length; $i++){
	$tw = substr($ref_seq, $i, $length);
	$tm = substr($mut_seq, $i, $length);
	$tag = substr($tw, 0, 3);
	print $tag "$tw\t$chr $pos $rf $alt tw\n";
	$tag = substr($tm, 0, 3);
	print $tag "$tm\t$chr $pos $rf $alt tm\n";
    }		
}
close(IN);
&closeTag;
&sortTag;

system("cat *.snp.sort.$number > target.snp.st.$number");
&waitFile("target.snp.st.$number");
system("rm *.snp.sort.$number");
system("zcat sort_uniq/*.gz 2> /dev/null |join - target.snp.st.$number | cut -d ' ' -f 2- > target.snp.$number");
&waitFile("target.snp.$number");
system("rm target.snp.st.$number");
if($ARGV[0] ne ""){
    &report("Verifying SNPs $number: Sorting for target.snp.count.$number");
}
system("sort -T . $sort_opt target.snp.$number| uniq -c > target.snp.count.$number");
&waitFile("target.snp.count.$number");
system("rm target.snp.$number");

if($ARGV[0] ne ""){
    &report("Verifying SNPs $number: Selecting reads containing SNP");
}
if ($control eq "default" or $control eq "" or $control eq $ref){
    $control = "$refdir/sort_uniq/*.gz";
}else{
    if ($tmpdir eq "."){
	$control = "$cwd/$control/sort_uniq/*.gz";
    }else{
	$control = "$controldir/sort_uniq/*.gz";
    }
}

open(IN, "zcat $control 2> /dev/null |");
while(<IN>){
    chomp;
    $clength = length($_);
    last;
}
close(IN);

&openTag;

if($type eq "vcf"){
    open(IN,  "$cwd/$target/$target.vcf");
}elsif($type eq "bi"){
    open(IN,  "$cwd/$target/$target.snp.$number");
}elsif($type eq "kmer"){
    open(IN,  "sort $cwd/$target/$target.map.$number |");
}else{
    print $usage;
    exit;
}
while(<IN>){
    chomp;
    @row = split;
    if ($type eq "vcf"){
	if(length($row[3]) == 1 and length($row[4]) == 1){
	    $chr   = $row[0];
	    $pos   = $row[1];
	    $rf    = $row[3];
	    $alt   = $row[4];
	    $count = 20;
	}
    }elsif ($type eq "bi"){
	    $chr   = $row[0];
	    $pos   = $row[1];
	    $rf    = $row[2];
	    $alt   = $row[3];
	    $count = $row[4];
    }elsif ($type eq "kmer"){
	    $chr   = $row[0];
	    $pos   = $row[1];
	    $rf    = $row[4];
	    $alt   = $row[5];
	    $alt =~ s/$rf//;
	    $count = 5;
    }
    $chr_prev = $chr_name;
    ($chr_name = $chr) =~ s/^0+//;
    if ($chr_name ne $chr_prev){
	open(CHR, "$cwd/$ref/chr$chr_name");
	$chr_seq = <CHR>;
	close(CHR);
    }
    $ref_seq = substr($chr_seq, $pos - $clength, $clength * 2 -1);
    next if length($ref_seq) != $clength * 2 -1;
    $head = substr($ref_seq, 0, $clength-1);
    $tail = substr($ref_seq, $clength, $clength);
    $mut_seq = $head . $alt . $tail;
    for($h = 0; $h < 10; $h++){
	$k = $j + $h -5;
	if ($k >= 0 and $k <= $#dat){
	    ($ichr, $ipos, $irf, $ialt, $icount) = split(' ', $dat[$k]);
	    if ($icount > 2 and abs($pos - $ipos) < $clength and $pos != $ipos){
		$i = $clength - ($pos - $ipos) -1;
		$head = substr($mut_seq, 0, $i);
		$tail = substr($mut_seq, $i +1);
		$mut_seq = $head . $ialt . $tail;
	    }
	}
    }
    for($i = 0; $i < $clength; $i++){
	$cw = substr($ref_seq, $i, $clength);
	$cm = substr($mut_seq, $i, $clength);
	$tag = substr($cw, 0, 3);
	print $tag "$cw\t$chr $pos $rf $alt cw\n";
	$tag = substr($cm, 0, 3);
	print $tag "$cm\t$chr $pos $rf $alt cm\n";
    }		
}
close(IN);
&closeTag;
&sortTag;

if($ARGV[0] ne ""){
    &report("Verifying SNPs $number: Selecting reads with control allele");
}
system("cat *.snp.sort.$number > snp.sort.$number");
system("zcat $control 2> /dev/null | join - snp.sort.$number| cut -d ' ' -f 2- > control.snp.$number");
&waitFile("control.snp.$number");
system("rm *.snp.sort.$number");
if($ARGV[0] ne ""){
    &report("Verifying SNPs $number: Sorting for control.snp.count.$number");
}
system("sort -T . $sort_opt control.snp.$number| uniq -c > control.snp.count.$number");
&waitFile("control.snp.count.$number");
system("rm control.snp.$number");

open(IN, "target.snp.count.$number");
while(<IN>){
    chomp;
    @row = split;
    if($row[5] eq "tw"){
	$tw{"$row[1]\t$row[2]\t$row[3]\t$row[4]"} = $row[0];
    }elsif($row[5] eq "tm"){
	$tm{"$row[1]\t$row[2]\t$row[3]\t$row[4]"} = $row[0];
    }
}
close(IN);

system("rm target.snp.count.$number");

open(IN, "control.snp.count.$number");
while(<IN>){
    chomp;
    @row = split;
    if($row[5] eq "cw"){
	$cw{"$row[1]\t$row[2]\t$row[3]\t$row[4]"} = $row[0];
    }elsif($row[5] eq "cm"){
	$cm{"$row[1]\t$row[2]\t$row[3]\t$row[4]"} = $row[0];
    }
}

system("rm control.snp.count.$number");

if ($type eq "vcf"){
    open(OUT, "> $workdir/$target.vcf.verify");
    open(IN, "$cwd/$target/$target.vcf");
}elsif ($type eq "bi"){
    open(OUT, "> $workdir/$target.snp.verify.$number");
    open(IN, "$cwd/$target/$target.snp.$number");
}elsif ($type eq "kmer"){
    open(OUT, "> $workdir/$target.kmer.verify.$number");
    open(IN, "$cwd/$target/$target.map.$number");
}
while(<IN>){
    chomp;
    @row = split;
    if ($type eq "vcf"){
	if ($row[0] =~ /^chr/){
	    $row[0] =~ s/chr//;
	    $row[0] += 0;
	    $cw = $cw{"$row[0]\t$row[1]\t$row[3]\t$row[4]"};
	    $cm = $cm{"$row[0]\t$row[1]\t$row[3]\t$row[4]"};
	    $tw = $tw{"$row[0]\t$row[1]\t$row[3]\t$row[4]"};
	    $tm = $tm{"$row[0]\t$row[1]\t$row[3]\t$row[4]"};
	}
    }elsif ($type eq "bi"){
	$cw = $cw{"$row[0]\t$row[1]\t$row[2]\t$row[3]"};
	$cm = $cm{"$row[0]\t$row[1]\t$row[2]\t$row[3]"};
	$tw = $tw{"$row[0]\t$row[1]\t$row[2]\t$row[3]"};
	$tm = $tm{"$row[0]\t$row[1]\t$row[2]\t$row[3]"};
    }elsif ($type eq "kmer"){
	$rf = $row[4];
	$alt = $row[5];
	$alt =~ s/$rf//;
	$cw = $cw{"$row[0]\t$row[1]\t$rf\t$alt"};
	$cm = $cm{"$row[0]\t$row[1]\t$rf\t$alt"};
	$tw = $tw{"$row[0]\t$row[1]\t$rf\t$alt"};
	$tm = $tm{"$row[0]\t$row[1]\t$rf\t$alt"};
    }
    $cw += 0;
    $cm += 0;
    $tw += 0;
    $tm += 0;
    $genotype = "";
    if ($cw >= 5 and $cm <= 1){
	if ($tm >= 5 and $tw <= 1){
	    $genotype = 'M';
	}elsif(($tm + $tw) > 0 and $tm / ($tm + $tw) > 0.3 and $tm / ($tm + $tw) < 0.7 and $cm <= 1 and $tw >= 5){
	    $genotype = "H";
	}
    }
    if ($type eq "vcf"){
	if($cw > 0){
	    print OUT "$_\t$cw\t$cm\t$tw\t$tm\t$genotype\n";
	}else{
	    print OUT "$_\n";
	}
    }elsif ($type eq "bi"){
	$row[0] =~ s/^0+//;
	$row[1] =~ s/^0+//;
	print OUT "$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$cw\t$cm\t$tw\t$tm\t$genotype\n";
    }elsif ($type eq "kmer"){
	$row[0] =~ s/^0+//;
	$row[1] =~ s/^0+//;
	print OUT "$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\t$row[9]\t$row[10]\t$row[11]\t$row[12]\t$row[13]\t$row[14]\t$cw\t$cm\t$tw\t$tm\t$genotype\n";
    }
}
close(OUT);

if ($tmpdir ne "."){
    if ($type eq "vcf"){
	system("cp $target.vcf.verify $cwd/$target");
    }elsif ($type eq "bi"){
	system("cp $target.snp.verify.$number $cwd/$target");
    }elsif ($type eq "kmer"){
	system("cp $target.kmer.verify.$number $cwd/$target");
    }
    system("rm -r $workdir");
    system("rm -r $refdir");
    system("rm -r $controldir") if $controldir ne "";
}

if (-e "$cwd/$target/snp.sort.$number"){
    system("rm $cwd/$target/snp.sort.$number");
}
if (-e "$cwd/$target/$target.snp.$number"){
    system("rm $cwd/$target/$target.snp.$number");
}

if($ARGV[0] ne ""){
    &report("Verifying SNPs $number: Done");
}

sub openTag{
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "> $workdir/$tag.snp.tmp.$number");
	    }
	}
    }
}

sub sortTag{
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		if($ARGV[0] ne ""){
		    &report("Verifying SNPs $number: Sorting for $tag.snp.sort.$number");
		}
		system("sort -T . $sort_opt $tag.snp.tmp.$number > $tag.snp.sort.$number");
		&waitFile("$tag.snp.sort.$number");
		system("rm $tag.snp.tmp.$number");
	    }
	}
    }
}
