#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

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
    $cwd = `pwd`;
    chomp($cwd);
}elsif($ENV{target} ne ""){
    $target    = $ENV{target};
    $control   = $ENV{control};
    $ref       = $ENV{ref};
    $number    = $ENV{number};
    $type      = $ENV{type};
    $tmpdir    = $ENV{tmpdir};
    $cwd       = $ENV{PBS_O_WORKDIR};
}else{
    print $usage;
    exit;
}

if($tmpdir eq ""){
    $tmpdir = ".";
    $ref_path = "$cwd/$ref";
    $workdir = "$cwd/$target";
}else{
    if (! -d $tmpdir){
	print "$tmpdir is not directory.

$usage";
	exit;
    }
    $ref_path = "$tmpdir/$ref";
    if (! -e $ref_path){
	system("mkdir $ref_path");
    }
    system("/usr/bin/rsync -a $cwd/$ref/$ref.sort_uniq $cwd/$ref/chr* $ref_path");
    $workdir = "$tmpdir/$target";
    if (-d $workdir){
	system("rm -r $workdir");
    }
    system("mkdir $workdir");
}

if ($control eq "default" or $control eq ""){
    $control = "$ref_path/$ref.sort_uniq";
}else{
    $control = "$cwd/$control/$control.sort_uniq";
}

if (-e "$cwd/$target/$target.sort_uniq.gz"){
    $file_type = "gz";
    $target_sort_uniq = "$cwd/$target/$target.sort_uniq.gz";
}

$number = "01" if $number eq "";

@nuc = ('A', 'C', 'G', 'T');

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

foreach $i (@chr){
    my $chr_file = "$ref_path/chr$i";
    open (IN, $chr_file);
    ($chr{$i} = <IN>) =~ y/a-z/A-Z/;
    close(IN);
}

if ($file_type eq "gz"){
    open(IN, "zcat $target_sort_uniq 2> /dev/null |");
}else{
    open(IN, "$cwd/$target/$target.sort_uniq");
}
while(<IN>){
    chomp;
    $length = length($_);
    last;
}
close(IN);

open(IN, $control);
while(<IN>){
    chomp;
    $clength = length($_);
    last;
}
close(IN);

&openTag;
if ($type eq "vcf"){
    open(IN,  "$cwd/$target/$target.vcf");
}elsif ($type eq "bi"){
    open(IN, "$cwd/$target/$target.snp.$number");
 }elsif ($type eq "kmer"){
    open(IN, "$cwd/$target/$target.map.$number");
}
while(<IN>){
    chomp;
    @row = split;
    if ($type eq "vcf"){
	next if  length($row[3]) != 1 or length($row[4]) != 1;
	$row[0] =~ s/chr//g;
	$row[0] += 0;
	$chr = $row[0];
	$pos = $row[1];
	$rf = $row[3];
	$alt = $row[4];
	$count = 20;
    }elsif ($type eq "bi"){
	$chr = $row[0];
	$pos = $row[1];
	$rf = $row[2];
	$alt = $row[3];
	$count = $row[4];
    }elsif ($type eq "kmer"){
	$chr = $row[0];
	$pos = $row[1];
	$rf = $row[4];
	$alt = $row[5];
	$alt =~ s/$rf//;
	$count = 5;
    }

    $ref_seq = substr($chr{$chr}, $pos - $length, $length * 2 -1);
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
if ($file_type eq "gz"){    
    system("zcat $target_sort_uniq 2> /dev/null | join - target.snp.st.$number | cut -d ' ' -f 2- > target.snp.$number && rm *.snp.sort.$number target.snp.st.$number");
}else{
    system("join $cwd/$target/$target.sort_uniq target.snp.st.$number | cut -d ' ' -f 2- > target.snp.$number && rm *.snp.sort.$number target.snp.st.$number");
}
system("sort -T $tmpdir target.snp.$number| uniq -c > target.snp.count.$number && rm target.snp.$number");
    
&openTag;
if ($type eq "vcf"){
    open(IN,  "$cwd/$target/$target.vcf");
}elsif ($type eq "bi"){
    open(IN, "$cwd/$target/$target.snp.$number");
 }elsif ($type eq "kmer"){
    open(IN, "$cwd/$target/$target.map.$number");
}
while(<IN>){
    chomp;
    @row = split;
    if ($type eq "vcf"){
	next if  length($row[3]) != 1 or length($row[4]) != 1;
	$row[0] =~ s/chr//g;
	$row[0] += 0;
	$chr = $row[0];
	$pos = $row[1];
	$rf = $row[3];
	$alt = $row[4];
	$count = 20;
    }elsif ($type eq "bi"){
	$chr = $row[0];
	$pos = $row[1];
	$rf = $row[2];
	$alt = $row[3];
	$count = $row[4];
    }elsif ($type eq "kmer"){
	$chr = $row[0];
	$pos = $row[1];
	$rf = $row[4];
	$alt = $row[5];
	$alt =~ s/$rf//;
	$count = 5;
    }
    $ref_seq = substr($chr{$chr}, $pos - $clength, $clength * 2 -1);
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

system("cat *.snp.sort.$number | join $control - | cut -d ' ' -f 2- > control.snp.$number");
system("sort -T $tmpdir control.snp.$number| uniq -c > control.snp.count.$number && rm *.snp.sort.$number control.snp.$number");

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
	print OUT "$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$cw\t$cm\t$tw\t$tm\t$genotype\n";
    }elsif ($type eq "kmer"){
	print OUT "$_\t$cw\t$cm\t$tw\t$tm\t$genotype\n";
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

sub closeTag{
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		close($tag);
	    }
	}
    }
}

sub sortTag{
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		system("sort -T $tmpdir $workdir/$tag.snp.tmp.$number > $workdir/$tag.snp.sort.$number && rm $workdir/$tag.snp.tmp.$number");
	    }
	}
    }
}

sub complement{
    my $seq = shift;
    my @seq = split('', $seq);
    my $i = length($seq);
    my $out = "";
    while($i > 0){
        $i--;
        if ($seq[$i] eq "A"){
            $out .= "T";
        }elsif($seq[$i] eq "C"){
            $out .= "G";
        }elsif($seq[$i] eq "G"){
            $out .= "C";
        }elsif($seq[$i] eq "T"){
            $out .= "A";
        }else{
            $out .= "N";
        }
    }
    return $out;
}

sub bynumber{
    $a <=> $b;
}
