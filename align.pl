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
     align.pl - bidirectional alignment program. 

     For example,
     perl align.pl target reference margin tmpdir

     qsub -v target=ERR194147,ref=hg38,margin=5,tmpdir=/mnt/ssd align.pl

     margin and tmpdir is optional, can be ommitted.

     If each computer node have a local disk, specify the tmpdir to the local
     disk is recommended.
     High speed disk like as SSD for local disk is prefered.
     The size of local disk should be 1TB or more for analysis of human genome.

     Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] ne ""){
    $target = $ARGV[0];
    $ref = $ARGV[1];
    $margin = $ARGV[2];
    $tmpdir = $ARGV[3];
    $workdir = $target;
}elsif($ENV{target} ne ""){
    $target    = $ENV{target};
    $ref       = $ENV{ref};
    $margin    = $ENV{margin};
    $tmpdir    = $ENV{tmpdir};
    $workdir = "$cwd/$target";
}else{
    print $usage;
    exit;
}

if ($ref eq ""){
    print "reference is empty.
$usage";
    exit;
}

$margin = 0 if $margin eq "";

if ($tmpdir ne ""){
    $refdir = $tmpdir . "/pedtmp." . int(rand 1000000);
    system("$rsync -a $cwd/$ref/ $refdir");
    $targetdir = $tmpdir . "/pedtmp." . int(rand 1000000);
    system("mkdir $targetdir");
    system("$rsync -a $cwd/$target/sort_uniq $targetdir");
}else{
    $targetdir = "$cwd/$target";
    $refdir = "$cwd/$ref";
}

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

if ($chr[0] eq ""){
    opendir(REF, $refdir);
    foreach(readdir(REF)){
	chomp;
	if (/^chr/){
	    ($chr = $_) =~ s/^chr//; 
	    push(@chr, $chr);
	}
    }
    closedir(REF);
}

if($ARGV[0] ne ""){
    &report("Aligning of $target sequence to $ref genome. margin = $margin: Reading genome sequeces");
}
foreach $i (@chr){
    my $chr_file = "$refdir/chr$i";
    open (IN, $chr_file);
    ($chr{$i} = <IN>) =~ y/a-z/A-Z/;
    close(IN);
}

open(IN, "zcat $targetdir/sort_uniq/*.gz 2> /dev/null |");
while(<IN>){
    chomp;
    $length = length($_);
    last;
}
close(IN);

&openTag;

if($ARGV[0] ne ""){
    &report("Aligning of $target sequence to $ref genome. margin = $margin: Dividing the sort_uniq sequence");
}

open(IN, "zcat $targetdir/sort_uniq/*.gz 2> /dev/null |");
while(<IN>){
    chomp;
    $head_pos = $margin + 1;
    $head = substr($_, $head_pos -1, 20);
    $tag = substr($head, 0, 3);
    print $tag "$head $_ $head_pos\n";
}    

&map;

&openTag;

if($ARGV[0] ne ""){
    &report("Aligning of $target sequence to $ref genome. margin = $margin: Dividing the 5'-mapped sequences");
}

open(IN, "cat $targetdir/*.map.$margin|");
while(<IN>){
    chomp;
    $tail_pos = $length - 20 - $margin + 1;
    $tail = substr($_, $tail_pos - 1, 20);
    $tag = substr($tail, 0, 3);
    print $tag "$tail $_ $tail_pos\n";
}    

&map;

if($ARGV[0] ne ""){
    &report("Aligning of $target sequence to $ref genome. margin = $margin: Output alignments");
}

&analysis;

if ($tmpdir ne ""){ 
    system("cp $targetdir/$target.aln.$margin $cwd/$target");
    system ("rm -r $targetdir $refdir");
}

if($ARGV[0] ne ""){
    &report("Aligning of $target sequence to $ref genome. margin = $margin: Done");
}

sub openTag{
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "> $targetdir/$tag.tmp.$margin")
	    }
	}
    }
}

sub map{
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		close($tag);
		if($ARGV[0] ne ""){
		    &report("Aligning of $target sequence to $ref genome. margin = $margin: Mapping for the $tag.$margin subfile");
		}
		system("sort $sort_opt -T $targetdir $targetdir/$tag.tmp.$margin > $targetdir/$tag.$margin");
		&waitFile("$targetdir/$tag.$margin");
		system("zcat $refdir/ref20_uniq.$tag.gz | join $targetdir/$tag.$margin - |cut -c 22- > $targetdir/$tag.map.$margin");
		&waitFile("$targetdir/$tag.map.$margin");
		system("rm $targetdir/$tag.tmp.$margin $targetdir/$tag.$margin");
	    }
	}
    }
}

sub analysis{
    $total = 0;
    open(OUT, "> $targetdir/$target.aln.$margin");
    open(IN, "cat $targetdir/*.map.$margin|");
    while(<IN>){
	chomp;
	($seq, $hpos, $hchr, $head_pos, $head_direction, $tpos, $tchr, $tail_pos, $tail_direction) = split;
	$length = length($seq);
	$head = substr($_, $hpos -1, 20);
	$tail = substr($_, $tpos -1, 20);
	$hhit = 0;
	$thit = 0;
	
	if ($head_direction eq "f"){
	    $head_seq = substr($chr{$hchr}, $head_pos - $margin -1, $length);
	    next if $head_seq eq $seq;
	    if ($hchr eq $tchr and $tail_direction eq "f"){
		$distance = $tail_pos - $head_pos;
		if ($distance == $tpos - $hpos){
		    &snp;
		    next;
		}
	    }else{
		$distance = "";
	    }
	    if ($tail_direction eq "f"){
		$tail_seq = substr($chr{$tchr}, $tail_pos - $length + 20 + $margin -1, $length);
	    }else{
		$tail_seq = &complement(substr($chr{$tchr}, $tail_pos - 20 -$margin, $length));
	    }
	    
	    $head_bar = "";
	    $tail_bar = "";
	    @head = split('', $head_seq);
	    @seq = split('', $seq);
	    @tail = split('', $tail_seq);
	    @head_bar = ();
	    @tail_bar = ();
	    $head_space = "";
	    $tail_space = "";
	    for($i = 0; $i < $length; $i++){
		if ($seq[$i] eq $head[$i]){
		    $head_bar .= "|";
		    $head_bar[$i] = "|";
		}else{
		    $head_bar .= " ";
		}
		if ($seq[$i] eq $tail[$i]){
		    $tail_bar .= "|";
		    $tail_bar[$i] = "|";
		}else{
		    $tail_bar .= " ";
		}
	    }
	    
	    $head_junction = "";
	    $head_fail = 0;
	    $head_space = "        ";
	    for($i = 10; $i < $length; $i++){
		$unmatch = 0;
		$head_space .= " ";
		if ($head_bar[$i] eq ""){
		    for ($j = 1; $j < 6; $j++){
			if ($head_bar[$i + $j] eq ""){
			    $unmatch ++
			}
		    }
		    if ($unmatch >= 2){
			$hcount = 0;
			for ($k = 0; $k < $i; $k++){
			    if($head_bar[$k] eq "|"){
				$hcount++;
			    }
			}
			if ($hcount > $i -1 and $i > $hpos + 20 + 5){
			    $head_junction = $head_pos + $i - $margin;
			}
			last;
		    }
		}
	    }
	    $tail_junction = "";
	    $tail_fail = 0;
	    for($i = 0; $i < $length - 10; $i++){
		$tail_space .= " ";
	    }
	    for($i = $length - 10; $i > 0; $i--){
		$unmatch = 0;
		chop($tail_space);
		if ($tail_bar[$i] eq ""){
		    for ($j = 1; $j < 6; $j++){
			if ($tail_bar[$i - $j] eq ""){
			    $unmatch ++
			}
		    }
		    if ($unmatch >= 2){
			$tcount = 0;
			for ($k = $length; $k > $i; $k--){
			    if($tail_bar[$k] eq "|"){
				$tcount++
			    }
			}
			if ($i < $length -20 - $margin - 5){
			    if ($tail_direction eq "f"){
				$tail_junction = $tail_pos - $length + $i + 20 + $margin;
			    }else{
				$tail_junction = $tail_pos + $length - $i - 20 - $margin;
			    }
			}
			last;
		    }
		}
	    }
	    
	    $type = "";
	    if ($head_junction ne "" and $tail_junction ne ""){
		if ($distance ne ""){
		    $distance = $length - 20 - $distance - $margin * 2;
		    if ($distance > 0 ){
			$type = "insertion";
		    }elsif ($distance < 0){
			$type = "deletion";
		    }
		}else{
		    if ($hchr ne $tchr){
			$type = "translocation";
		    }elsif ($tail_direction eq "r"){
			$type = "inversion";
		    }
		}
		print  OUT "# $hchr $head_junction $tchr $tail_junction $tail_direction $type $distance
$head $tail $margin
$head_space Chr$hchr $head_junction
$head_space |
$head_seq
$head_bar
$seq
$tail_bar
$tail_seq
$tail_space |
$tail_space Chr$tchr $tail_junction

";
		$total ++;
		if ($total % 1000000 == 0){
		    report("Aligning of $target sequence to $ref genome. margin = $margin: Output alignments. $total alignments are processed.");
		}
	    }
	}
    }
    close(IN);
    close(OUT);
    report("Aligning of $target sequence to $ref genome. margin = $margin: Output alignments. $total alignments are processed.");
    system("rm $targetdir/*.map.$margin");
}

sub snp{
    $head_bar = "";
    @head = split('', $head_seq);
    @seq = split('', $seq);
    $mcount = 0;
    $out = "";
    for($i = 0; $i < $length; $i++){
	if ($seq[$i] eq $head[$i]){
	    $head_bar .= "|";
	}else{
	    $head_bar .= " ";
	    $pos = $head_pos + $i - $margin;
	    if ($i >= $margin and $i < $length - $margin){
		$out .= "# $hchr $pos snp $head[$i] $seq[$i]\n";
		$mcount++;
	    }
	}
    }
    if ($mcount > 0 and $mcount <= 5){
	print OUT $out . "$head_seq
$head_bar
$seq

";
	$total ++;
	if ($total % 1000000 == 0){
	    report("Aligning of $target sequence to $ref genome. margin = $margin: Output alignments. $total alignments are processed.");
	}
    }
}
