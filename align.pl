#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     align.pl - bi-directional alignment program. 

e.g. perl align.pl target reference margin tmpdir

     qsub -v target=ERR194147,ref=hg38,margin=5,tmpdir=/mnt/ssd align.pl

     margin and tmpdir can be ommitted.

If each computer node have a local disk, specify the tmpdir to the local
disk is recommended. High speed disk like as SSD for local disk is prefered.
The size of local disk should be 1TB or more for analysis of human genome.
Use of local disk for temporary working directory in each node is recommended. 

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[0] ne ""){
    $target = $ARGV[0];
    $ref = $ARGV[1];
    $margin = $ARGV[2];
    $tmpdir = $ARGV[3];
    $cwd = `pwd`;
    chomp($cwd);
    $workdir = $target;
}elsif($ENV{target} ne ""){
    $target    = $ENV{target};
    $ref       = $ENV{ref};
    $margin    = $ENV{margin};
    $tmpdir    = $ENV{tmpdir};
    $cwd       = $ENV{PBS_O_WORKDIR};
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
    system("/usr/bin/rsync -a $cwd/$ref $tmpdir");
    $ref_path = "$tmpdir/$ref";
    $tmpdir = "$tmpdir/$target";
    if (-d $tmpdir){
	system("rm -r $tmpdir");
    }
    system("mkdir $tmpdir");
}else{
    $tmpdir = ".";
    $ref_path = "$cwd/$ref";
}

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

open(IN, "$target.sort_uniq");
while(<IN>){
    chomp;
    $length = length($_);
    last;
}
close(IN);

&openTag;

open(IN, "$target.sort_uniq");
while(<IN>){
    chomp;
    $head_pos = $margin + 1;
    $head = substr($_, $head_pos -1, 20);
    $tag = substr($head, 0, 3);
    print $tag "$head $_ $head_pos\n";
}    

&map;

&openTag;

open(IN, "cat $tmpdir/*.map.$margin|");
while(<IN>){
    chomp;
    $tail_pos = $length - 20 - $margin + 1;
    $tail = substr($_, $tail_pos - 1, 20);
    $tag = substr($tail, 0, 3);
    print $tag "$tail $_ $tail_pos\n";
}    

&map;

&analysis;

if ($tmpdir ne "."){
    if (-d $tmpdir){
	system ("rm -r $tmpdir");
    }
}else{
    system ("rm *.map.$margin");
}

sub openTag{
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "> $tmpdir/$tag.tmp.$margin")
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
		system("sort -T $tmpdir $tmpdir/$tag.tmp.$margin > $tmpdir/$tag.$margin && rm $tmpdir/$tag.tmp.$margin");
		system("zcat $ref_path/ref20_uniq.$tag.gz | join $tmpdir/$tag.$margin - |cut -c 22- > $tmpdir/$tag.map.$margin && rm $tmpdir/$tag.$margin");
	    }
	}
    }
}

sub analysis{
    open(OUT, "> $tmpdir/$target.aln.$margin");
    open(IN, "cat $tmpdir/*.map.$margin|");
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
	    }
	}
    }
    close(IN);
    close(OUT);
    if ($tmpdir ne "."){
	system("cp $tmpdir/$target.aln.$margin $workdir");
    }
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
