#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

@nuc = ('A', 'C', 'G', 'T');

if ($ARGV[0] ne "retry"){
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		open($tag, "> ref20.$tag");
	    }
	}
    }
    
    foreach $i (1 .. 19, 'X', 'Y'){
	&mk20mer($i);
    }
    
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		close($tag);
	    }
	}
    }
}

foreach $nuc (@nuc){
    $tag[0] = $nuc;
    foreach $nuc (@nuc){
	$tag[1] = $nuc;
	foreach $nuc (@nuc){
	    $tag[2] = $nuc;
	    $tag = join('', @tag);
	    &mkUniq($tag);
	}
    }
}

sub mkUniq{
    my $tag = shift;
    system("sort -T ./  ref20.$tag > ref20_sort.$tag");
    open(OUT, "|gzip -f > ref20_uniq.$tag.gz");
    open(IN, "ref20_sort.$tag");
    while(<IN>){
	chomp;
	@row = split;
	if ($prev ne "" and $prev ne $row[0]){
	    print OUT "$pline\n" if $count == 1;
	    $count =0;
	}
	$prev = $row[0];
	$pline = $_;
	$count++;
    }
    close(IN);
    close(OUT);
    system("rm ref20.$tag ref20.$tag.gz ref20_sort.$tag");
}

sub mk20mer{
    my $i;
    my $chr = shift;
    $file = "chr$chr";
    open(IN, $file);
    while(<IN>){
	$data = $_;
    }
    $length = length($data);
    for ($i = 0; $i <= $length - 20; $i++){
	$forward = substr($data, $i, 20);
	if ($forward !~ /[MRWSYKVHDBN]/){
	    $reverse = complement($forward);
	    $fpos = $i + 1;
	    $rpos = $fpos + 20 -1;
	    $tag = substr($forward, 0, 3);
	    print $tag "$forward\t$chr\t$fpos\tf\n";
	    $tag = substr($reverse, 0, 3);
	    print $tag "$reverse\t$chr\t$rpos\tr\n";
	}
    }
}

sub complement{
    my $seq = shift;
    my @seq = split('', $seq);
    my $length = length($seq);
    my $i = 0;
    my $out = "";
    for ($i = $length -1 ; $i >= 0; $i--){
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
