#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$target = $ARGV[0];
$ref   = $ARGV[1];
$opt   = $ARGV[2];

$wd = ".";

if ($target eq "" or $ref eq ""){
    print "
e.g. perl snp2vcf.pl SRR5886855 hg38

SRR5886855.snp.verify files in target directoy will be converted to SRR5886855.vcf.

" ;
    exit;
}

open(ALN, "$target/$target.aln");
open(OUT, "> $target/$target.vcf");
open(TMP, "|sort -S 1M -T $wd/$target > $wd/$target/$target.tmp");

print OUT "##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth)\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">
##INFO=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##INFO=<ID=SVLEN,Number=.,Type=Float,Description=\"Difference in length between REF and ALT alleles\">
##INFO=<ID=SVTYPE,Number=1,Type=Integer,Description=\"Type of structural variant\">
##ALT=<ID=DEL,Description=\"Deletion\">
##ALT=<ID=INS,Description=\"Insertion of novel sequence\">
##ALT=<ID=INV,Description=\"Inversion of reference sequence\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$target\n";


open(IN, "cat $target/$target.bi.snp|");
while(<IN>){
    next if $opt eq "" and ! /M|H/;
    chomp;
    @row = split('\t', $_);
    $dp = $row[$#row] -1;
    $chr= $row[0];
    $pos = $row[1];
    if ($chr =~ /\[0-9]/){
	$chr = "000$chr";
	$chr = substr($chr, length($chr) - 3, 3);
    }
    $pos = "00000000000" . $pos;
    $pos = substr($pos, length($pos) - 11, 11);
    
    $dp = $row[6] + $row[7];
    next if $dp == 0;
    $af = int(1000 * $row[7]/$dp)/1000;
    next if $row[7] <= 2; 
    if (/M/){
	print TMP "$chr\t$pos\t.\t$row[2]\t$row[3]\t1000\tPASS\tGT=1/1;AF=$af;DP=$dp\tGT:AD:DP\t1/1:$row[6],$row[7]:$dp\n";
    }elsif(/H/){
	print TMP "$chr\t$pos\t.\t$row[2]\t$row[3]\t1000\tPASS\tGT=0/1;AF=$af;DP=$dp\tGT:AD:DP\t0/1:$row[6],$row[7]:$dp\n";
    }else{
	$qual = $row[7] * 10;
	print TMP "$chr\t$pos\t.\t$row[2]\t$row[3]\t$qual\t.\tAF=$af;DP=$dp\tAD:DP\t$row[6],$row[7]:$dp\n";
    }
}


open(IN, "cat $target/$target.sv|");
while(<IN>){
    next if $opt eq "" and ! /M|H/;
    chomp;
    
    $info = "";
    
    @row = split('\t', $_);
    if ($prev_chr ne $row[0]){
	open(CHR, "$wd/$ref/chr$row[0]");
	binmode($fchr);
    }
    $prev_chr = $row[0];
    
    $dp = $row[$#row] -1;
    if (length($row[0]) <= 3 or $row[0] =~ /^[0-9]+$/){
	$chr =  $row[0];
    }
    
    if ($row[5] eq "deletion"){
	$pos = $row[1] - 1;
	seek(CHR, $pos -1, 0);
	read(CHR, $seq, abs($row[6]) + 1);
	$alt = substr($seq, 0, 1);
	$reference = $seq;
	($homseq = $row[12]) =~ y/\_//d;
	$homlen = length($homseq);
	$end = $row[3] + 1;
	if ($homlen > 0){
	    $info = "SVTYPE=DEL;END=$end;HOMLEN=$homlen;HOMSEQ=$homseq;SVLEN=$row[6];";
	}else{
	    $info = "SVTYPE=DEL;END=$end;SVLEN=$row[6];";
	}
	if ($row[6] <= -80){
	    $reference = $alt;
	    $alt = "<DEL>";
	}
    }elsif ($row[5] eq "insertion"){
	$alt = &searchInsertion($row[0], $row[1]);
	$pos = $row[1] - 1;
	seek(CHR, $pos -1, 0);
	read(CHR, $reference, 1);
	$end = $row[3] + 1;
	$info = "SVTYPE=INS;END=$end;SVLEN=$row[6];";	
	if ($alt =~/[0-9]/){
	    $alt = "<INS>";
	}
    }elsif ($row[5] eq "inversion"){
	$pos = $row[1] - 1;
	seek(CHR, $pos -1, 0);
	read(CHR, $reference, 1);
	$info = "SVTYPE=INV;END=$row[3];";	
	    $alt = "<INV>";
    }elsif ($row[5] eq "translocation"){
	$pos = $row[1] - 1;
	seek(CHR, $pos -1, 0);
	read(CHR, $reference, 1);
	$info = "SVTYPE=BND;";	
	$alt = $reference . "]$row[2]:$row[3]";
    }
    
    $dp = $row[9] + $row[10];
    next if $dp == 0;
    $af = int(1000 * $row[10]/$dp)/1000;
    next if $row[10] <= 2;
    if ($chr =~ /\[0-9]/){
	$chr = "000$chr";
	$chr = substr($chr, length($chr) - 3, 3);
    }
    $pos = "00000000000" . $pos;
    $pos = substr($pos, length($pos) - 11, 11);
    if (/M/){
	print TMP  "$chr\t$pos\t.\t$reference\t$alt\t1000\tPASS\t$info" . "GT=1/1;AF=$af;DP=$dp\tGT:AD:DP\t1/1:$row[9],$row[10]:$dp\n";
    }elsif(/H/){
	print TMP  "$chr\t$pos\t.\t$reference\t$alt\t1000\tPASS\t$info" . "GT=0/1;AF=$af;DP=$dp\tGT:AD:DP\t0/1:$row[9],$row[10]:$dp\n";
    }else{
	$qual = $row[10] * 10;
	print TMP  "$chr\t$pos\t.\t$reference\t$alt\t$qual\t.\t$info" . "AF=$af;DP=$dp\tAD:DP\t$row[9],$row[10]:$dp\n";
    }
    
}

close(TMP);
open(IN, "$target/$target.tmp");
while(<IN>){
    @row = split;
    $row[0] =~ s/^0*//;
    $row[1] =~ s/^0*//;
    $out = join("\t", @row);
    print OUT "$out\n";
}
close(IN);
close(OUT);
close(ALN);
system("rm $target/$target.tmp");

sub searchInsertion{
    my ($chr, $pos) = @_;
    my (@row, $size, $top, $bottom, $middle, $data, $ichr, $ipos, $hit);
    if ($chr =~ /^[0-9]*$/){
	$chr = "000$chr";
	$chr = substr($chr, length($chr) - 3, 3);
    }
    
    $pos = "00000000000" . $pos;
    $pos = substr($pos, length($pos) - 11, 11);
    
    if ($wd eq ""){
	$wd = ".";
    }
    
    $size = -s "$wd/$target/$target.index";
    open(INDEX, "$wd/$target/$target.index");
    binmode(INDEX);
    $top = 0;
    $bottom = $size;
    $middle = int($size / 2);
    while($bottom - $top > 1){
	seek(INDEX, $middle, 0);
	read(INDEX, $data, 1000);
	foreach (split('\n', $data)){
	    @row = split;
	    if ($row[1] =~ /^[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]$/){
		$ichr = $row[0];
		$ipos = $row[1];
		if ($chr eq $ichr){
		    if ($pos eq $ipos){
			$insert = &getInsert($row[3]);
			return $insert;
		    }
		}
	    }
	}
	if ($ichr gt $chr){
	    $bottom = $middle;
	}elsif($ichr eq $chr){
	    if ($ipos gt $pos){
		$bottom = $middle;
	    }else{
		$top = $middle;
	    }
	}else{
	    $top = $middle;
	}
	$middle = int(($top + $bottom) / 2);
    }

}

sub getInsert{
    my $address = shift;
    my $length = 0;
    my $count = 0;
    my $insert;
    seek(ALN, $address, 0);
    while(<ALN>){
	chomp;
	$count ++;
	if (/^#/){
	    $count = 1;
	    $insert_length = (split("\ ", $_))[7];
	}
	$pos = (split("\ ", $_))[2] if $count == 1;
	$length = length($_) if $count == 4;
	if ($count == 7){
	    $insert = substr($_, $length - 3, $insert_length + 1);
	    if (length($insert) < $insert_length){
		$insert = $insert_length;
	    }
	    return $insert;
	}
    }
}
