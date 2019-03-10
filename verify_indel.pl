#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     verify_indel.pl - verification program for structual valiations.

e.g. perl verify_indel.pl target control reference number type tmpdir

     perl verify_indel.pl ERR194147 default hg38 01 bi /mnt/ssd

     qsub -v target=ERR194147,control=default,ref=hg38,number=01,type=bi,tmpdir=/mnt/ssd verify_indel.pl

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
	if ($row[2] != 0){
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
if ($ftype eq "vcf"){
    open(IN, "$cwd/$target/$target.vcf");
}else{
    open(IN, "$cwd/$target/$target.indel.$number");
}
while(<IN>){
    chomp;
    if ($ftype eq "vcf"){
	($hchr, $hpos, $nop, $tref, $talt) = split('\t', $_);
	$hchr =~ s/chr//;
	$hchr += 0;
	$tref =~ y/a-z/A-Z/;
	$talt =~ y/a-z/A-Z/;
	$lref = length($tref);
	$lalt = length($talt);
	if($lref > 1 and $lalt == 1){
	    $head = substr($chr{$hchr}, $hpos - $length + 1, $length -1);
	    next if length($head) != $length -1;
	    $tail = substr($chr{$hchr}, $hpos + $lref -1, $length - 1);
	    next if length($tail) != $length - 1
	    &printTargetVcf($talt);
	}elsif($talt =~ /,/){
	    $head = substr($chr{$hchr}, $hpos - $length, $length -1);
	    next if length($head) != $length -1;
	    foreach (split(',', $talt)){
		$tail = $_ . substr($chr{$hchr}, $hpos + $lref -1, $length - 1);
		&printTargetVcf($_);
	    }
	}elsif($lref == 1 and $lalt >1){
	    $head = substr($chr{$hchr}, $hpos - $length, $length -1);
	    next if length($head) != $length -1;
	    $tail = $talt . substr($chr{$hchr}, $hpos, $length - 1);
	    &printTargetVcf($talt);
	}
    }else{
	@row = split('\t', $_);
	($hchr, $hpos, $tchr, $tpos, $direction, $type, $size) = (split(' ', $row[0]))[1.. 7];
	$current = "$hchr $hpos $tchr $tpos";
	next if $current eq $prev;
	$prev = $current;
	$posa = length($row[2]);
	$posb = length($row[8]);
	if ($posa < $posb){
	    $inside = substr($row[5], $posa -1, $posb - $posa + 1);
	    $head = substr($chr{$hchr}, $hpos - $length + ($posb - $posa -1), $length - ($posb - $posa -1) -1);
	    next if length($head) != $length - ($posb - $posa -1) -1;
	    if ($direction eq "f"){
		$tail = substr($chr{$tchr}, $tpos, $length - ($posb - $posa -1) -1);
		next if length($tail) != $length - ($posb - $posa -1) -1;
	    }else{
		$tail = &complement(substr($chr{$tchr}, $tpos - $length + ($posb - $posa -1), $length -($posb - $posa -1)-1));
		next if length($tail) != $length -($posb - $posa -1)-1;
	    }
	    $ref_seq = substr($chr{$hchr}, $hpos - $length, $length * 2 -1);
	    next if length($ref_seq) != $length * 2 -1;
	    $mut_seq = $head . $inside . $tail;
	}elsif($posa == $posb){
	    $inside = substr($row[5], $posa -1, $posb - $posa + 1);
	    $head = substr($chr{$hchr}, $hpos - $length + ($posb - $posa), $length - ($posb - $posa -1) -2);
	    next if length($head) != $length - ($posb - $posa -1) -2;
	    if ($direction eq "f"){
		$tail = substr($chr{$tchr}, $tpos, $length - ($posb - $posa) -1);
		next if length($tail) != $length - ($posb - $posa) -1;
	    }else{
		$tail = &complement(substr($chr{$tchr}, $tpos - $length + ($posb - $posa), $length -($posb - $posa)-1));
		next if length($tail) != $length -($posb - $posa)-1;
	    }
	    $ref_seq = substr($chr{$hchr}, $hpos - $length, $length * 2 -1);
	    next if length($ref_seq) != $length * 2 -1;
	    $mut_seq = $head . $inside . $tail;
	}else{
	    $inside = substr($row[5], $posb, $posa - $posb -1);
	    $head = substr($chr{$hchr}, $hpos - $length, $length - ($posa - $posb -1) -1);
	    next if length($head) != $length - ($posa - $posb -1) -1;
	    if ($direction eq "f"){
		$tail = substr($chr{$tchr}, $tpos + ($posa - $posb -1), $length - ($posa - $posb -1) -1);
	    }else{
		$tail = &complement(substr($chr{$tchr}, $tpos - $length, $length - ($posa - $posb -1) -1));
	    }
	    next if length($tail) != $length - ($posa - $posb -1) -1;
	    $ref_seq = substr($chr{$hchr}, $hpos - $length - ($posa - $posb), $length * 2 -1);
	    next if length($ref_seq) != $length * 2 -1;
	    $mut_seq = $head . $inside . $tail;
	}

	$slength = length($ref_seq) - $length;
	for($i = 0; $i <= $slength; $i++){
	    $tw = substr($ref_seq, $i, $length);
	    $tag = substr($tw, 0, 3);
	    print $tag "$tw\t$hchr $hpos $tchr $tpos $direction $type $size $inside tw\n" ;
	}

	$slength = length($mut_seq) - $length;
	for($i = 0; $i <= $slength; $i++){
	    $tm = substr($mut_seq, $i, $length);
	    $tag = substr($tm, 0, 3);
	    print $tag "$tm\t$hchr $hpos $tchr $tpos $direction $type $size $inside tm\n" ;
	}
    }
}
&closeTag;
&sortTag;

system("cat *.indel_sort.$number > indel_target.st.$number && rm *.indel_sort.$number");
if (-e "$cwd/$target/$target.sort_uniq.gz"){
    system("zcat $cwd/$target/$target.sort_uniq.gz 2> /dev/null |join - indel_target.st.$number | cut -d ' ' -f 2- > indel_target.$number && rm indel_target.st.$number");
}else{
    system("join $cwd/$target/$target.sort_uniq indel_target.st.$number | cut -d ' ' -f 2- > indel_target.$number && rm indel_target.st.$number");
}
system("sort -T $tmpdir indel_target.$number| uniq -c > indel_target.count.$number && rm indel_target.$number");

&openTag;

if($ftype eq "vcf"){
    open(IN,  "$cwd/$target/$target.vcf");
}else{
    open(IN,  "$cwd/$target/$target.indel.$number");
}
while(<IN>){
    chomp;
    if($ftype eq "vcf"){
	($hchr, $hpos, $nop, $tref, $talt) = split('\t', $_);
	$hchr =~ s/chr//;
	$hchr += 0;
	$tref =~ y/a-z/A-Z/;
	$talt =~ y/a-z/A-Z/;
	$lref = length($tref);
	$lalt = length($talt);
	if($lref > 1 and $lalt == 1){
	    $head = substr($chr{$hchr}, $hpos - $clength + 1, $clength -1);
	    $tail = substr($chr{$hchr}, $hpos + $lref -1, $clength - 1);
	    next if length($head) != $clength -1;
	    next if length($tail) != $clength -1;
	    &printControlVcf($talt);
	}elsif($lref == 1 and $lalt >1){
	    $head = substr($chr{$hchr}, $hpos - $clength, $clength -1);
	    $tail = $talt . substr($chr{$hchr}, $hpos, $clength - 1);
	    next if length($head) != $clength -1;
	    &printControlVcf($talt);
	}elsif($talt =~ /,/){
	    $head = substr($chr{$hchr}, $hpos - $clength, $clength -1);
	    next if length($head) != $clength -1;
	    foreach (split(',', $talt)){
		$tail = $_ . substr($chr{$hchr}, $hpos + $lref -1, $clength - 1);
	    &printControlVcf($_);
	    }
	}
	$ref_seq = substr($chr{$hchr}, $hpos - $clength, $clength * 2 -1);
	next if length($ref_seq) != $clength * 2 -1;
	$mut_seq = $head . $tail;
    }else{
	@row = split('\t', $_);
	($hchr, $hpos, $tchr, $tpos, $direction, $type, $size) = (split(' ', $row[0]))[1.. 7];
	$current = "$hchr $hpos $tchr $tpos";
	next if $current eq $prev;
	$prev = $current;
	$posa = length($row[2]);
	$posb = length($row[8]);
	
	if ($posa <= $posb){
	    $inside = substr($row[5], $posa -1, $posb - $posa + 1);
	    $head = substr($chr{$hchr}, $hpos - $clength + ($posb - $posa -1), $clength - ($posb - $posa -1) -1);
	    next if length($head) != $clength - ($posb - $posa -1) -1;
	    if ($direction eq "f"){
		$tail = substr($chr{$tchr}, $tpos, $clength - ($posb - $posa -1) -1);
	    }else{
		$tail = &complement(substr($chr{$tchr}, $tpos - $clength + ($posb - $posa -1), $clength -($posb - $posa -1)-1));
	    }
	    next if length($tail) != $clength -($posb - $posa -1)-1;
	    $ref_seq = substr($chr{$hchr}, $hpos - $clength, $clength * 2 -1);
	    next if length($ref_seq) != $clength * 2 -1;
	    $mut_seq = $head . $inside . $tail;
	}elsif($posa == $posb){
	    $inside = substr($row[5], $posa -1, $posb - $posa + 1);
	    $head = substr($chr{$hchr}, $hpos - $clength + ($posb - $posa), $clength - ($posb - $posa -1) -2);
	    next if length($head) != $clength - ($posb - $posa -1) -2;
	    if ($direction eq "f"){
		$tail = substr($chr{$tchr}, $tpos, $clength - ($posb - $posa) -1);
	    }else{
		$tail = &complement(substr($chr{$tchr}, $tpos - $clength + ($posb - $posa), $clength -($posb - $posa)-1));
	    }
	    next if length($tail) != $clength -($posb - $posa)-1;
	    $ref_seq = substr($chr{$hchr}, $hpos - $clength, $clength * 2 -1);
	    next if length($ref_seq) != $clength * 2 -1;
	    $mut_seq = $head . $inside . $tail;
	}else{
	    $inside = substr($row[5], $posb, $posa - $posb -1);
	    next if $clength - ($posa - $posb -1) -1 <= 0;
	    $head = substr($chr{$hchr}, $hpos - $clength, $clength - ($posa - $posb -1) -1);
	    next if length($head) != $clength - ($posa - $posb -1) -1;
	    if ($direction eq "f"){
		$tail = substr($chr{$tchr}, $tpos + ($posa - $posb -1), $clength - ($posa - $posb -1) -1);
	    }else{
		next if $clength - ($posa - $posb -1) -1 <= 0;
		$tail = &complement(substr($chr{$tchr}, $tpos - $clength, $clength - ($posa - $posb -1) -1));
	    }
	    next if length($tail) != $clength - ($posa - $posb -1) -1;
	    $ref_seq = substr($chr{$hchr}, $hpos - $clength - ($posa - $posb), $clength * 2 -1);
	    next if length($ref_seq) != $clength * 2 -1;
	    $mut_seq = $head . $inside . $tail;
	}
	$slength = length($ref_seq) - $clength;
	for($i = 0; $i <= $slength; $i++){
	    $cw = substr($ref_seq, $i, $clength);
	    $tag = substr($cw, 0, 3);
	    print $tag "$cw\t$hchr $hpos $tchr $tpos $direction $type $size $inside cw\n" ;
	}
	$slength = length($mut_seq) - $clength;
	for($i = 0; $i <= $slength; $i++){
	    $cm = substr($mut_seq, $i, $clength);
	    $tag = substr($cm, 0, 3);
	    print $tag "$cm\t$hchr $hpos $tchr $tpos $direction $type $size $inside cm\n" ;
	}
    }
}
&closeTag;
&sortTag;

system("cat *.indel_sort.$number | join $control - | cut -d ' ' -f 2- > indel_control.$number && rm *.indel_sort.$number");
system("sort -T $tmpdir indel_control.$number| uniq -c > indel_control.count.$number && rm indel_control.$number");

open(IN, "indel_target.count.$number");
while(<IN>){
    chomp;
    @row = split;
    $dat = join("\t", @row[1 .. $#row -1]);
    if ($ftype eq "vcf"){
	$dat .= "\t" if $#row == 4;
    }else{
	$dat .= "\t" if $#row == 8;
    }
    $all{"$row[1]\t$row[2]\t$row[3]\t$row[4]"} = $dat;
    if($row[$#row] eq "tw"){
        $tw{"$row[1]\t$row[2]\t$row[3]\t$row[4]"} = $row[0];
    }elsif($row[$#row] eq "tm"){
        $tm{"$row[1]\t$row[2]\t$row[3]\t$row[4]"} = $row[0];
    }
}
close(IN);

open(IN, "indel_control.count.$number");
while(<IN>){
    chomp;
    @row = split;
    $dat = join("\t", @row[1 .. $#row -1]);
    if ($ftype eq "vcf"){
	$dat .= "\t" if $#row == 4;
    }else{
	$dat .= "\t" if $#row == 8;
    }
    $all{"$row[1]\t$row[2]\t$row[3]\t$row[4]"} = $dat;
    if($row[$#row] eq "cw"){
        $cw{"$row[1]\t$row[2]\t$row[3]\t$row[4]"} = $row[0];
    }elsif($row[$#row] eq "cm"){
        $cm{"$row[1]\t$row[2]\t$row[3]\t$row[4]"} = $row[0];
    }
}
close(IN);

system("rm indel_target.count.$number indel_control.count.$number");

open (OUT, "> $target.indel_result.$number");
foreach $dat (sort keys %all){
    $cw = $cw{$dat};
    $cm = $cm{$dat};
    $tw = $tw{$dat};
    $tm = $tm{$dat};
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
    @row = split('\t', $all{$dat});
    if ($ftype eq "vcf"){
	print OUT "$row[0]\t$row[1]\t$row[2]\t$row[3]\t$cw\t$cm\t$tw\t$tm\t$genotype\t$row[4]\n";
    }else{
	if ($row[5] =~ /inv|trans/){
	    $row[7] = $row[6];
	    $row[6] = "";
	}
	print OUT "$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$cw\t$cm\t$tw\t$tm\t$genotype\t$row[7]\n";
    }
    $chr_exist{$row[0]} = 1;
    $chr_exist{$row[2]} = 1;
}

if ( -e "$cwd/$target/$target.indel.verify.$number"){
    system("rm $cwd/$target/$target.indel.verify.$number");
}

foreach $chr (@chr){
    next if ! $chr_exist{$chr};
    open(IN, "$target.indel_result.$number");
    if ($ftype eq "vcf"){
	open(OUT, "|sort -T $tmpdir -k 2 -n >> $cwd/$target/$target.indel_vcf.verify.$number");
    }else{
	open(OUT, "|sort -T $tmpdir -k 2 -n >> $cwd/$target/$target.indel.verify.$number");
    }
    while(<IN>){
	@row = split;
	if ($row[0] eq $chr){
	    print OUT;
	}
    }
    close(IN);
    close(OUT);
}

system("rm $target.indel_result.$number");

if ($tmpdir ne "."){
    system("rm -r $workdir");
}

sub openTag{
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "> $workdir/$tag.indel_tmp.$number");
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
		system("sort -T $tmpdir $tag.indel_tmp.$number > $tag.indel_sort.$number && rm $tag.indel_tmp.$number");
	    }
	}
    }
}

sub printTargetVcf{
    my $talt = shift;
    $talt =~ y/a-z/A-Z/;
    $ref_seq = substr($chr{$hchr}, $hpos - $length, $length * 2 -1);
    $mut_seq = $head . $tail;
    $slength = length($ref_seq) - $length;
    for($i = 0; $i <= $slength; $i++){
	$tw = substr($ref_seq, $i, $length);
	$tag = substr($tw, 0, 3);
	print $tag "$tw\t$hchr $hpos $tref $talt tw\n" ;
    }
    $slength = length($mut_seq) - $length;
    for($i = 0; $i <= $slength; $i++){
	$tm = substr($mut_seq, $i, $length);
	$tag = substr($tm, 0, 3);
	print $tag "$tm\t$hchr $hpos $tref $talt tm\n" ;
    }
}

sub printControlVcf{
    my $talt = shift;
    $talt =~ y/a-z/A-Z/;
    $ref_seq = substr($chr{$hchr}, $hpos - $clength, $clength * 2 -1);
    $mut_seq = $head . $tail;
    $slength = length($ref_seq) - $clength;
    for($i = 0; $i <= $slength; $i++){
	$cw = substr($ref_seq, $i, $clength);
	$tag = substr($cw, 0, 3);
	print $tag "$cw\t$hchr $hpos $tref $talt cw\n" ;
    }
    $slength = length($mut_seq) - $clength;
    for($i = 0; $i <= $slength; $i++){
	$cm = substr($mut_seq, $i, $clength);
	$tag = substr($cm, 0, 3);
	print $tag "$cm\t$hchr $hpos $tref $talt cm\n" ;
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
