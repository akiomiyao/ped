#!/usr/bin/perl
#
# This file is a script for Polymorphic Edge Detection.
#
# Multithreaded version of all in one.
#
# Copyright (C) 2019 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#
# Author: Akio Miyao <miyao@affrc.go.jp>
#

use threads;
use Thread::Semaphore;

my $max_semaphore = 14;

$usage = '
 ped.pl - program for polymorphic edge detection. 

 e.g. perl ped.pl target=ERR194147,ref=hg38

 If you want to detect polymorphisms between target and control,
 e.g. perl ped.pl target=ERR3063487,control=ERR3063486,ref=WBcel235

 If you want to specify the working directory,
 e.g. perl ped.pl target=ERR194147,ref=hg38,wd=/home/you/work

 For setting of muximum number of threads,
 e.g. perl ped.pl target=ERR194147,ref=hg38,thread=16

 Results will be saved in target directory.

     Author: Akio Miyao <miyao@affrc.go.jp>

';

my $sort_opt = "-S 1M";
@nuc = ('A', 'C', 'G', 'T');
$ENV{LANG} = "C";

if ($ARGV[0] =~ /target|ref/){
    my @arg = split(',', $ARGV[0]);
    foreach (sort @arg){
	my ($name, $val) = split('=', $_);
	$$name = $val;
	$log .= "# $name: $val\n";
    }
}elsif ($ARGV[1] ne ""){
    $target = $ARGV[0];
    $control = $ARGV[1];
    $ref = $ARGV[2];
    $method = $ARGV[3];
}else{
    print $usage;
    exit;
}

if ($wd eq ""){
    $wd = `pwd`;
    chomp($wd);
}

$uname = `uname`;
chomp($uname);

open(REPORT, "> $wd/$target/$target.log");
print REPORT "# Log of ped.pl
$log";
report("Job begin");

if ($thread ne ""){
    $max_semaphore = $thread;
}else{
    if ($uname eq "FreeBSD"){
	open(IN, "sysctl kern.smp.cpus |");
	while(<IN>){
	    chomp;
	    $processor = (split(': ', $_))[1];
	}
	close(IN);
    }elsif($uname eq "Linux"){
	open(IN, "/proc/cpuinfo");
	while(<IN>){
	    $processor ++  if /processor/;
	}
	close(IN);
    }
    if ($processor > 6){
	$processor -= 2;
	$max_semaphore = $processor;
	$semaphore4sort = 4;
    }
}
my $semaphore = Thread::Semaphore->new($max_semaphore);

$refdir = "$wd/$ref";
system("mkdir $wd/$target/tmp") if ! -e "$wd/$target/tmp";
$tmpdir = "$wd/$target/tmp";

$control = $ref if $control eq "" or $control eq "default";
open(IN, "$wd/config");
while(<IN>){
    chomp;
    @row = split('\t', $_);
    if ($row[1] eq "description"){
	$desc{$row[0]} = $row[2];
    }elsif($row[1] eq "wget"){
	$wget{$row[0]} = $row[2];
    }elsif($row[1] eq "file"){
	$file{$row[0]} = $row[2];
    }elsif($row[0] eq $ref && $row[1] eq "chromosome"){
	@row = split;
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

if (! -e "$wd/$target/sort_uniq/$target.sort_uniq.TTT.gz"){
    &mkSortUniq($target);
}

if ($control ne $ref && ! -e "$wd/$control/sort_uniq/$control.sort_uniq.TTT.gz"){
    &mkSortUniq($control);
}

if (! -e "$wd/$ref/sort_uniq/$ref.sort_uniq.TTT.gz" and ! -e "$wd/$ref/sort_uniq/$ref.TTT.gz"){
    &mkRef;
}

open(IN, "zcat $wd/$target/sort_uniq/*.gz 2> /dev/null |");
while(<IN>){
    chomp;
    $length = length($_);
    last;
}
close(IN);

open(IN, "zcat $wd/$control/sort_uniq/*.gz 2> /dev/null |");
while(<IN>){
    chomp;
    $clength = length($_);
    last;
}
close(IN);

&mkData4MapF;

$semaphore->down($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
&sortData4Map;

$semaphore->up($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
&map;

&mkData4MapR;

$semaphore->down($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
&sortData4Map;

$semaphore->up($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
&map;
&align;

system("cat $tmpdir/$target.aln.* > $wd/$target/$target.aln && rm $wd/$target/tmp/*");

$semaphore->down;
threads->new(\&index);

$semaphore->down($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
threads->new(\&svSort);
&joinAll;

$semaphore->down;
threads->new(\&svMkT);
&joinAll;

&sortSeq;

$semaphore->up($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
&joinTarget;

$semaphore->down;
threads->new(\&svMkC);
&joinAll;

$semaphore->down($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
&sortSeq;

$semaphore->up($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
&joinControl;

&svReadCount;

$semaphore->down;
threads->new(\&snpMkT);
&joinAll;

$semaphore->down($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
&sortSeq;

$semaphore->up($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
&joinTarget;

$semaphore->down;
threads->new(\&snpMkC);
&joinAll;

$semaphore->down($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
&sortSeq;

$semaphore->up($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
&joinControl;

&snpReadCount;

&toVcf;
system("rm -r $wd/$target/tmp");
report("Job completed");

sub mkRef{
    my (@row, $i, $remote_file);
    if (! -e "$wd/$ref"){
	system("mkdir $wd/$ref");
    }
    chdir "$wd/$ref";
    if ($file eq ""){
	@row = split('/', $wget{$target});
	$remote_file = $row[$#row];
	if (! -e $remote_file and $wget{$ref} ne ""){
	    &report("Downloading $wget{$ref}");
	    system("wget $wget{$ref}");
	}
	&mkChr;
    }else{
	&mkChrFromFile($file);
    }
    &mk20;
    &mkControlRead;
    system("rm -r $wd/$ref/tmp");
}

sub mkChrFromFile{
    &report("Making chromosome file.");
    my $out;
    die "$file is not found in $wd/$ref." if ! -e "$wd/$ref/$file";
    @chr = ();
    if ($file =~ /gz$/){
	open(IN, "zcat $wd/$ref/$file|");
    }elsif($file =~ /bz2$/){
	open(IN, "bzcat $wd/$ref/$file|");
    }else{
	open(IN, "$wd/$ref/$file");
    }
    while(<IN>){
	chomp;
	if (/^>/){
	    close(OUT);
	    $out = (split)[0];
	    $out =~ s/>//;
	    $out =~ s/^chr//i;
	    $out += 0 if $out =~ /^[0-9]*$/;
	    push(@chr, $out) if ! $chr_flag;
	    $out = "chr$out";
	    open(OUT, "> $out");
	}else{
	    y/a-z/A-Z/;
	    print OUT;
	}
    }
    close(IN);
    close(OUT);
}

sub mkControlRead{
    &report("Making Control Read.");
    if (! -e "$wd/$ref/sort_uniq"){
	system("mkdir $wd/$ref/sort_uniq");
    }

    for (@chr){
	next if $_ eq "NOP";
	$semaphore->down;
	&report("Processing Chr$_");
	threads->new(\&mkControlReadChr, $_);
    }
    &joinAll;

    $semaphore->down($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		&report("Making Control Read. $tag");
		threads->new(\&mkControlReadSub, $tag);
	    }
	}
    }
    &joinAll;
    $semaphore->up($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
}

sub mkControlReadChr{
    my $chr = shift;
    my ($tag, $nuca, $nucb, $nucc, $fin, $fout, $read, $data, $i);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$fout = "$tag.$chr";
		open($fout, "> $wd/$ref/tmp/$tag.$chr");
	    }
	}
    }
    $fin = $chr;
    open($fin, "$wd/$ref/chr$chr");
    $data = <$fin>;
    close(IN);
    $i = 0;
    while(1){
	$read = substr($data, $i, 100);
	last if length($read) != 100;
	if ($read !~ /[MRWSYKVHDBN]/){
	    $tag = substr($read, 0, 3) . ".$chr";
	    
	    print $tag "$read\n";
	    $read = &complement($read);
	    $tag = substr($read, 0, 3) . ".$chr";
		print $tag "$read\n";
	}
	$i += 2;
    }
    close($fin);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		close($fout);
	    }
	}
    }
   
    $semaphore->up;
}

sub mkControlReadSub{
    my $tag = shift;
    system("cat $wd/$ref/tmp/$tag.* | sort -T $wd/$ref/tmp/ $sort_opt  | uniq | gzip > $wd/$ref/sort_uniq/$ref.sort_uniq.$tag.gz && rm $wd/$ref/tmp/$tag.*");
    $semaphore->up;
}

sub mkUniq{
    my $tag = shift;
    &report("Making ref20_uniq.$tag.gz");
    system("cat $wd/$ref/tmp/ref20.$tag.* | sort -T $tmpdir $sort_opt  > $wd/$ref/tmp/ref20_sort.$tag");
    &waitFile("$wd/$ref/tmp/ref20_sort.$tag");
    open(OUT, "|gzip -f > $wd/$ref/ref20_uniq.$tag.gz");
    open(IN, "$wd/$ref/tmp/ref20_sort.$tag");
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
    &waitFile("$wd/$ref/ref20_uniq.$tag.gz");
    system("rm $wd/$ref/tmp/ref20.$tag.* $wd/$ref/tmp/ref20_sort.$tag");
    $semaphore->up;
}

sub mk20mer{
    my ($fin, $fout, $i, @tag, $tag, $data, $forward, $length, $comp, $fpos, $rpos, );
    my $chr = shift;

    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		$fout = "$tag.$chr";
		open($fout, "> $wd/$ref/tmp/ref20.$tag.$chr");
	    }
	}
    }

    my $file = "chr$chr";
    open(IN, $file);
    while(<IN>){
	$data = $_;
    }
    $length = length($data);
    for ($i = 0; $i <= $length - 20; $i++){
	$forward = substr($data, $i, 20);
	if ($forward !~ /[MRWSYKVHDBN]/){
	    $comp = complement($forward);
	    $fpos = $i + 1;
	    $rpos = $fpos + 20 -1;
	    $tag = substr($forward, 0, 3) . ".$chr";
	    print $tag "$forward\t$chr\t$fpos\tf\n";
	    $tag = substr($comp, 0, 3) . ".$chr";
	    print $tag "$comp\t$chr\t$rpos\tr\n";
	}
    }
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag) . ".$chr";
		close($tag);
	    }
	}
    }
    $semaphore->up;
}

sub mk20{
    system("mkdir $wd/$ref/tmp");
    &report("Making 20mer position file.");
    
    foreach $i (@chr){
	next if $i eq "NOP";
	$semaphore->down;
	&report("Processing Chr$i");
	threads->new(\&mk20mer, $i);
    }
    &joinAll;
    
    $semaphore->down($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		$semaphore->down;
		threads->new(\&mkUniq, $tag);
	    }
	}
    }    
    &joinAll;
    $semaphore->up($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
}

sub mkChr{
    my $output;
    my @file = split('/', $wget{$ref});
    my $file = $file[$#file];
    $file =~ s/\ +$//g;
    my $i = 0;
    &report("Making chromosome file from $file.");
    if ($ref eq "hg38"){
	open(IN, "zcat $file|");
	while(<IN>){
	    chomp;
	    if (/^>/){
		close(OUT);
		$output = (split)[0];
		$output =~ s/>//;
		if ($output =~ /\_/){
		    $flag = 1;
		}else{
		    $flag = 0;
		    $output =~ s/chr//i;
		    $output += 0 if $output =~ /^[0-9]*$/;
		    $output = "chr$output";
		    open(OUT, "> $output");
		}
	    }elsif(! $flag){
		y/a-z/A-Z/;
		print OUT;
	    }
	}
	close(IN);
	close(OUT);
    }else{
	if ($ref eq "IWGSC1.0"){
	    system("unzip iwgsc_refseqv1.0_all_chromosomes.zip && mv iwgsc_refseqv1.0_all_chromosomes/$file{$ref} .");
	}elsif($ref eq "SL3"){
	    system("tar xvfz $file && rm $file");
	    $file =~ s/\.tar\.gz//;
 	}elsif($ref =~ /^B73/){
	    for ($i = 1; $i <= 10; $i++){
		$file = $wget{$ref} . "chr$i.fna.gz";
		&report("Downloading $file");
		system("wget $file");
		$file = "chr$i.fna.gz";
		open(IN, "zcat $file|");
		open(OUT, "> chr$i");
		while(<IN>){
		    chomp;
		    if (! /^>/){
			y/a-z/A-Z/;
			print OUT;
		    }
		}
		close(IN);
		close(OUT);
	    }
	    return;
	}
	
	if ($file =~ /gz$/){
	    open(IN, "zcat $file|");
	}elsif ($file =~ /bz2$/){
	    open(IN, "bzcat $file|");
	}elsif ($file =~ /xz$/){
	    open(IN, "xzcat $file|");
	}else{
	    $file = $file{$ref} if $file{$ref} ne "";
	    open(IN, $file);
	}
	while(<IN>){
	    chomp;
	    if (/^>/){
		close(OUT);
		$output = "chr" . $chr[$i];
		last if $chr[$i] eq "";
		open(OUT, "> $output");
		$i++;
	    }else{
		y/a-z/A-Z/;
		print OUT;
	    }
	}
	close(IN);
	close(OUT);
    }
}

sub mkSortUniq{
    my $target = shift;
    report("Making sort_uniq files for $target");
    opendir(DIR, "$wd/$target/read");
    foreach (sort readdir(DIR)){
	if (/gz$/){
	    $gz_file .= "$_ ";
	}elsif(/bz2$/){
	    $bz_file .= "$_ ";
	}elsif(/xz$/){
	    $xz_file .= "$_ ";
	}elsif(/q$/){
	    $fq_file .= "$_ ";
	}
    }
    close(DIR);
    system("mkdir $wd/$target/sort_uniq") if ! -e "$wd/$target/sort_uniq";
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "> $wd/$target/sort_uniq/$tag.tmp");
	    }
	}
    }
    if ($gz_file ne ""){
	$cmd = "cd $wd/$target/read && zcat $gz_file |";
	&sortUniqSub($cmd);
    }
    if($bz_file ne ""){
	$cmd = "cd $wd/$target/read && bzcat $gb_file |";
	&sortUniqSub($cmd);
    }
    if($xz_file ne ""){
	$cmd = "cd $wd/$target/read && xzcat $gb_file |";
	&sortUniqSub($cmd);
    }
    if($fq_file ne ""){
	$cmd = "cd $wd/$target/read && cat $fq_file |";
	&sortUniqSub($cmd);
    }
    &closeTag;
    $semaphore->down($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		threads->new(\&sortUniqSort, $tag, $target);
	    }
	}
    }
    &joinAll;
    $semaphore->up($max_semaphore - $semaphore4sort) if $semaphore4sort != 0;
}

sub sortUniqSort{
    my ($tag, $target) = @_;
    report("Making $target.sort_uniq files: Sorting for $target.sort_uniq.$tag.gz");
    system("sort -T $wd/$target/sort_uniq/ $sort_opt $wd/$target/sort_uniq/$tag.tmp | uniq | gzip > $wd/$target/sort_uniq/$target.sort_uniq.$tag.gz && rm $wd/$target/sort_uniq/$tag.tmp");
    $semaphore->up;
}

sub sortUniqSub{
    my $cmd = shift;
    my $count = 0;
    open(IN, $cmd);
    while(<IN>){
	if ($count == 1 and !/N/){
	    chomp;
	    $tag = substr($_, 0, 3);
	    print $tag "$_\n";
	    $complement = &complement($_);
	    $tag = substr($complement, 0, 3);
	    print $tag "$complement\n";
	    $total ++;
	    if ($total % 1000000 == 0){
		report("Making $target.sort_uniq files: Split to subfiles. $total reads processed");
	    }
	}elsif($count == 4){
	    $count = 0;
	}
	$count++;
    }
    close(IN);
}

sub joinAll{
    foreach my $thr (threads->list){
	$thr->join;
    }
}

sub toVcf{
    report("Convert to vcf format");
    my $fin = "in-vcf";
    my $fout ="out-vcf";
    open($fin, "cat $wd/$target/$target.bi.snp|");
    open($fout, "> $wd/$target/$target.bi.vcf");
    print $fout "##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth)\">
##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$target\n";
    while(<$fin>){
	chomp;
	@row = split('\t', $_);
	$dp = $row[$#row] -1;
	if (length($row[0]) <= 3 or $row[0] =~ /^[0-9]+$/){
	    $row[0] = "chr" . $row[0];
	}
	$dp = $row[7] + $row[8];
	next if $dp == 0;
	$af = int(1000 * $row[8]/$dp)/1000;
	if (/M/){
	    print $fout "$row[0]\t$row[1]\t.\t$row[2]\t$row[3]\t1000\tPASS\tGT=1/1;AF=$af;DP=$dp\tGT:AD:DP\t1/1:$row[7],$row[8]:$dp\n";
	}elsif(/H/){
	    print $fout "$row[0]\t$row[1]\t.\t$row[2]\t$row[3]\t1000\tPASS\tGT=0/1;AF=$af;DP=$dp\tGT:AD:DP\t0/1:$row[7],$row[8]:$dp\n";
	}
    }
}

sub snpMkT{
    report("Making data for verification of snp. target");
    my (@row, $chr, $pos, $ref, $alt, $count, $ref_seq, $mut_seq, $head, $tail, $i, $ipos, $iref, $ialt, $tpos, $tw, $tm, $tag, $nuca, $nucb, $nucc, @dat);
    my $s = {};

    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "> $tmpdir/$tag.tmp");
	    }
	}
    }
    my $fin = "in-mkt";
    my $fout = "out-mkt";
    open($fout, "|sort $sort_opt -T $tmpdir | uniq -c > $wd/$target/tmp/snp.list");
    open($fin, "$wd/$target/$target.aln");
    while(<$fin>){
	chomp;
	if (/^#/ and /snp/){
	    @row = split;
	    if ($row[1] !~/[IVXYZ]/){
		$row[1] = "000" . $row[1];
		$row[1] = substr($row[1], length($row[1]) -3, 3);
	    }
	    $row[2] = "000000000000" . $row[2];
	    $row[2] = substr($row[2], length($row[2]) -11, 11);
	    print $fout "$row[1] $row[2] $row[4] $row[5]\n";
	}
    }
    close($fin);
    close($fout);

    open($fin, "$wd/$target/tmp/snp.list");
    while(<$fin>){
	chomp;
	($count, $chr, $pos, $ref, $alt) = split;
	$chr += 0 if $chr !~ /[IVXYZ]/;
	if ($prev_chr ne $chr){
	    my $chr_file = "$refdir/chr$chr";
	    open (IN, $chr_file);
	    ($chr_seq = <IN>) =~ y/a-z/A-Z/;
	    close(IN);
	    @dat = ();
	    report("Making data for verification of snp. target chr$chr");
	}
	$count = int(1 + $count/4);
	$pos += 0;
	if ($count >= 5){
	    push(@dat, "$pos $ref $alt");
	    if ($#dat >= 10){
		shift(@dat);
	    }
	}
	$ref_seq = substr($chr_seq, $pos - $length, $length * 2 -1);
	next if length($ref_seq) != $length * 2 -1;
	$head = substr($ref_seq, 0, $length-1);
	$tail = substr($ref_seq, $length, $length);
	$mut_seq = $head . $alt . $tail;
	foreach (@dat){
	    ($tpos, $iref, $ialt) = split;
	    $i = $tpos - ($pos - $length) -1;
	    if ($i > 0 and $i < $length * 2){
		$head = substr($mut_seq, 0, $i);
		$tail = substr($mut_seq, $i + 1);
		$mut_seq = $head . $ialt . $tail;
	    }
	}
	for($i = 0; $i < $length; $i++){
	    $tw = substr($ref_seq, $i, $length);
	    $tm = substr($mut_seq, $i, $length);
	    $tag = substr($tw, 0, 3);
	    print $tag "$tw\t$chr $pos $ref $alt tw\n";
	    $tag = substr($tm, 0, 3);
	    print $tag "$tm\t$chr $pos $ref $alt tm\n";
	}
	$prev_chr = $chr;
    }
    close($fin);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		close($tag);
	    }
	}
    }
    $semaphore->up;
}

sub snpMkC{
    my (@row, $chr, $pos, $ref, $alt, $count, $chr_seq, $ref_seq, $mut_seq, $head, $tail, $i, $ipos, $iref, $ialt, $tpos, $cw, $cm, $tag, $nuca, $nucb, $nucc, @dat);
    my $s = {};

    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "> $tmpdir/$tag.tmp");
	    }
	}
    }
    my $fin = "in-mkc";

    open($fin, "$wd/$target/tmp/snp.list");
    while(<$fin>){
	chomp;
	($count, $chr, $pos, $ref, $alt) = split;
	$chr += 0 if $chr !~ /[IVXYZ]/;
	if ($prev_chr ne $chr){
	    my $chr_file = "$refdir/chr$chr";
	    open (IN, $chr_file);
	    ($chr_seq = <IN>) =~ y/a-z/A-Z/;
	    close(IN);
	    @dat = ();
	    report("Making data for verification of snp. control chr$chr");
	}
	$count = int(1 + $count/4);
	$pos += 0;
	if ($count >= 5){
	    push(@dat, "$pos $ref $alt");
	    if ($#dat >= 10){
		shift(@dat);
	    }
	}
	$ref_seq = substr($chr_seq, $pos - $clength, $clength * 2 -1);
	next if length($ref_seq) != $clength * 2 -1;
	$head = substr($ref_seq, 0, $clength-1);
	$tail = substr($ref_seq, $clength, $clength);
	$mut_seq = $head . $alt . $tail;

	foreach (@dat){
	    ($tpos, $iref, $ialt) = split;
	    $i = $tpos - ($pos - $clength) -1;
	    if ($i > 0 and $i < $clength * 2){
		$head = substr($mut_seq, 0, $i);
		$tail = substr($mut_seq, $i + 1);
		$mut_seq = $head . $ialt . $tail;
	    }
	}
	for($i = 0; $i < $clength; $i++){
	    $cw = substr($ref_seq, $i, $clength);
	    $cm = substr($mut_seq, $i, $clength);
	    $tag = substr($cw, 0, 3);
	    print $tag "$cw\t$chr $pos $ref $alt cw\n";
	    $tag = substr($cm, 0, 3);
	    print $tag "$cm\t$chr $pos $ref $alt cm\n";
	}
	$prev_chr = $chr;
    }
    close($fin);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		close($tag);
	    }
	}
    }
    $semaphore->up;
}

sub svReadCount{
    my (@row, $dat, %count, $cm, $cw, $tm, $tw, $count, @prev, $prev, $genotype);
    open(OUT, "> $wd/$target/$target.sv");
    open(IN, "cat $tmpdir/*.target $tmpdir/*.control |");
    while(<IN>){
	chomp;
	@row = split;
	if ($row[0] !~/[IVXYZ]/){
	    $row[0] = "000" . $row[0];
	    $row[0] = substr($row[0], length($row[0]) -3, 3);
	}
	$row[1] = "000000000000" . $row[1];
	$row[1] = substr($row[1], length($row[1]) -11, 11);
	$dat = join(" ", @row);
	$count{$dat} ++;
    }
    close(IN);
    $cm = 0;
    $cw = 0;
    $tm = 0;
    $tw = 0;
    foreach (sort keys %count){
	$count = $count{$_};
	@row = split;
	if($row[0] !~/[IVXYZ]/){
	    $row[0] += 0;
	}
	$row[1] += 0;
	$dat = join("\t", @row[0 .. $#row -1]);
	if($dat ne $prev){
	    $genotype = "";
	    if ($cw >= 5 and $cm <= 1){
		if ($tm >= 5 and $tw <= 1){
		    $genotype = 'M';
		}elsif(($tm + $tw) > 0 and $tm / ($tm + $tw) > 0.3 and $tm / ($tm + $tw) < 0.7 and $cm <= 1 and $tw >= 5){
		    $genotype = "H";
		}
	    }
	    @prev = split('\t', $prev);
	    if ($prev[5] =~~/version/){
		$prev[7] = $prev[6];
		$prev[6] = "_";
	    }
	    print OUT "$prev[0]\t$prev[1]\t$prev[2]\t$prev[3]\t$prev[4]\t$prev[5]\t$prev[6]\t$cw\t$cm\t$tw\t$tm\t$genotype\t$prev[7]\n" if $prev ne "";
	    $cm = 0;
	    $cw = 0;
	    $tm = 0;
	    $tw = 0;
	}	
	$cm = $count if $row[$#row] eq "cm";
	$cw = $count if $row[$#row] eq "cw";
	$tm = $count if $row[$#row] eq "tm";
	$tw = $count if $row[$#row] eq "tw";
	$prev = $dat;
    }
}

sub snpReadCount{
    my (@row, $dat, %count, $cm, $cw, $tm, $tw, $count, @prev, $prev, $genotype);
    open(OUT, "> $wd/$target/$target.bi.snp");
    open(IN, "cat $tmpdir/*.target $tmpdir/*.control |");
    while(<IN>){
	chomp;
	@row = split;
	if ($row[0] !~/[IVXYZ]/){
	    $row[0] = "000" . $row[0];
	    $row[0] = substr($row[0], length($row[0]) -3, 3);
	}
	$row[1] = "000000000000" . $row[1];
	$row[1] = substr($row[1], length($row[1]) -11, 11);
	$dat = join(" ", @row);
	$count{$dat} ++;
    }
    close(IN);
    $cm = 0;
    $cw = 0;
    $tm = 0;
    $tw = 0;
    foreach (sort keys %count){
	$count = $count{$_};
	@row = split;
	if($row[0] !~/[IVXYZ]/){
	    $row[0] += 0;
	}
	$row[1] += 0;
	$dat = join("\t", @row[0 .. $#row -1]);
	if($dat ne $prev){
	    $genotype = "";
	    if ($cw >= 5 and $cm <= 1){
		if ($tm >= 5 and $tw <= 1){
		    $genotype = 'M';
		}elsif(($tm + $tw) > 0 and $tm / ($tm + $tw) > 0.3 and $tm / ($tm + $tw) < 0.7 and $cm <= 1 and $tw >= 5){
		    $genotype = "H";
		}
	    }
	    print OUT "$prev\t$cw\t$cm\t$tw\t$tm\t$genotype\n" if $prev ne "";
	    $cm = 0;
	    $cw = 0;
	    $tm = 0;
	    $tw = 0;
	}	
	$cm = $count if $row[$#row] eq "cm";
	$cw = $count if $row[$#row] eq "cw";
	$tm = $count if $row[$#row] eq "tm";
	$tw = $count if $row[$#row] eq "tw";
	$prev = $dat;
    }
}

sub joinControlFunc{
    my $tag = shift;
    report("Selecting sequence data for verify. $tag");
    if (-e "$wd/$control/sort_uniq/$control.$tag.gz"){
	system("zcat $wd/$control/sort_uniq/$control.$tag.gz 2> /dev/null | join - $tmpdir/$tag | cut -d ' ' -f 2- > $tmpdir/$tag.control");
    }else{
	system("zcat $wd/$control/sort_uniq/$control.sort_uniq.$tag.gz 2> /dev/null | join - $tmpdir/$tag | cut -d ' ' -f 2- > $tmpdir/$tag.control");
    }
    &waitFile("$tmpdir/$tag.control");
    system("rm $tmpdir/$tag");
    $semaphore->up;
}

sub joinControl{
    my ($nuca, $nucb, $nucc, $tag);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		threads->new(\&joinControlFunc, $tag);
	    }
	}
    }
    foreach my $thr (threads->list){
	$thr->join;
    }
}

sub joinTargetFunc{
    my $tag = shift;
    report("Selecting sequence data for verify. $tag");
    system("zcat $wd/$target/sort_uniq/$target.sort_uniq.$tag.gz 2> /dev/null | join - $tmpdir/$tag | cut -d ' ' -f 2- > $tmpdir/$tag.target");
    &waitFile("$tmpdir/$tag.target");
    system("rm $tmpdir/$tag");
    $semaphore->up;
}

sub joinTarget{
    my ($nuca, $nucb, $nucc, $tag);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		threads->new(\&joinTargetFunc, $tag);
	    }
	}
    }
    foreach my $thr (threads->list){
	$thr->join;
    }
}

sub sortSeqFunc{
    my $tag = shift;
    report("Sorting sequence data for verify. $tag");
    system("cat $tmpdir/$tag.tmp | sort $sort_opt -T $tmpdir  > $tmpdir/$tag");
    &waitFile("$tmpdir/$tag");
    system("rm $tmpdir/$tag.tmp");
    $semaphore->up;
}

sub sortSeq{
    my ($nuca, $nucb, $nucc, $tag);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		threads->new(\&sortSeqFunc, $tag);
	    }
	}
    }
    foreach my $thr (threads->list){
	$thr->join;
    }
}

sub svMkC{
    report("Making data for verification of sv");
    my ($nuca, $nucb, $nucc, $tag, $hchr, $hpos, $tchr, $tpos, $direction, $type, $size, @row, $current, $prev, $posa, $posb, $inside, $head, $tail, $ref_seq, $ut_seq, $slength, $tw);
    foreach $i (@chr){
	my $chr_file = "$refdir/chr$i";
	open (IN, $chr_file);
	($chr{$i} = <IN>) =~ y/a-z/A-Z/;
	close(IN);
    }

    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "> $tmpdir/$tag.tmp");
	    }
	}
    }
    my $fin = "in-sva";
    open($fin, "$tmpdir/$target.sv.sort");
    while(<$fin>){
	chomp;
	@row = split('\t', $_);
	($hchr, $hpos, $tchr, $tpos, $direction, $type, $size) = (split(' ', $row[0]))[1.. 7];
	$hchr =~ s/^0+//;
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
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		close($tag);
	    }
	}
    }
    $semaphore->up;
    report("Making data for verification of sv complete");
}

sub svMkT{
    report("Making data for verification of sv");
    my ($nuca, $nucb, $nucc, $tag, $hchr, $hpos, $tchr, $tpos, $direction, $type, $size, @row, $current, $prev, $posa, $posb, $inside, $head, $tail, $ref_seq, $ut_seq, $slength, $tw);
    foreach $i (@chr){
	my $chr_file = "$refdir/chr$i";
	open (IN, $chr_file);
	($chr{$i} = <IN>) =~ y/a-z/A-Z/;
	close(IN);
    }

    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "> $tmpdir/$tag.tmp");
	    }
	}
    }
    my $fin = "in-sva";
    open($fin, "$tmpdir/$target.sv.sort");
    while(<$fin>){
	chomp;
	 @row = split('\t', $_);
	($hchr, $hpos, $tchr, $tpos, $direction, $type, $size) = (split(' ', $row[0]))[1.. 7];
	$hchr =~ s/^0+//;
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
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		close($tag);
	    }
	}
    }
    $semaphore->up;
    report("Making data for verification of sv complete");
}
    
sub svSort{
    report("Making $target.sv.sort");
    my ($flag, @row, $chr, $out, $count);
    my $fin = "in-sv";
    my $fout = "out-sv";
    open($fin, "$wd/$target/$target.aln");
    open($fout, "|sort $sort_opt -T $tmpdir -k 1 -k 2 -n |uniq > $tmpdir/$target.sv.sort");
    while(<$fin>){
	chomp;
	if (/ion/){
	    $flag = 1;
	    @row = split;
	    $chr = "00$row[1]";
	    $chr = substr($chr, length($chr) - 3, 3);
	    if ($row[7] eq ""){
            $out .= "$row[0] $chr $row[2] $row[3] $row[4] $row[5] $row[6]\t";
        }else{
            $out .= "$row[0] $chr $row[2] $row[3] $row[4] $row[5] $row[6] $row[7]\t";
        }
	    $count = 0;
	}elsif($flag and $count != 1){
	    if ($_ eq ""){
		print $fout "$out\n";
		$flag = 0;
		$out = "";
	    }else{
		$out .= "$_\t";
	    }
	}
	$count++;
    }
    close($fin);
    close($fout);
    report("Making $target.sv.sort complete");
    $semaphore->up;
}

sub index{
    report("Making index");
    my $pos  = 0;
    my ($length, @row);
    open(INDEXOUT, "|sort $sort_opt -T $tmpdir > $wd/$target/$target.index");
    open(INDEXIN, "$wd/$target/$target.aln");
    while(<INDEXIN>){
	$length = length($_);
	if(/^#/){
	    @row = split;
	    if ($row[1] =~ /^[0-9]*$/){
		$row[1] = "000$row[1]";
		$row[1] = substr($row[1], length($row[1]) -3, 3);
	    }
	    $row[2] = "0000000000$row[2]";
	    $row[2] = substr($row[2], length($row[2]) - 10, 10);
	    print INDEXOUT "$row[1] $row[2] S $pos\n";
	}
	$pos += $length;
    }
    close(INDEXIN);
    close(INDEXOUT);
    report("Making index complete");
    $semaphore->up;
}

sub sortData4Map{ 
    my ($nuca, $nucb, $nucc, $tag);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		threads->new(\&sortData4MapFunc, $tag);
	    }
	}
    }
    foreach my $thr (threads->list){
	$thr->join;
    }
}

sub sortData4MapFunc{
    my $tag = shift;
    report("Sorting sequence data. $tag");
    system("cat $tmpdir/$tag.tmp.* | sort $sort_opt -T $tmpdir  > $tmpdir/$tag");
    &waitFile("$tmpdir/$tag");
    system("rm $tmpdir/$tag.tmp.*");
    $semaphore->up;
}


sub mkData4MapF{
    my ($nuca, $nucb, $nucc, $tag);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		threads->new(\&mkData4MapFFunc, $tag);
		report("Making data for first mapping. $tag");
	    }
	}
    }
    foreach my $thr (threads->list){
	$thr->join;
    }
}

sub mkData4MapFFunc{
    my $tag = shift;
    my ($nuca, $nucb, $nucc, $subtag, $margin, $head_pos, $head);
    my $fin = $tag . "IN";
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		open($subtag, "> $tmpdir/$subtag.tmp.$tag")
	    }
	}
    }
    open($fin, "zcat $wd/$target/sort_uniq/$target.sort_uniq.$tag.gz 2> /dev/null |");
    while(<$fin>){
	chomp;
	foreach $margin ('0', '5', '10', '15'){
	    $head_pos = $margin + 1;
	    $head = substr($_, $head_pos -1, 20);
	    $subtag = substr($head, 0, 3);
	    print $subtag "$head $_ $head_pos\n";
	}
    }
    close($fin);  
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		close($subtag)
	    }
	}
    }
    $semaphore->up;
}

sub mkData4MapR{
    my ($nuca, $nucb, $nucc, $tag);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		report("Making data for second mapping. $tag");
		threads->new(\&mkData4MapRFunc, $tag);
	    }
	}
    }
    foreach my $thr (threads->list){
	$thr->join;
    }
}

sub mkData4MapRFunc{
    my $tag = shift;
    my ($nuca, $nucb, $nucc, $subtag, $margin, $tail_pos, $tail);
    my $fin = $tag . "IN";
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag =  $nuca . $nucb . $nucc;
		open($subtag, "> $tmpdir/$subtag.tmp.$tag")
	    }
	}
    }
    open($fin, "$tmpdir/$tag.map");
    while(<$fin>){
	chomp;
	$margin = (split)[1] - 1;
	$tail_pos = $length - 20 - $margin + 1;
	$tail = substr($_, $tail_pos - 1, 20);
	$subtag = substr($tail, 0, 3);
	print $subtag "$tail $_ $tail_pos\n";
    }
    close($fin);  
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		close($subtag)
	    }
	}
    }
     $semaphore->up;
}

sub map{
    my ($nuca, $nucb, $nucc, $subtag);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		close($tag);
		$semaphore->down;
		report("Mapping. $tag");
		threads->new(\&mapFunc, $tag);
	    }
	}
    }
    foreach my $thr (threads->list){
	$thr->join;
    }
}

sub mapFunc{
    my $tag = shift;
    system("zcat $refdir/ref20_uniq.$tag.gz | join $tmpdir/$tag - |cut -c 22- > $tmpdir/$tag.map");
    &waitFile("$tmpdir/$tag.map");
    system("rm $tmpdir/$tag");
    $semaphore->up;
}

sub align{
    my ($nuca, $nucb, $nucc, $subtag);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		close($tag);
		$semaphore->down;
		report("Bidirectional alignment. $tag");
		threads->new(\&alignFunc, $tag);
	    }
	}
    }
    foreach my $thr (threads->list){
	$thr->join;
    }
}

sub alignFunc{
    my $tag = shift;
    my ($seq, $hpos, $hchr, $head_pos, $head_direction, $tpos, $tchr, $tail_pos, $tail_direction, $length, $head, $tail, $hhit, $thit, $margin, $head_seq, $distance, $head_bar, $tail_bar, $head_space, $tail_space, @head, @tail, @seq, $mcount, $out, $i, $j, $k, $head_junction, $head_fail, $unmatch, $hcount, $head_junction, $tail_junction, $tail_fail, $tcount, $tail_direction, $type, $distance, $fin, $fout, $chr_file);
    $fin = $tag;
    $fout = "$tag-out";
    $chr_file = "$tag-chr";
    open($fout, "> $tmpdir/$target.aln.$tag");
    open($fin, "$tmpdir/$tag.map");
    while(<$fin>){
	chomp;
	($seq, $hpos, $hchr, $head_pos, $head_direction, $tpos, $tchr, $tail_pos, $tail_direction) = split;
	$length = length($seq);
	$head = substr($_, $hpos -1, 20);
	$tail = substr($_, $tpos -1, 20);
	$hhit = 0;
	$thit = 0;
	$margin = $hpos -1;

	open($chr_file, "$wd/$ref/chr$hchr");
	if ($head_direction eq "f"){
	    seek($chr_file, $head_pos - $margin -1, 0);
	    read($chr_file, $head_seq, $length);
	    $head_seq =~ y/a-z/A-Z/;
	    next if $head_seq eq $seq;
	    if ($hchr eq $tchr and $tail_direction eq "f"){
		$distance = $tail_pos - $head_pos;
		if ($distance == $tpos - $hpos){
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
			print $fout $out . "$head_seq
$head_bar
$seq

";
		    }
		    next;
		}
	    }else{
		$distance = "";
	    }
	    if ($hchr ne $tchr){
		open($chr_file, "$wd/$ref/chr$tchr");
	    }
	    if ($tail_direction eq "f"){
		seek($chr_file, $tail_pos - $length + 20 + $margin -1, 0);
		read($chr_file, $tail_seq, $length);
		$tail_seq =~ y/a-z/A-Z/;
	    }else{
		seek($chr_file, $tail_pos - 20 -$margin, 0);
		read($chr_file, $tail_seq, $length);
		$tail_seq =~ y/a-z/A-Z/;
		$tail_seq = &complement($tail_seq);
	    }
	    close($chr_file);
	    
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
		print  $fout "# $hchr $head_junction $tchr $tail_junction $tail_direction $type $distance
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
    close($fin);
    close($fout);
    $semaphore->up;
    system("rm $tmpdir/$tag.map");
}

sub report{
    my $message = shift;
    my $now = `date`;
    chomp($now);
    $message = "$now $message\n";
    print $message;
    print REPORT $message;
}

sub waitFile{
    my $file = shift;
    while(1){
	last if -e $file;
	sleep 10;
    }
    while(1){
	$mtime = (stat($file))[9];
	if (time > $mtime + 5){
	    return;
	}
	sleep 1;
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
