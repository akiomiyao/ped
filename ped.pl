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
# Author: MIYAO Akio <miyao@affrc.go.jp>
#

use threads;
use Thread::Semaphore;

$usage = '
 ped.pl - program for polymorphic edge detection. 

 e.g. perl ped.pl target=ERR194147,ref=hg38

 If you want to detect polymorphisms between target and control,
 perl ped.pl target=ERR3063487,control=ERR3063486,ref=WBcel235

 If you want to specify the working directory,
 perl ped.pl target=ERR194147,ref=hg38,wd=/home/you/work

 If you want to specify both the working directory and tmp directory,
 perl ped.pl target=ERR194147,ref=hg38,wd=/home/you/work,tmpdir=/mnt/ssd

 If you want to set muximum number of threads,
 perl ped.pl target=ERR194147,ref=hg38,thread=14
 In the case of file open error, reduce muximum number of threads.
 Because default setting of ulimit is 1024, 14 threads is muximun for
 defalt setting of OS. 

 For kmer method,
 perl ped.pl target=ERR3063487,control=ERR3063486,ref=WBcel235,method=kmer

 Results will be saved in the target directory.

';

my $sort_opt = "-S 1M";
@nuc = ('A', 'C', 'G', 'T');
$ENV{LANG} = "C";

if ($ARGV[0] =~ /target|ref/){
    my @arg = split(',', $ARGV[0]);
    foreach (sort @arg){
	next if $_ eq "";
	my ($name, $val) = split('=', $_);
	$$name = $val;
	$log .= "# $name: $val\n";
    }
}elsif ($ARGV[1] ne ""){
    $target = $ARGV[0];
    $control = $ARGV[1];
    $ref = $ARGV[2];
    $method = $ARGV[3];
}elsif ($ENV{target} ne ""){
    $target = $ENV{target};
    $control = $ENV{control};
    $ref = $ENV{ref};
    $method = $ENV{method};
    $wd = $ENV{wd};
    $tmpdir = $ENV{tmpdir};
    $max_semaphore = $ENV{thread};
}else{
    print $usage;
}

if ($wd eq ""){
    $wd = `pwd`;
    chomp($wd);
}

$uname = `uname`;
chomp($uname);

if ($method eq ""){
    if ($ref eq ""){
	$method = "kmer";
    }else{
	$method = "bidirectional";
    }
}

$zcat = "zcat";

if ($uname eq "FreeBSD"){
    die "curl not found. Please install curl." if ! -e "/usr/local/bin/curl";
    open(IN, "sysctl kern.smp.cpus |");
    while(<IN>){
	chomp;
	$processor = (split(': ', $_))[1];
    }
    close(IN);
}elsif($uname eq "Darwin"){
    $zcat = "gzcat";
    open(IN, "sysctl hw.logicalcpu |");
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
}else{
    die "Operating system is unknown.\n";
}

if ($uname eq "Darwin"){
    $max_semaphore = 1;
}elsif ($processor >= 8 ){
    $max_semaphore = 8;
}else{
    $max_semaphore = $processor;
}
$max_semaphore = 3 if $max_semaphore eq "";
$max_semaphore = $thread if $thread ne "";
$max_semaphore = 14 if $max_semaphore > 14;
$semaphore = Thread::Semaphore->new($max_semaphore);

$refdir = "$wd/$ref";
if ($tmpdir eq ""){
    system("mkdir $wd/$target/tmp") if ! -e "$wd/$target/tmp";
    $tmpdir = "$wd/$target/tmp";
}
$control = $ref if $control eq "" or $control eq "default";

if (! -e "config"){
    system("curl -O https://raw.githubusercontent.com/akiomiyao/ped/master/config");
}

open(IN, "config");
while(<IN>){
    chomp;
    my @row = split('\t', $_);
    if ($row[1] eq "description"){
	$desc{$row[0]} = $row[2];
    }elsif($row[1] eq "wget"){
	$curl{$row[0]} = $row[2];
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

if ($ARGV[0] eq ""){
    print " Currently, reference genomes listed below are suported.

  Name           Description\n";
    foreach (sort keys %desc){
	my $name = "$_                  ";
	$name = substr($name, 0, 14);
	print "  $name $desc{$_}\n";
    }
    print "
 Author: MIYAO Akio <miyao\@affrc.go.jp>

";
    exit;
}

open(REPORT, "> $wd/$target/$target.log");
print REPORT "# Log of ped.pl
$log";
report("Job begin: $method method");


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

if ($ref ne "" and ! -e "$wd/$ref/sort_uniq/$ref.sort_uniq.TTT.gz" and ! -e "$wd/$ref/sort_uniq/$ref.TTT.gz"){
    &mkRef;
}

if (! -e "$wd/$target/sort_uniq/$target.sort_uniq.TTT.gz"){
    &mkSortUniq($target);
}

if ($control ne $ref && ! -e "$wd/$control/sort_uniq/$control.sort_uniq.TTT.gz"){
    &mkSortUniq($control);
}

open(IN, "$zcat $wd/$target/sort_uniq/*.gz 2> /dev/null |");
while(<IN>){
    chomp;
    $length = length($_);
    last;
}
close(IN);

open(IN, "$zcat $wd/$control/sort_uniq/*.gz 2> /dev/null |");
while(<IN>){
    chomp;
    $clength = length($_);
    last;
}
close(IN);

if ($method eq "kmer"){
    &kmer;
}else{
    &bidirectional;
}

if ($tmpdir eq "$wd/$target/tmp"){
    system("rm -r $tmpdir");
}else{
    system("rm $tmpdir/*");
}
system("rm -r $wd/$control/tmp") if -e "$wd/$control/tmp";
report("Job completed: $method method.");

sub kmer{
    &countKmer($target);
    &countKmer($control);
    &joinKmer;
    if($ref eq ""){
	system("cat $tmpdir/$target.snp.* > $wd/$target/$target.kmer.snp");
	system("rm $tmpdir/$target.snp.*");
	return;
    }
    &mapKmer;
    &snpMkT;
    &sortSeq;
    &joinTarget;
    &snpMkC;
    &sortSeq;
    &joinControl;
    &kmerReadCount;
    &toVcf;
    &primer;
}

sub bidirectional{
    &mkData4MapF;
    &sortData4Map;
    &map;
    &mkData4MapR;
    &sortData4Map;
    &map;
    &align;

    $semaphore->down;
    threads->new(\&index);
    $semaphore->down;
    threads->new(\&svSort);
    &joinAll;

    &svMkT;
    &sortSeq;
    &joinTarget;
    &svMkC;
    &sortSeq;
    &joinControl;
    &svReadCount;
    &snpMkT;
    &sortSeq;
    &joinTarget;
    &snpMkC;
    &sortSeq;
    &joinControl;
    &snpReadCount;
    &toVcf;
    &primer(sv);
    &primer;
}

sub primer{
    my $type = shift;
    my $count = "0001";
    my $line = 0;
    my $num;
    if ($type eq "sv"){
	open(IN, "$wd/$target/$target.sv");
	open(OUT, "> $wd/$target/tmp/$target.$count");
    }else{
	if ($method eq "kmer"){
	    open(IN, "$wd/$target/$target.kmer.snp");
	    open(OUT, "> $wd/$target/tmp/$target.$count");
	    $type = "kmer";
	}else{
	    open(IN, "$wd/$target/$target.bi.snp");
	    open(OUT, "> $wd/$target/tmp/$target.$count");
	    $type = "bi";
	}
    }
    while(<IN>){
	next if ! /H|M/;
	$line++;
	print OUT;
	if ($line == 10000){
	    $line = 0;
	    close(OUT);
	    $count++;
	    $num = "000$count";
	    $num = substr($num, length($num) - 4, 4);
	    open(OUT, "> $wd/$target/tmp/$target.$num");
	}
    }
    close(OUT);
    for ($i = 1; $i <= $count; $i++){
	$semaphore->down;
	$num = "000$i";
	$num = substr($num, length($num) - 4, 4);
	&report("Output primer sequences. $num");
	threads->new(\&primerFunc, $num, $type);
    }
    &joinAll;
    system("cat $wd/$target/tmp/primer.* > $wd/$target/$target.$type.primer && rm $wd/$target/tmp/*");
}

sub primerFunc{
    my(@row, $seq, $i, $f, $r, @f, @r, $gc, $flag, $hpos, $tpos, $fpos, $rpos, $length, $tg, $head, $tail, $minimum_length, $out, $fin, $fout, $fchr, $fchrt);
    my ($num, $type) = @_;
    $fin = "in.$num";
    $fout = "out.$num";
    $fchr = "chr.$num";
    $fchrt = "chrt.$num";
    open($fin, "$wd/$target/tmp/$target.$num");
    open($fout, "> $wd/$target/tmp/primer.$num");
    while(<$fin>){
	@row = split('\t', $_);
	$dat = $_;
	chomp($dat);
	if ($prev_chr ne $row[0]){
	    open($fchr, "$wd/$ref/chr$row[0]");
	    binmode($fchr);
	}
	$prev_chr = $row[0];
	@f = ();
	@r = ();
	
	if ($type eq "sv"){
	    if($row[5] =~ /ins|del/){
		if ($row[1] <= $row[3]){
		    $hpos = $row[1];
		    $tpos = $row[3]
		}else{
		    $hpos = $row[3];
		    $tpos = $row[1]
		}
	    }else{
		    $hpos = $row[1];
		    $tpos = $row[3]
	    }
	    seek($fchr, $hpos - 250, 0);
	    read($fchr, $head, 250);
	    if ($row[0] ne $row[2]){
		open($fchrt, "$wd/$ref/chr$row[2]");
		binmode($fchrt);
		seek($fchrt, $row[3], 0);
		read($fchrt, $tail, 250);
		close($fchrt);
	    }elsif($row[5] =~ /inversion/){
		seek($fchr, $tpos, 0);
		read($fchr, $tail, 250);
		$tail = &complement($tail);
	    }else{
		seek($fchr, $tpos, 0);
		read($fchr, $tail, 250);
	    }
	    for($i = 0; $i < 230; $i++){
		$f = substr($head, $i, 20);
		($gc = $f) =~ y/AT//d;
		if (length($gc) == 11){
		    push(@f, $f);
		}
	    }
	    for($i = 0; $i < 230; $i++){
		$r = substr($tail, $i, 20);
		($gc = $r) =~ y/AT//d;
		if (length($gc) == 11){
		    push(@r, &complement($r));
		}
	    }
	}else{
	    seek($fchr, $row[1] - 250, 0);
	    read($fchr, $seq, 500);
	    
	    for($i = 0; $i < 230; $i++){
		$f = substr($seq, $i, 20);
		($gc = $f) =~ y/AT//d;
		if (length($gc) == 11){
		    push(@f, $f);
		}
	    }
	    for($i = 250; $i < 480; $i++){
		$r = substr($seq, $i, 20);
		($gc = $r) =~ y/AT//d;
		if (length($gc) == 11){
		    push(@r, &complement($r));
		}
	    }
	}
	$minimum_length = 10000;
	$out = "";
	$minimum_count = 0;
	foreach $f (reverse @f){
	    foreach $r (@r){
		if (&checkDimer($f, $r)){
		    if ($type eq "sv"){
			$fpos = index($head, $f, 0);
			$rpos = index($tail, complement($r), 0);
			$length = 250 - $fpos + $rpos + length($row[12]);
			$tg = substr($head, 229, 20) . " " . substr($tail, 0, 20);
			next if $length < 70;
			if ($minimum_length > $length){
			    $minimum_length = $length;
			    $out = "$dat\t$f\t$r\t$length\t$tg\n";
			    $minimum_count ++;
			    last if $minimum_count >= 3;
			}
		    }else{
			$fpos = index($seq, $f, 0);
			$rpos = index($seq, complement($r), 0);
			$length = $rpos - $fpos + 20;
			next if $length < 70;
			@row = split('\t', $dat);
			$tg = substr($seq, 99, 150) . "[$row[3]/$row[2]]" . substr($seq, 250, 149);
		    }
		    if ($minimum_length > $length){
			$minimum_length = $length;
			$out = "$dat\t$f\t$r\t$length\t$tg\n";
			$minimum_count ++;
			last if $minimum_count >= 3;
		    }
		}	   
	    }
	    last if $minimum_count >= 3;
	}
	if ($out ne ""){
	    print $fout $out;
	}
    }
    close($fin);
    close($fchr);
    close($fout);
    $semaphore->up;
}

sub checkDimer{
    my ($f, $r) = @_;
    my @f = split('', $f);
    my @r = split('', $r);
    for ($i = 3; $i < 20; $i++){
        my $match = 0;
        for ($j = 0; $j < $i; $j++){
            my $fn = $f[20 - $i + $j];
            my $rn = &complement($r[19 - $j]);
            if ($fn eq $rn){
                $match ++;
            }
        }
        if ($match >= $i - 1){
            return 0;
        }
    }
    return 1;
}

sub mapKmer{
    my ($tag, $nuca, $nucb, $nucc);
   
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		&report("Mapping SNPs between $target and $control. $tag");
		threads->new(\&mapKmerSub, $tag);
	    }
	}
    }
    &joinAll;
}

sub mapKmerSub{
    my $tag = shift;
    my $fin = "map.in.$tag";
    my $fout = "map.out.$tag";
    my $fref = "map.ref.$tag";
    my (@row, %dat, $seq, $nuc, $pos, @dat);
    open($fin, "$tmpdir/$target.snp.$tag");
    while(<$fin>){
	chomp;
	@row = split;
	$seq = shift(@row);
	if (length($row[0]) == 1){
	    $dat{$seq} = "$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\t$row[9]";
	}
    }
    close($fin);
    open($fout, "> $tmpdir/$target.map.$tag");
    open($fref, "$zcat $wd/$ref/ref20_uniq.$tag.gz|");
    while(<$fref>){
	chomp;
	@row = split;
	$seq = substr($row[0], 0, 19);
	$nuc = substr($row[0], 19, 1);
	if ($dat{$seq}){
	    @dat = split('\t', $dat{$seq});
	    if($row[3] eq "f"){
		$pos = $row[2] + 19;
		if ($nuc eq $dat[0]){
		    print $fout "$row[1]\t$pos\t$seq\t$nuc\t$dat[0]\t$dat[1]\t$row[3]\t$dat[2]\t$dat[3]\t$dat[4]\t$dat[5]\t$dat[6]\t$dat[7]\t$dat[8]\t$dat[9]\n";
		}
	    }else{
		$pos = $row[2] - 19;
		$dat[0] = complement($dat[0]);
		$dat[1] = complement($dat[1]);
		$seq = complement($seq);
		$nuc = complement($nuc);
		if ($nuc eq $dat[0]){
		    print $fout "$row[1]\t$pos\t$seq\t$nuc\t$dat[0]\t$dat[1]\t$row[3]\t$dat[2]\t$dat[3]\t$dat[4]\t$dat[5]\t$dat[6]\t$dat[7]\t$dat[8]\t$dat[9]\n";
		}
	    }
	}
    }
    close($fout);
    close($fref);
    $semaphore->up;
}

sub joinKmer{
    my ($tag, $nuca, $nucb, $nucc);
   
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		&report("Detection SNPs between $target and $control. $tag");
		threads->new(\&joinKmerSub, $tag);
	    }
	}
    }
    &joinAll;
}

sub joinKmerSub{
    my $tag = shift;
    $cutoff = 10 if $cutoff eq "";
    my $fin = "join.in.$tag";
    my $fout = "join.out.$tag";
    system("$zcat $wd/$target/lbc/$target.lbc.$tag.gz > $tmpdir/$target.lbc.$tag");
    system("$zcat $wd/$control/lbc/$control.lbc.$tag.gz > $tmpdir/$control.lbc.$tag");

    open($fout, "> $tmpdir/$target.snp.$tag");
    open($fin, "join $tmpdir/$control.lbc.$tag $tmpdir/$target.lbc.$tag |");
    while(<$fin>){
	my $rep = 0;
	my $a = "";
	my $b = "";
	my $nohita = 0;
	my $nohitb = 0;
	my $pola = 0;
	my $polb = 0;
	my @row = split;
	my $i;
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
	    print $fout "$row[0]\t$a\t$b\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\n";
	}
    }
    close($fin);
    close($fout);
    system("rm $tmpdir/$target.lbc.$tag $tmpdir/$control.lbc.$tag");
    $semaphore->up;
}

sub countKmer{
    my $target = shift;
    my ($tag, $nuca, $nucb, $nucc);
    return if -e "$wd/$target/lbc/$target.lbc.TTT.gz";
    if (! -e "$tmpdir"){
	system("mkdir $tmpdir");
    }
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		&report("Making kmer for $target. $tag");
		&countKmerSub($target, $tag);
	    }
	}
    }

    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		&report("Making kmer for $target. Sorting of $tag subfiles");
		threads->new(\&countKmerSort, $target, $tag);
	    }
	}
    }
    &joinAll;

     foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		&report("Joining kmer for $target. $tag");
		threads->new(\&countKmerMerge, $target, $tag);
	    }
	}
    }
    &joinAll;

    system("mkdir $wd/$target/lbc") if ! -e "$wd/$target/lbc";
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		&report("Making last base counts for $target. $tag");
		threads->new(\&countKmerLbc, $target, $tag);
	    }
	}
    }
    &joinAll;
}

sub countKmerLbc{
    my ($target, $tag) = @_;
    my ($fin, $fout, $seq, $count, $head, $prev, $nuc, %count);

    $fin = "$tag.lbc.in";
    $fout = "$tag.lbc.out";
    open($fout, "|gzip > $wd/$target/lbc/$target.lbc.$tag.gz");
    open($fin, "$tmpdir/$tag.count");
    while(<$fin>){
	chomp;
	($seq, $count) = split;
        $head = substr($seq, 0, 19);
        $nuc = substr($seq, 19, 1);
        
        if ($head ne $prev and $prev ne ""){
            print $fout "$prev\t$count{A}\t$count{C}\t$count{G}\t$count{T}\n";
            $count{A} = 0;
            $count{C} = 0;
            $count{G} = 0;
            $count{T} = 0;
        }
        $count{$nuc} = $count;
        
        $prev = $head;
    }
    print $fout "$prev\t$count{A}\t$count{C}\t$count{G}\t$count{T}\n";
    close($fin);
    close($fout);
    system("rm $tmpdir/$tag.count");
    $semaphore->up;
}

sub countKmerMerge{
    my ($target, $parent) = @_;
    my ($nuca, $nucb, $nucc, $tag, $fin, $fout, @row, $count);
    system("touch $tmpdir/$parent.count");
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		&report("Making kmer for $target. $parent Merging of $tag subfile");
		$fin = "$parent.$tag.join.in";
		$fout = "$parent.$tag.join.out";
		open($fout, "> $tmpdir/$parent.count.tmp");
		open($fin, "$zcat $tmpdir/$parent.$tag.count.gz |join -a 1 -a 2 $tmpdir/$parent.count - |");
		while(<$fin>){
		    chomp;
		    @row = split;
		    $count = $row[1] + $row[2];
		    print $fout "$row[0]\t$count\n";
		}
		close($fin);
		close($fout);
		system("mv $tmpdir/$parent.count.tmp $tmpdir/$parent.count && rm $tmpdir/$parent.$tag.count.gz");
	    }
	}
    }
    $semaphore->up;
}

sub countKmerSort{
    my ($target, $parent) = @_;
    my ($nuca, $nucb, $nucc, $tag);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		&report("Making kmer for $target. Sorting of $parent.$tag subfile");
		system("$zcat $tmpdir/$parent.$tag.gz | sort $sort_opt -T $tmpdir |uniq -c | awk '{print \$2 \"\t\" \$1}'| gzip > $tmpdir/$parent.$tag.count.gz && rm $tmpdir/$parent.$tag.gz");
	    }
	}
    }
    $semaphore->up;
}

sub countKmerSub{
    my ($target, $parent) = @_;
    my ($nuca, $nucb, $nucc, $fin, $fout, $length, $seq, $ktag, $i);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$fout = "$tag.$parent.out";
		open($fout, "|gzip > $tmpdir/$tag.$parent.gz");
	    }
	}
    }
    $fin = "$parent.in";
    if (-e "$wd/$target/sort_uniq/$target.sort_uniq.$parent.gz"){
	open($fin, "$zcat $wd/$target/sort_uniq/$target.sort_uniq.$parent.gz 2> /dev/null |");
    }else{
	open($fin, "$zcat $wd/$target/sort_uniq/$target.$parent.gz 2> /dev/null |");
    }
    while(<$fin>){
	chomp;
	$length = length($_);
	for ($i = 0; $i <= $length - 20; $i++){
	    $seq = substr($_, $i, 20);
	    if ($seq !~ /N/){
		$ktag = substr($seq, 0, 3);
		$fout = "$ktag.$parent.out";
		print $fout "$seq\n";
	    }
	}
    }
    close($fin);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$fout = "$tag.$parent.out";
		close($fout);
	    }
	}
    }
}

sub mkRef{
    my (@row, $i, $remote_file);
    if (! -e "$wd/$ref"){
	system("mkdir $wd/$ref");
    }
    if ($file eq ""){
	@row = split('/', $curl{$ref});
	$remote_file = $row[$#row];
	if (! -e "$wd/$ref/$remote_file" and $curl{$ref} ne ""){
	    &report("Downloading $curl{$ref}");
	    system("cd $wd/$ref && curl -O $curl{$ref} && cd $wd");
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
    chdir "$wd/$ref";
    die "$file is not found in $wd/$ref." if ! -e "$wd/$ref/$file";
    @chr = ();
    if ($file =~ /gz$/){
	open(IN, "$zcat $wd/$ref/$file|");
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
	    $out += 0 if $out =~ /^[0-9]+$/;
	    push(@chr, $out) if ! $chr_flag;
	    $out = "chr$out";
	    open(OUT, "> $out");
	}else{
	    y/a-z/A-Z/;
	    y/ACGT/N/c;
	    print OUT;
	}
    }
    close(IN);
    close(OUT);
    chdir $wd;
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
    binmode($fin);
    $i = 0;
    while(1){
	seek($fin, $i, 0);
	read($fin, $read, 100);
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
		$fout = "$tag.$chr";
		close($fout);
	    }
	}
    }
    $semaphore->up;
}

sub mkControlReadSub{
    my $tag = shift;
    my ($nuca, $nucb, $nucc, $subtag, $fin, $ftmp, $fout);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		$ftmp = $tag . $subtag;
		open($ftmp, "> $wd/$ref/tmp/tmp.$ftmp");
	    }
	}
    }
    $fin = $tag;
    open($fin, "cat $wd/$ref/tmp/$tag.* |");
    while(<$fin>){
	$ftmp = substr($_, 0, 6);
	print $ftmp $_;
    }
    close($fin);
    $fout = "tmpout.$tag";
    open($fout, "|uniq |gzip -f > $wd/$ref/sort_uniq/$ref.sort_uniq.$tag.gz");
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		$ftmp = $tag . $subtag;
		close($ftmp);
		$fin = "tmp.$tag";
		open($fin, "sort $sort_opt -T $wd/$ref/tmp $wd/$ref/tmp/tmp.$ftmp |");
		while(<$fin>){
		    print $fout $_;
		}
		close($fin);
		system("rm $wd/$ref/tmp/tmp.$ftmp");
	    }
	}
    }
    close($fout);
    system("rm $wd/$ref/tmp/$tag.*");
    $semaphore->up;
}

sub mkUniq{
    my $tag = shift;
    my ($nuca, $nucb, $nucc, $subtag, $fin, $fout, $funiq);
    &report("Making ref20_uniq.$tag.gz");
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		$fout = $tag . $subtag;
		open($fout, "> $wd/$ref/tmp/tmp.$fout");
	    }
	}
    }
    $fin = $tag;
    open($fin, "cat $wd/$ref/tmp/ref20.$tag.* |");
    while(<$fin>){
	$fout = substr($_, 0, 6);
	print $fout $_;
    }
    close($fin);
    $funiq = "tmpout.$tag";
    open($funiq, "|gzip -f > $wd/$ref/ref20_uniq.$tag.gz");
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		$fout = $tag . $subtag;
		close($fout);
		$count = 0;
		$fin = "tmp.$tag";
		open($fin, "sort $sort_opt -T $tmpdir $wd/$ref/tmp/tmp.$fout |");
		while(<$fin>){
		    chomp;
		    @row = split;
		    if ($prev ne "" and $prev ne $row[0]){
			print $funiq "$pline\n" if $count == 1;
			$count =0;
		    }
		    $prev = $row[0];
		    $pline = $_;
		    $count++;
		}
		print $funiq "$pline\n" if $count == 1;
		close($fin);
		system("rm $wd/$ref/tmp/tmp.$fout");
	    }
	}
    }
    close($funiq);
    system("rm $wd/$ref/tmp/ref20.$tag.*");
    $semaphore->up;
}

sub mk20mer{
    my ($fin, $fout, $i, @tag, $tag, $forward, $length, $comp, $fpos, $rpos);
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

    my $file = "$wd/$ref/chr$chr";
    open(IN, $file);
    binmode(IN);
    $length = -s $file;
    for ($i = 0; $i <= $length - 20; $i++){
	seek(IN, $i, 0);
	read(IN, $forward, 20);
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
    close(IN);
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		$fout = "$tag.$chr";
		close($fout);
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
}

sub mkChr{
    my $output;
    my @file = split('/', $curl{$ref});
    my $file = $file[$#file];
    $file =~ s/\ +$//g;
    my $i = 0;
    chdir "$wd/$ref";
    &report("Making chromosome file from $file");
    if ($ref eq "hg38"){
	open(IN, "$zcat $file|");
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
		    $output += 0 if $output =~ /^[0-9]+$/;
		    $output = "chr$output";
		    open(OUT, "> $output");
		}
	    }elsif(! $flag){
		y/a-z/A-Z/;
		y/ACGT/N/c;
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
		$file = $curl{$ref} . "chr$i.fna.gz";
		&report("Downloading $file");
		system("curl -O $file");
		$file = "chr$i.fna.gz";
		open(IN, "$zcat $file|");
		open(OUT, "> chr$i");
		while(<IN>){
		    chomp;
		    if (! /^>/){
			y/a-z/A-Z/;
			y/ACGT/N/c;
			print OUT;
		    }
		}
		close(IN);
		close(OUT);
	    }
	    chdir $wd;
	    return;
	}
	
	if ($file =~ /gz$/){
	    open(IN, "$zcat $file|");
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
		y/ACGT/N/c;
		print OUT;
	    }
	}
	close(IN);
	close(OUT);
    }
    chdir $wd;
}

sub mkSortUniq{
    my $subject = shift;
    my ($cmd, $gz_file, $bz_file, $xz_file, $fq_file);
    report("Making sort_uniq files for $subject");
    opendir(DIR, "$wd/$subject/read");
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

    if ($gz_file ne ""){
	$cmd = "cd $wd/$subject/read && $zcat $gz_file |";
	&sortUniqSub($cmd, $subject);
    }
    if($bz_file ne ""){
	$cmd = "cd $wd/$subject/read && bzcat $bz_file |";
	&sortUniqSub($cmd, $subject);
    }
    if($xz_file ne ""){
	$cmd = "cd $wd/$subject/read && xzcat $xz_file |";
	&sortUniqSub($cmd, $subject);
    }
    if($fq_file ne ""){
	$cmd = "cd $wd/$subject/read && cat $fq_file |";
	&sortUniqSub($cmd, $subject);
    }

    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		$semaphore->down;
		threads->new(\&sortUniqSort, $tag, $subject);
	    }
	}
    }
    &joinAll;
}

sub sortUniqSort{
    my ($tag, $subject) = @_;
    my ($nuca, $nucb, $nucc, $subtag);
    report("Making $subject.sort_uniq files: Sorting for $subject.sort_uniq.$tag.gz");
    system("rm $subject.sort_uniq.$tag") if -e "$subject.sort_uniq.$tag";
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $tag . $nuca . $nucb . $nucc;
		system("sort -T $wd/$subject/sort_uniq/ $sort_opt $wd/$subject/sort_uniq/$subtag.tmp | uniq >> $wd/$subject/sort_uniq/$subject.sort_uniq.$tag && rm $wd/$subject/sort_uniq/$subtag.tmp");
	    }
	}
    }
    system("gzip $wd/$subject/sort_uniq/$subject.sort_uniq.$tag");
    $semaphore->up;
}

sub sortUniqSub{
    my ($cmd, $subject) = @_;
    my $count = 0;
    my $total = 0;
    my ($nuca, $nucb, $nucc, $nucd, $nuce, $nucf, $tag, $subtag, $count, $complement, $total);
    system("mkdir $wd/$subject/sort_uniq") if ! -e "$wd/$subject/sort_uniq";
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, " > $wd/$subject/sort_uniq/$tag.tmp");
	    }
	}
    }
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
		report("Making $subject.sort_uniq files: Split to subfiles. $total reads processed");
	    }
	}elsif($count == 4){
	    $count = 0;
	}
	$count++;
    }
    close(IN);
    &closeTag;
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "$wd/$subject/sort_uniq/$tag.tmp");
		foreach $nucd (@nuc){
		    foreach $nuce (@nuc){
			foreach $nucf (@nuc){
			    $subtag = $tag . $nucd . $nuce . $nucf;
			    open ($subtag, "> $wd/$subject/sort_uniq/$subtag.tmp");
			}
		    }
		}
		while(<$tag>){
		    $subtag = substr($_, 0, 6);
		    print $subtag $_;
		}
		close($tag);
		system("rm $wd/$subject/sort_uniq/$tag.tmp");
		foreach $nucd (@nuc){
		    foreach $nuce (@nuc){
			foreach $nucf (@nuc){
			    $subtag = $tag. $nucd . $nuce . $nucf;
			    close($subtag);
			}
		    }
		}
	    }
	}
    }
}

sub toVcf{
    report("Convert to vcf format");
    my $fin = "in-vcf";
    my $fout ="out-vcf";
    if ($method eq "kmer"){
	open($fin, "$wd/$target/$target.kmer.snp");
	open($fout, "> $wd/$target/$target.kmer.vcf");
    }else{
	open($fin, "$wd/$target/$target.bi.snp");
	open($fout, "> $wd/$target/$target.bi.vcf");
    }
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
	if($method eq "kmer"){
	    $dp = $row[17] + $row[18];
	    next if $dp == 0;
	    $af = int(1000 * $row[18]/$dp)/1000;
	    if (/M/){
		print $fout "$row[0]\t$row[1]\t.\t$row[4]\t$row[5]\t1000\tPASS\tGT=1/1;AF=$af;DP=$row[18]\tGT:AD:DP\t1/1:$row[17],$row[18]:$dp\n";
	    }elsif(/H/){
		print $fout "$row[0]\t$row[1]\t.\t$row[4]\t$row[5]\t1000\tPASS\tGT=0/1;AF=$af;DP=$row[18]\tGT:AD:DP\t0/1:$row[17],$row[18]:$dp\n";
	    }
	}else{
	    $dp = $row[6] + $row[7];
	    next if $dp == 0;
	    $af = int(1000 * $row[7]/$dp)/1000;
	    if (/M/){
		print $fout "$row[0]\t$row[1]\t.\t$row[2]\t$row[3]\t1000\tPASS\tGT=1/1;AF=$af;DP=$dp\tGT:AD:DP\t1/1:$row[6],$row[7]:$dp\n";
	    }elsif(/H/){
		print $fout "$row[0]\t$row[1]\t.\t$row[2]\t$row[3]\t1000\tPASS\tGT=0/1;AF=$af;DP=$dp\tGT:AD:DP\t0/1:$row[6],$row[7]:$dp\n";
	    }
	}
    }
}

sub snpMkT{
    report("Making data for verification of snp. target");
    my (@row, $chr, $pos, $ref, $alt, $count, $ref_seq, $mut_seq, $head, $tail, $i, $ipos, $iref, $ialt, $tpos, $tw, $tm, $tag, $nuca, $nucb, $nucc, @dat);

    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "|gzip > $tmpdir/$tag.tmp.gz");
	    }
	}
    }
    my $fin = "in-mkt";
    my $fout = "out-mkt";
    open($fout, "|sort $sort_opt -T $tmpdir | uniq -c > $tmpdir/snp.list");
    if ($method eq "kmer"){
	open($fin, "cat $tmpdir/$target.map.* |");
    }else{
	open($fin, "$wd/$target/$target.aln");
    }
    while(<$fin>){
	chomp;
	if ($method eq "kmer"){
	    @row = split;
	    if ($row[0] =~ /^[0-9]+$/){
		$row[0] = "000" . $row[0];
		$row[0] = substr($row[0], length($row[0]) -3, 3);
	    }
	    $row[1] = "000000000000" . $row[1];
	    $row[1] = substr($row[1], length($row[1]) -11, 11);
	    $ref = $row[4];
	    $alt = $row[5];
	    $alt =~ s/$ref//;
	    print $fout "$row[0] $row[1] $ref $alt\n";
	}elsif (/^#/ and /snp/){
	    @row = split;
	    if ($row[1] =~ /^[0-9]+$/){
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

    open($fin, "$tmpdir/snp.list");
    while(<$fin>){
	chomp;
	if ($method eq "kmer"){
	    ($count, $chr, $pos, $ref, $alt) = split;
	    $count = 25;
	}else{
	    ($count, $chr, $pos, $ref, $alt) = split;
	}
	$chr =~ s/^0+//g;
	if ($prev_chr ne $chr){
	    my $chr_file = "$refdir/chr$chr";
	    open (IN, $chr_file);
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
	seek(IN, $pos - $length, 0);
	read(IN, $ref_seq, $length * 2 -1);
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
    close(IN);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		close($tag);
	    }
	}
    }
}

sub snpMkC{
    my (@row, $chr, $pos, $ref, $alt, $count, $chr_seq, $ref_seq, $mut_seq, $head, $tail, $i, $ipos, $iref, $ialt, $tpos, $cw, $cm, $tag, $nuca, $nucb, $nucc, @dat);

    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "|gzip > $tmpdir/$tag.tmp.gz");
	    }
	}
    }
    my $fin = "in-mkc";

    open($fin, "$tmpdir/snp.list");
    while(<$fin>){
	chomp;
	if ($method eq "kmer"){
	    ($count, $chr, $pos, $ref, $alt) = split;
	    $count = 25;
	}else{
	    ($count, $chr, $pos, $ref, $alt) = split;
	}
	$chr =~ s/^0+//g;
	if ($prev_chr ne $chr){
	    my $chr_file = "$refdir/chr$chr";
	    open (IN, $chr_file);
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
	seek(IN, $pos - $clength, 0);
	read(IN, $ref_seq, $clength * 2 -1);
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
    close(IN);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		close($tag);
	    }
	}
    }
}

sub svReadCount{
    my (@row, $dat, %count, $cm, $cw, $tm, $tw, $count, @prev, $prev, $genotype);
    open(OUT, "> $wd/$target/$target.sv");
    open(IN, "cat $tmpdir/*.target $tmpdir/*.control |");
    while(<IN>){
	chomp;
	@row = split;
	if ($row[0] =~ /^[0-9]+$/){
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
	$row[0] =~ s/^0+//g;
	$row[1] += 0;
	$dat = join("\t", @row[0 .. $#row -1]);
	if($dat ne $prev){
	    if ($cw >= 5 and $cm <= 1){
		if ($tm >= 5 and $tw <= 1){
		    $genotype = 'M';
		}elsif(($tm + $tw) > 0 and $tm / ($tm + $tw) > 0.3 and $tm / ($tm + $tw) < 0.7 and $cm <= 1 and $tw >= 5){
		    $genotype = "H";
		}else{
		    $genotype = "R";
		}
	    }else{
		$genotype = "N";
	    }
	    @prev = split('\t', $prev);
	    if ($prev[5] =~ /translocation|inversion/){
		$prev[7] = $prev[6];
		$prev[6] = "N";
	    }
	    $prev[7] = "_" if $prev[7] eq "";
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
    if ($cw >= 5 and $cm <= 1){
	if ($tm >= 5 and $tw <= 1){
	    $genotype = 'M';
	}elsif(($tm + $tw) > 0 and $tm / ($tm + $tw) > 0.3 and $tm / ($tm + $tw) < 0.7 and $cm <= 1 and $tw >= 5){
	    $genotype = "H";
	}else{
	    $genotype = "R";
	}
	$genotype = "N";
    }
    @prev = split('\t', $prev);
    if ($prev[5] =~ /translocation|inversion/){
	$prev[7] = $prev[6];
	$prev[6] = "N";
    }
    $prev[7] = "_" if $prev[7] eq "";
    print OUT "$prev[0]\t$prev[1]\t$prev[2]\t$prev[3]\t$prev[4]\t$prev[5]\t$prev[6]\t$cw\t$cm\t$tw\t$tm\t$genotype\t$prev[7]\n";
    close(OUT);
}

sub kmerReadCount{
    my (@row, $dat, %count, $cm, $cw, $tm, $tw, $count, @prev, $prev, $genotype);
    open(OUT, "| sort $sort_opt -T $tmpdir > $tmpdir/$target.map");
    open(IN, "cat $tmpdir/$target.map.* |");
    while(<IN>){
	chomp;
	@row = split;
	if ($row[0] =~ /^[0-9]+$/){
	    $row[0] = "000" . $row[0];
	    $row[0] = substr($row[0], length($row[0]) -3, 3);
	}
	$row[1] = "000000000000" . $row[1];
	$row[1] = substr($row[1], length($row[1]) -11, 11);
	print OUT "$row[0]:$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\t$row[9]\t$row[10]\t$row[11]\t$row[12]\t$row[13]\t$row[14]\n";
    }
    close(IN);
    close(OUT);

    open(OUT, "> $tmpdir/$target.verify");
    open(IN, "cat $tmpdir/*.target $tmpdir/*.control |");
    while(<IN>){
	chomp;
	@row = split;
	if ($row[0] =~ /^[0-9]+$/){
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
	$row[0] =~ s/^0+//g;
	$row[1] += 0;
	$dat = join("\t", @row[0 .. $#row -1]);
	if($dat ne $prev){
	    if ($cw >= 5 and $cm <= 1){
		if ($tm >= 5 and $tw <= 1){
		    $genotype = 'M';
		}elsif(($tm + $tw) > 0 and $tm / ($tm + $tw) > 0.3 and $tm / ($tm + $tw) < 0.7 and $cm <= 1 and $tw >= 5){
		    $genotype = "H";
		}else{
		    $genotype = "R";
		}
	    }else{
		$genotype = "N";
	    }
	    $_ = $prev;
	    @prev = split;
	    if ($prev[0] =~ /^[0-9]+$/){
		$prev[0] = "000" . $prev[0];
		$prev[0] = substr($prev[0], length($prev[0]) -3, 3);
	    }
	    $prev[1] = "000000000000" . $prev[1];
	    $prev[1] = substr($prev[1], length($prev[1]) -11, 11);
	    print OUT "$prev[0]:$prev[1]\t$cw\t$cm\t$tw\t$tm\t$genotype\n" if $prev ne "";
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
    if ($cw >= 5 and $cm <= 1){
	if ($tm >= 5 and $tw <= 1){
	    $genotype = 'M';
	}elsif(($tm + $tw) > 0 and $tm / ($tm + $tw) > 0.3 and $tm / ($tm + $tw) < 0.7 and $cm <= 1 and $tw >= 5){
	    $genotype = "H";
	}else{
	    $genotype = "R";
	}
    }else{
	$genotype = "N";
    }
    $_ = $prev;
    @prev = split;
    if ($prev[0] =~ /^[0-9]+$/){
	$prev[0] = "000" . $prev[0];
	$prev[0] = substr($prev[0], length($prev[0]) -3, 3);
    }
    $prev[1] = "000000000000" . $prev[1];
    $prev[1] = substr($prev[1], length($prev[1]) -11, 11);
    print OUT "$prev[0]:$prev[1]\t$cw\t$cm\t$tw\t$tm\t$genotype\n";
    close(OUT);
    
    open(IN, "join $tmpdir/$target.map $tmpdir/$target.verify|");
    open(OUT, "> $wd/$target/$target.kmer.snp");
    while(<IN>){
	chomp;
	@row = split;
	($chr, $pos) = split(':', $row[0]);
	$chr =~ s/^0+//g;
	$pos += 0;
	print OUT "$chr\t$pos\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\t$row[9]\t$row[10]\t$row[11]\t$row[12]\t$row[13]\t$row[14]\t$row[15]\t$row[16]\t$row[17]\t$row[18]\n";
    }
    close(IN);
    close(OUT);
}

sub snpReadCount{
    my (@row, $dat, %count, $cm, $cw, $tm, $tw, $count, @prev, $prev, $genotype);
    open(OUT, "| sort $sort_opt -T $tmpdir | uniq -c > $tmpdir/bi.snp");
    open(IN, "cat $tmpdir/*.target $tmpdir/*.control |");
    while(<IN>){
	chomp;
	@row = split;
	if ($row[0] =~ /^[0-9]+$/){
	    $row[0] = "000" . $row[0];
	    $row[0] = substr($row[0], length($row[0]) -3, 3);
	}
	$row[1] = "000000000000" . $row[1];
	$row[1] = substr($row[1], length($row[1]) -11, 11);
	$dat = join(":", @row);
	print OUT "$dat\n";
    }
    close(IN);
    close(OUT);
    $cm = 0;
    $cw = 0;
    $tm = 0;
    $tw = 0;
    &waitFile("$tmpdir/bi.snp");
    open(OUT, "> $wd/$target/$target.bi.snp");
    open(IN, "$tmpdir/bi.snp");
    while(<IN>){
	chomp;
	($count, $dat) = split;
	@row = split(':', $dat);
	$row[0] =~ s/^0+//g;
	$row[1] += 0;
	$dat = join("\t", @row[0 .. $#row -1]);
	if($dat ne $prev){
	    if ($cw >= 5 and $cm <= 1){
		if ($tm >= 5 and $tw <= 1){
		    $genotype = 'M';
		}elsif(($tm + $tw) > 0 and $tm / ($tm + $tw) > 0.3 and $tm / ($tm + $tw) < 0.7 and $cm <= 1 and $tw >= 5){
		    $genotype = "H";
		}else{
		    $genotype = "R";
		}
	    }else{
		$genotype = "N";
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
    if ($cw >= 5 and $cm <= 1){
	if ($tm >= 5 and $tw <= 1){
	    $genotype = 'M';
	}elsif(($tm + $tw) > 0 and $tm / ($tm + $tw) > 0.3 and $tm / ($tm + $tw) < 0.7 and $cm <= 1 and $tw >= 5){
	    $genotype = "H";
	}else{
	    $genotype = "R";
	}
    }else{
	$genotype = "N";
    }
    print OUT "$prev\t$cw\t$cm\t$tw\t$tm\t$genotype\n";
    close(OUT);
    system("rm $tmpdir/bi.snp");
}

sub joinControlFunc{
    my $tag = shift;
    report("Selecting sequence data for verify. $tag");
    if (-e "$wd/$control/sort_uniq/$control.$tag.gz"){
	system("bash -c \"join <($zcat $wd/$control/sort_uniq/$control.$tag.gz) <($zcat $tmpdir/$tag.gz)$tmpdir/$tag | cut -d ' ' -f 2- > $tmpdir/$tag.control\"");
    }else{
	system("bash -c \"join <($zcat $wd/$control/sort_uniq/$control.sort_uniq.$tag.gz) <($zcat $tmpdir/$tag.gz) | cut -d ' ' -f 2- > $tmpdir/$tag.control\"");
    }
    &waitFile("$tmpdir/$tag.control");
    system("rm $tmpdir/$tag.gz");
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
    &joinAll;
}

sub joinTargetFunc{
    my $tag = shift;
    report("Selecting sequence data for verify. $tag");
    system("bash -c \"join <($zcat $wd/$target/sort_uniq/$target.sort_uniq.$tag.gz) <($zcat $tmpdir/$tag.gz) | cut -d ' ' -f 2- > $tmpdir/$tag.target\"");
    &waitFile("$tmpdir/$tag.target");
    system("rm $tmpdir/$tag.gz");
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
    &joinAll;
}

sub sortSeqFunc{
    my $tag = shift;
    my ($nuca, $nucb, $nucc, $subtag, $fin, $fout);
    report("Sorting sequence data for verify. $tag");
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		$fout = $tag . $subtag;
		open($fout, "> $tmpdir/tmp.$fout");
	    }
	}
    }
    $fin = $tag;
    open($fin, "$zcat $tmpdir/$tag.tmp.gz |");
    while(<$fin>){
	$subtag = substr($_, 0, 6);
	print $subtag $_;
    }
    close($fin);
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		$fout = $tag . $subtag;
		close($fout);
		system("sort $sort_opt -T $tmpdir $tmpdir/tmp.$fout > $tmpdir/tmpsort.$fout && rm $tmpdir/tmp.$fout");
	    }
	}
    }
    system("cat $tmpdir/tmpsort.$tag* |gzip > $tmpdir/$tag.gz");
    system("rm $tmpdir/$tag.tmp.gz $tmpdir/tmpsort.$tag*");
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
    &joinAll;
}

sub svMkC{
    report("Making control data for verification of sv");
    my ($nuca, $nucb, $nucc, $tag, $hchr, $hpos, $tchr, $tpos, $direction, $type, $size, @row, $current, $prev, $prev_hchr, $posa, $posb, $inside, $head, $tail, $ref_seq, $mut_seq, $slength, $cm, $cw, $hchr_seq);

    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "|gzip > $tmpdir/$tag.tmp.gz");
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
	$hpos += 0;
	$current = "$hchr $hpos $tchr $tpos";
	next if $current eq $prev;
	if ($hchr ne $prev_hchr){
	    my $chr_file = "$refdir/chr$hchr";
	    open (HIN, $chr_file);
	}
	$prev = $current;
	$prev_hchr = $hchr;
	$posa = length($row[2]);
	$posb = length($row[8]);
	if ($posa < $posb){
	    $inside = substr($row[5], $posa -1, $posb - $posa + 1);
	    seek(HIN, $hpos - $clength + ($posb - $posa -1), 0);
	    next if $clength - ($posb - $posa -1) -1 < 0;
	    read(HIN, $head, $clength - ($posb - $posa -1) -1);
	    next if length($head) != $clength - ($posb - $posa -1) -1;
	    if ($direction eq "f"){
		if ($tchr eq $hchr){
		    seek(HIN, $tpos, 0);
		    next if $clength - ($posb - $posa -1) -1 < 0;
		    read(HIN, $tail, $clength - ($posb - $posa -1) -1);
		}else{
		    open(IN, "$refdir/chr$tchr");
		    seek(IN, $tpos, 0);
		    next if $clength - ($posb - $posa -1) -1 < 0;
		    read(IN, $tail, $clength - ($posb - $posa -1) -1);
		    close(IN);
		}
		next if length($tail) != $clength - ($posb - $posa -1) -1;
	    }else{
		if ($tchr eq $hchr){
		    seek(HIN, $tpos - $clength + ($posb - $posa -1), 0);
		    next if $clength - ($posb - $posa -1) -1 < 0;
		    read(HIN, $tail, $clength - ($posb - $posa -1) -1);
		    $tail = &complement($tail);
		}else{
		    open(IN, "$refdir/chr$tchr");
		    seek(IN, $tpos - $clength + ($posb - $posa -1), 0);
		    next if $clength - ($posb - $posa -1) -1 < 0;
		    read(IN, $tail, $clength - ($posb - $posa -1) -1);
		    close(IN);
		    $tail = &complement($tail);
		}
		next if length($tail) != $clength - ($posb - $posa -1) -1;
	    }
	    seek(HIN, $hpos - $clength, 0);
	    read(HIN, $ref_seq, $clength * 2 -1);
	    next if length($ref_seq) != $clength * 2 -1;
	    $mut_seq = $head . $inside . $tail;
	}elsif($posa == $posb){
	    $inside = substr($row[5], $posa -1, $posb - $posa + 1);
	    seek(HIN, $hpos - $clength + ($posb - $posa), 0);
	    next if $clength - ($posb - $posa -1) -2 < 0;
	    read(HIN, $head, $clength - ($posb - $posa -1) -2);
	    next if length($head) != $clength - ($posb - $posa -1) -2;
	    if ($direction eq "f"){
		if ($tchr eq $hchr){
		    seek(HIN, $tpos, 0);
		    next if $clength - ($posb - $posa) -1 < 0;
		    read(HIN, $tail, $clength - ($posb - $posa) -1);
		}else{
		    open(IN, "$refdir/chr$tchr");
		    seek(IN, $tpos, 0);
		    next if $clength - ($posb - $posa) -1 < 0;
		    read(IN, $tail, $clength - ($posb - $posa) -1);
		    close(IN);
		}
		next if length($tail) != $clength - ($posb - $posa) -1;
	    }else{
		if ($tchr eq $hchr){
		    seek(HIN, $tpos- $clength + ($posb - $posa), 0);
		    next if $clength - ($posb - $posa) -1 < 0;
		    read(HIN, $tail, $clength - ($posb - $posa) -1);
		    $tail = &complement($tail);
		}else{
		    open(IN, "$refdir/chr$tchr");
		    seek(IN, $tpos - $clength + ($posb - $posa), 0);
		    next if $clength - ($posb - $posa) -1;
		    read(IN, $tail, $clength - ($posb - $posa) -1);
		    close(IN);
		    $tail = &complement($tail);
		}
		next if length($tail) != $clength - ($posb - $posa) -1;
	    }
	    seek(HIN, $hpos - $clength, 0);
	    read(HIN, $ref_seq, $clength * 2 -1);
	    next if length($ref_seq) != $clength * 2 -1;
	    $mut_seq = $head . $inside . $tail;
	}else{
	    $inside = substr($row[5], $posb, $posa - $posb -1);
	    seek(HIN, $hpos - $clength, 0);
	    next if $clength - ($posa - $posb -1) -1 < 0;
	    read(HIN, $head, $clength - ($posa - $posb -1) -1);
	    next if length($head) != $clength - ($posa - $posb -1) -1;
	    if ($direction eq "f"){
		if ($tchr eq $hchr){
		    seek(HIN, $tpos + ($posa - $posb -1), 0);
		    next if $clength - ($posa - $posb -1) -1 < 0;
		    read(HIN, $tail, $clength - ($posa - $posb -1) -1);
		}else{
		    open(IN, "$refdir/chr$tchr");
		    seek(IN, $tpos + ($posa - $posb -1), 0);
		    next if $clength - ($posa - $posb -1) -1 < 0;
		    read(IN, $tail, $clength - ($posa - $posb -1) -1);
		    close(IN);
		}
	    }else{
		if ($tchr eq $hchr){
		    seek(HIN, $tpos - $clength, 0);
		    next if $clength - ($posa - $posb -1) -1 < 0;
		    read(HIN, $tail, $clength - ($posa - $posb -1) -1);
		    $tail = &complement($tail);
		}else{
		    open(IN, "$refdir/chr$tchr");
		    seek(IN, $tpos - $clength, 0);
		    next if  $clength - ($posa - $posb -1) -1 < 0;
		    read(IN, $tail, $clength - ($posa - $posb -1) -1);
		    close(IN);
		    $tail = &complement($tail);
		}
	    }
	    next if length($tail) != $clength - ($posa - $posb -1) -1;
	    seek(HIN, $hpos - $clength - ($posa - $posb), 0);
	    read(HIN, $ref_seq, $clength * 2 -1);
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
    report("Making control data for verification of sv complete");
}

sub svMkT{
    report("Making target data for verification of sv");
    my ($nuca, $nucb, $nucc, $tag, $hchr, $hpos, $tchr, $tpos, $direction, $type, $size, @row, $current, $prev, $prev_hchr, $posa, $posb, $inside, $head, $tail, $ref_seq, $mut_seq, $slength, $tm, $tw, $hchr_seq);

    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "|gzip > $tmpdir/$tag.tmp.gz");
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
	$hpos += 0;
	$current = "$hchr $hpos $tchr $tpos";
	next if $current eq $prev;
	if ($hchr ne $prev_hchr){
	    my $chr_file = "$refdir/chr$hchr";
	    open (HIN, $chr_file);
	}
	$prev = $current;
	$prev_hchr = $hchr;
	$posa = length($row[2]);
	$posb = length($row[8]);
	if ($posa < $posb){
	    $inside = substr($row[5], $posa -1, $posb - $posa + 1);
	    seek(HIN, $hpos - $length + ($posb - $posa -1), 0);
	    next if $length - ($posb - $posa -1) -1 < 0;
	    read(HIN, $head, $length - ($posb - $posa -1) -1);
	    next if length($head) != $length - ($posb - $posa -1) -1;
	    if ($direction eq "f"){
		if ($tchr eq $hchr){
		    seek(HIN, $tpos, 0);
		    next if $length - ($posb - $posa -1) -1 < 0;
		    read(HIN, $tail, $length - ($posb - $posa -1) -1);
		}else{
		    open(IN, "$refdir/chr$tchr");
		    seek(IN, $tpos, 0);
		    next if $length - ($posb - $posa -1) -1 < 0;
		    read(IN, $tail, $length - ($posb - $posa -1) -1);
		    close(IN);
		}
		next if length($tail) != $length - ($posb - $posa -1) -1;
	    }else{
		if ($tchr eq $hchr){
		    seek(HIN, $tpos - $length + ($posb - $posa -1), 0);
		    next if $length - ($posb - $posa -1) -1 < 0;
		    read(HIN, $tail, $length - ($posb - $posa -1) -1);
		    $tail = &complement($tail);
		}else{
		    open(IN, "$refdir/chr$tchr");
		    seek(IN, $tpos - $length + ($posb - $posa -1), 0);
		    next if $length - ($posb - $posa -1) -1 < 0;
		    read(IN, $tail, $length - ($posb - $posa -1) -1);
		    close(IN);
		    $tail = &complement($tail);
		}
		next if length($tail) != $length - ($posb - $posa -1) -1;
	    }
	    seek(HIN, $hpos - $length, 0);
	    read(HIN, $ref_seq, $length * 2 -1);
	    next if length($ref_seq) != $length * 2 -1;
	    $mut_seq = $head . $inside . $tail;
	}elsif($posa == $posb){
	    $inside = substr($row[5], $posa -1, $posb - $posa + 1);
	    seek(HIN, $hpos - $length + ($posb - $posa), 0);
	    next if $length - ($posb - $posa -1) -2 < 0;
	    read(HIN, $head, $length - ($posb - $posa -1) -2);
	    next if length($head) != $length - ($posb - $posa -1) -2;
	    if ($direction eq "f"){
		if ($tchr eq $hchr){
		    seek(HIN, $tpos, 0);
		    next if $length - ($posb - $posa) -1 < 0;
		    read(HIN, $tail, $length - ($posb - $posa) -1);
		}else{
		    open(IN, "$refdir/chr$tchr");
		    seek(IN, $tpos, 0);
		    next if $length - ($posb - $posa) -1 < 0;
		    read(IN, $tail, $length - ($posb - $posa) -1);
		    close(IN);
		}
		next if length($tail) != $length - ($posb - $posa) -1;
	    }else{
		if ($tchr eq $hchr){
		    seek(HIN, $tpos- $length + ($posb - $posa), 0);
		    next if $length - ($posb - $posa) -1 < 0;
		    read(HIN, $tail, $length - ($posb - $posa) -1);
		    $tail = &complement($tail);
		}else{
		    open(IN, "$refdir/chr$tchr");
		    seek(IN, $tpos - $length + ($posb - $posa), 0);
		    next if $length - ($posb - $posa) -1 < 0;
		    read(IN, $tail, $length - ($posb - $posa) -1);
		    close(IN);
		    $tail = &complement($tail);
		}
		next if length($tail) != $length - ($posb - $posa) -1;
	    }
	    seek(HIN, $hpos - $length, 0);
	    read(HIN, $ref_seq, $length * 2 -1);
	    next if length($ref_seq) != $length * 2 -1;
	    $mut_seq = $head . $inside . $tail;
	}else{
	    $inside = substr($row[5], $posb, $posa - $posb -1);
	    seek(HIN, $hpos - $length, 0);
	    next if $length - ($posa - $posb -1) -1 < 0;
	    read(HIN, $head, $length - ($posa - $posb -1) -1);
	    next if length($head) != $length - ($posa - $posb -1) -1;
	    if ($direction eq "f"){
		if ($tchr eq $hchr){
		    seek(HIN, $tpos + ($posa - $posb -1), 0);
		    next if $length - ($posa - $posb -1) -1 < 0;
		    read(HIN, $tail, $length - ($posa - $posb -1) -1);
		}else{
		    open(IN, "$refdir/chr$tchr");
		    seek(IN, $tpos + ($posa - $posb -1), 0);
		    next if $length - ($posa - $posb -1) -1 < 0;
		    read(IN, $tail, $length - ($posa - $posb -1) -1);
		    close(IN);
		}
	    }else{
		if ($tchr eq $hchr){
		    seek(HIN, $tpos - $length, 0);
		    next if $length - ($posa - $posb -1) -1 < 0;
		    read(HIN, $tail, $length - ($posa - $posb -1) -1);
		    $tail = &complement($tail);
		}else{
		    open(IN, "$refdir/chr$tchr");
		    seek(IN, $tpos - $length, 0);
		    next if $length - ($posa - $posb -1) -1 < 0;
		    read(IN, $tail, $length - ($posa - $posb -1) -1);
		    close(IN);
		    $tail = &complement($tail);
		}
	    }
	    next if length($tail) != $length - ($posa - $posb -1) -1;
	    seek(HIN, $hpos - $length - ($posa - $posb), 0);
	    read(HIN, $ref_seq, $length * 2 -1);
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
    report("Making target data for verification of sv complete");
}
    
sub svSort{
    report("Making $target.sv.sort");
    my ($flag, @row, $chr, $out, $count);
    my $fin = "in-sv";
    my $fout = "out-sv";
    open($fin, "$wd/$target/$target.aln");
    open($fout, "|sort $sort_opt -T $tmpdir |uniq > $tmpdir/$target.sv.sort");
    while(<$fin>){
	chomp;
	if (/ion/){
	    $flag = 1;
	    @row = split;
	    if ($row[1] =~ /^[0-9]+$/){
		$chr = "000$row[1]";
		$chr = substr($chr, length($chr) - 3, 3);
	    }else{
		$chr = $row[1];
	    }
	    $pos = "00000000000" . $row[2];
	    $pos = substr($pos, length($pos) - 11, 11);
	    if ($row[7] eq ""){
            $out .= "$row[0] $chr $pos $row[3] $row[4] $row[5] $row[6]\t";
        }else{
            $out .= "$row[0] $chr $pos $row[3] $row[4] $row[5] $row[6] $row[7]\t";
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
	    if ($row[1] =~ /^[0-9]+$/){
		$row[1] = "000$row[1]";
		$row[1] = substr($row[1], length($row[1]) -3, 3);
	    }
	    $row[2] = "0000000000$row[2]";
	    $row[2] = substr($row[2], length($row[2]) - 11, 11);
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
    &joinAll;
}

sub sortData4MapFunc{
    my $tag = shift;
    my ($nuca, $nucb, $nucc, $subtag, $fin, $ftmp, $fout);
    report("Sorting sequence data. $tag");
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		$ftmp = $tag . $subtag;
		open($ftmp, "> $tmpdir/tmp.$ftmp");
	    }
	}
    }
    $fin = $tag;
    open($fin, "cat $tmpdir/$tag.tmp.* |");
    while(<$fin>){
	$ftmp = substr($_, 0, 6);
	print $ftmp $_;
    }
    close($fin);
    $fout = "tmpout.$tag";
    open($fout, "|gzip -f > $tmpdir/$tag.gz");
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		$ftmp = $tag . $subtag;
		close($ftmp);
		$fin = "tmp.$tag";
		open($fin, "sort $sort_opt -T $tmpdir $tmpdir/tmp.$ftmp |");
		while(<$fin>){
		    print $fout $_;
		}
		close($fin);
		system("rm $tmpdir/tmp.$ftmp");
	    }
	}
    }
    close($fout);
    &waitFile("$tmpdir/$tag.gz");
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
		report("Bidirectional alignment: Making data for first mapping. $tag");
		threads->new(\&mkData4MapFFunc, $tag);
	    }
	}
    }
    &joinAll;
}

sub mkData4MapFFunc{
    my $tag = shift;
    my ($nuca, $nucb, $nucc, $subtag, $margin, $head_pos, $head, $fout);
    my $fin = $tag . "IN";
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		$fout = $subtag . $tag;
		open($fout, "> $tmpdir/$subtag.tmp.$tag")
	    }
	}
    }
    open($fin, "$zcat $wd/$target/sort_uniq/$target.sort_uniq.$tag.gz |");
    while(<$fin>){
	chomp;
	foreach $margin ('0', '5', '10', '15'){
	    $head_pos = $margin + 1;
	    $head = substr($_, $head_pos -1, 20);
	    $subtag = substr($head, 0, 3);
	    $fout = $subtag . $tag;
	    print $fout "$head $_ $head_pos\n";
	}
    }
    close($fin);  
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		$fout = $subtag . $tag;
		close($fout)
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
    &joinAll;
}

sub mkData4MapRFunc{
    my $tag = shift;
    my ($nuca, $nucb, $nucc, $subtag, $margin, $tail_pos, $tail, $fout);
    my $fin = $tag . "IN";
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag =  $nuca . $nucb . $nucc;
		$fout = $subtag . $tag;
		open($fout, "> $tmpdir/$subtag.tmp.$tag")
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
	$fout = $subtag . $tag;
	print $fout "$tail $_ $tail_pos\n";
    }
    close($fin);  
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$subtag = $nuca . $nucb . $nucc;
		$fout = $subtag . $tag;
		close($fout);
	    }
	}
    }
    system("rm $tmpdir/$tag.map");
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
    &joinAll;
}

sub mapFunc{
    my $tag = shift;
    system("bash -c \"join <($zcat $tmpdir/$tag.gz) <($zcat $refdir/ref20_uniq.$tag.gz) |cut -c 22- > $tmpdir/$tag.map\"");
    &waitFile("$tmpdir/$tag.map");
    system("rm $tmpdir/$tag.gz");
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
    &joinAll;
    system("cat $tmpdir/$target.aln.* > $wd/$target/$target.aln && rm $tmpdir/*");
}

sub alignFunc{
    my $tag = shift;
    my ($seq, $hpos, $hchr, $head_pos, $head_direction, $tpos, $tchr, $tail_pos, $tail_direction, $length, $head, $tail, $hhit, $thit, $margin, $head_seq, $distance, $head_bar, $tail_bar, $head_space, $tail_space, @head, @tail, @seq, $mcount, $out, $i, $j, $k, $head_junction, $head_fail, $unmatch, $hcount, $head_junction, $tail_junction, $tail_fail, $tcount, $tail_direction, $type, $distance, $fin, $fout, $chr_file, $nflag);
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
	$nflag = 0;
	$margin = $hpos -1;
	
	open($chr_file, "$wd/$ref/chr$hchr");
	if ($head_direction eq "f"){
	    seek($chr_file, $head_pos - $margin -1, 0);
	    read($chr_file, $head_seq, $length);
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
			    $nflag = 1 if $head[$i] eq "N";
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
	    }else{
		seek($chr_file, $tail_pos - 20 -$margin, 0);
		read($chr_file, $tail_seq, $length);
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
		    $nflag = 1 if $head[$i] eq "N";
		    $head_bar .= " ";
		}
		if ($seq[$i] eq $tail[$i]){
		    $tail_bar .= "|";
		    $tail_bar[$i] = "|";
		}else{
		    $nflag = 1 if $tail[$i] eq "N";
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
	    if ($head_junction ne "" and $tail_junction ne "" and $nflag == 0){
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
    close($chr_file);
    $semaphore->up;
    system("rm $tmpdir/$tag.map");
}

sub joinAll{
    foreach my $thr (threads->list){
	$thr->join;
    }
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
	sleep 1;
    }
    while(1){
	my $mtime = (stat($file))[9];
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
