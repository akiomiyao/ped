#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     qsub_kmer.pl - pipeline of kmer method for computer cluster. 

     run  perl qsub_kmer.pl target control reference
     e.g. perl qsub_kmer.pl SRR8181712 default TAIR10

     This script is for computer cluster.

     Before run the script, downloading target reads and reference genome data
     are required.

     To make the reference data,
     run  perl mkref.pl reference
     e.g. perl mkref.pl TAIR10

     Run mkref.pl without reference name, references in config file will be listed.
     You can select reference from the list.

     To download from SRA in NCBI,
     run  perl download.pl accession
     e.g. perl download.pl SRR8181712

     To analyze your data of short reads,
     mkdir mydata
     mkdir mydata/read
     and the copy your fastq data to mydata/read
     and then run perl qsub_kmer.pl mydata myconrtol reference

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[1] eq ""){
    print $usage;
    exit;
}

$target  = $ARGV[0];
$control = $ARGV[1];
$ref     = $ARGV[2];
$tmpdir  = $ARGV[3];
$control = $ref if $control eq "default";
 
open(IN, "config");
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

report("kmer.pl start.");

if (! -e "$target/$target.sort_uniq"){
    report("Making $target.sort_uniq.");
    $qsub = "-v target=$target sort_uniq.pl";
    @job = ();
    &doQsub($qsub);
}

if (! -e "$control/$control.sort_uniq"){
    report("Making $control.sort_uniq.");
    $qsub = "-v target=$control sort_uniq.pl";
    @job = ();
    &doQsub($qsub);
}
&holdUntilJobEnd;

if (! -e "$control/$control.lbc.AAA.gz"){
    report("Making kmer count of $control.");
    $qsub = "-v target=$control split_sort_uniq.pl";
    &doQsub($qsub);
    &holdUntilJobEnd;
    @job = ();
    opendir(DIR, "$control/tmp");
    foreach (sort readdir(DIR)){
	next if !/$control/;
	@row = split('\.', $_);
	if ($tmpdir ne ""){
	    $qsub = "-v target=$row[0],number=$row[2],tmpdir=$tmpdir count.pl";
	}else{
	    $qsub = "-v target=$row[0],number=$row[2] count.pl";
	}
	&doQsub($qsub);
    }
    &holdUntilJobEnd;
    @job = ();

    opendir(DIR, "$control");
    foreach (sort readdir(DIR)){
	@row = split('\.', $_);
	if ($row[1] eq "count"){
	    $tag{$row[2]} = 1;
	}
    }
    foreach $tag (sort keys %tag){
	if ($tmpdir ne ""){
	    $qsub = "-v target=$control,tag=$tag,tmpdir=$tmpdir merge.pl";
	}else{
	    $qsub = "-v target=$control,tag=$tag merge.pl";
	}
	&doQsub($qsub);
    }
    &holdUntilJobEnd;
    @job = ();
    foreach $tag (sort keys %tag){
	if ($tmpdir ne ""){
	    $qsub = "-v target=$control,tag=$tag,tmpdir=$tmpdir last_base_count.pl";
	}else{
	    $qsub = "-v target=$control,tag=$tag last_base_count.pl";
	}
	&doQsub($qsub);
    }
    &holdUntilJobEnd;
    @job = ();
}

report("Making kmer count of $target.");
$qsub = "-v target=$target split_sort_uniq.pl";
&doQsub($qsub);
&holdUntilJobEnd;

@job = ();
opendir(DIR, "$target/tmp");
foreach (sort readdir(DIR)){
    next if ! /$target/;
    @row = split('\.', $_);
    if ($tmpdir ne ""){
	$qsub = "-v target=$target,number=$row[2],tmpdir=$tmpdir count.pl";
    }else{
	$qsub = "-v target=$target,number=$row[2] count.pl";
    }
    &doQsub($qsub);
}
&holdUntilJobEnd;

@job = ();
report("Marge of kmer count.");
opendir(DIR, "$target");
foreach (sort readdir(DIR)){
    @row = split('\.', $_);
    if ($row[1] eq "count"){
	$tag{$row[2]} = 1;
    }
}
close(DIR);

foreach $tag (sort keys %tag){
    if ($tmpdir ne ""){
	$qsub = "-v target=$target,tag=$tag,tmpdir=$tmpdir merge.pl";
    }else{
	$qsub = "-v target=$target,tag=$tag merge.pl";
    }
    &doQsub($qsub);
}
&holdUntilJobEnd;

@job = ();
foreach $tag (sort keys %tag){
    if ($tmpdir ne ""){
	$qsub = "-v target=$target,tag=$tag,tmpdir=$tmpdir last_base_count.pl";
    }else{
	$qsub = "-v target=$target,tag=$tag last_base_count.pl";
    }
    &doQsub($qsub);
}
&holdUntilJobEnd;

report("Making snp");
@job = ();
opendir(DIR, $target);
foreach $tag (sort keys %tag){
    if ($tmpdir ne ""){
	$qsub = "-v target=$target,control=$control,tag=$tag,tmpdir=$tmpdir snp.pl";
    }else{
	$qsub = "-v target=$target,control=$control,tag=$tag snp.pl";
    }
    &doQsub($qsub);
}
&holdUntilJobEnd;

system("cat $target.snp.??? > $target.kmer");

report("Making map");
@job = ();
opendir(DIR, $target);
foreach $tag (sort keys %tag){
    if ($tmpdir ne ""){
	$qsub = "-v target=$target,ref=$ref,tag=$tag,tmpdir=$tmpdir map.pl";
    }else{
	$qsub = "-v target=$target,ref=$ref,tag=$tag map.pl";
    }
    &doQsub($qsub);
}
&holdUntilJobEnd;
report("Verifying");
@job = ();
opendir(DIR, $target);
foreach $tag (sort keys %tag){
    if ($tmpdir ne ""){
	$qsub = "-v target=$target,control=$control,ref=$ref,number=$tag,type=kmer,tmpdir=$tmpdir verify_snp.pl";
    }else{
	$qsub = "-v target=$target,control=$control,ref=$ref,number=$tag,type=kmer verify_snp.pl";
    }
    &doQsub($qsub);
}
&holdUntilJobEnd;

open(IN, "cat $target/$target.kmer.verify.*| sort -T $target |");
while(<IN>){
    @row = split;
    if ($row[0] ne $prev){
	close(OUT);
	open(OUT, "| sort -T $target -k 2 -n > $target/$target.kmer_chr.$row[0]");
    }
    print OUT;
    $prev = $row[0];
}
close(IN);
close(OUT);

open(OUT, "> $target/$target.kmer.snp");
foreach $chr (@chr){
    open(IN, "$target/$target.kmer_chr.$chr");
    while(<IN>){
	print OUT;
    }
    close(IN);
}
close(OUT);

report("Convert to vcf");
system("rm $target/$target.count.* $target/$target.snp.* $target/$target.map.* $target/$target.kmer.verify.* $target/$target.kmer_chr.* $target/$target.lbc.* ");
system("perl snp2vcf.pl $target kmer");
report("kmer.pl done.");


sub report{
    my $message = shift;
    my $now = `date`;
    chomp($now);
    print "$now $message\n";
}

sub doQsub{
    my $qsub = shift;
    my $job;
    open(IN, "qsub $qsub |");
    $job = <IN>;
    close(IN);
    chomp($job);
    push(@job, $job);
    sleep 1;
}

sub holdUntilJobEnd{
    my (@row, %stat, $job, $flag);
    while(1){
	$flag = 0;
	open(IN, "qstat |");
	while(<IN>){
	    @row = split;
	    $stat{$row[0]} = $row[4];
	}
	close(IN);
	foreach $job(@job){
	    $flag = 1 if $stat{$job} ne "C";
	}
	last if ! $flag;
	sleep 10;
    }
}
