#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     qsub_bidirectional.pl - pipeline for bidirectional method for computer cluster. 

     run  perl qsub_bidirectional.pl target control reference tmpdir
     e.g. perl qsub_bidirectional.pl SRR8181712 default TAIR10 /mnt/ssd

     if you want to compare between target and control samples,
     specify control, e.g., SRR1581142.
     Otherwise, set default for control, reference data will be used as contol.

     This script is for coumputer cluster.

     Before run the script, downloading target reads and reference genome data
     are required.

     To make the reference data,
     run  perl mkref.pl reference
     e.g. perl mkref.pl TAIR10
     or   qsub -v target=TAIR10 mkref.pl
 
     Run without reference name, references in config file will be listed.
     You can select reference from the list.

     To download from SRA in NCBI,
     run  perl download.pl accession
     e.g. perl download.pl SRR8181712

     To analyze your data of short reads,
     mkdir mydata
     mkdir mydata/read
     and the copy your fastq data to mydata/read
     and then run perl qsub_bidirectional.pl mydata deault reference
               or perl qsub_bidirectional.pl mydata mycontrol reference

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

report("Aligning of $target sequence to $ref genome.");
@job = ();
foreach $margin (0, 5, 10, 15){
    if ($tmpdir ne ""){
	$qsub = "-v target=$target,ref=$ref,margin=$margin,tmpdir=$tmpdir align.pl";
    }else{
	$qsub = "-v target=$target,ref=$ref,margin=$margin align.pl";
    }
    &doQsub($qsub);
}
&holdUntilJobEnd;

@job = ();
report("Splitting of Indel alignment. ");
$qsub = "-v target=$target,ref=$ref split_indel.pl";
&doQsub($qsub);

report("Splitting of SNP alignment.");
$qsub = "-v target=$target,ref=$ref split_snp.pl";
&doQsub($qsub);
&holdUntilJobEnd;

@job = ();
opendir(DIR, $target);
foreach(sort readdir(DIR)){
    if (/$target.indel.[0-9][0-9]$/){
	$number = (split('\.', $_))[2];
	report("Verifying indel $number");
	if ($tmpdir ne ""){
	    $qsub = "-v target=$target,control=$control,ref=$ref,number=$number,type=bi,tmpdir=$tmpdir verify_indel.pl";
	}else{
	    $qsub = "-v target=$target,control=$control,ref=$ref,number=$number,type=bi verify_indel.pl";
	}
	&doQsub($qsub);
    }elsif(/$target.snp.[0-9][0-9]$/){
	$number = (split('\.', $_))[2];
	report("Verifying SNPs $number");
	if ($tmpdir ne ""){
	    $qsub = "-v target=$target,control=$control,ref=$ref,number=$number,type=bi,tmpdir=$tmpdir verify_snp.pl";
	}else{
	    $qsub = "-v target=$target,control=$control,ref=$ref,number=$number,type=bi verify_snp.pl";
	}
	&doQsub($qsub);
    }
}
closedir(DIR);
&holdUntilJobEnd;

report("Making vcf file of SNP");
system("perl snp2vcf.pl $target");
system("cat $target/$target.indel.verify.* > $target/$target.indel && rm $target/$target.indel.verify.* $target/$target.indel.??");
system("cat $target/$target.snp.verify.* > $target/$target.snp && rm $target/$target.snp.verify.* $target/$target.snp.??");
report("bidirectional.pl complete.");

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
