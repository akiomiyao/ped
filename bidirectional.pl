#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     bidirectional.pl - pipeline for bidirectional method for stand alone macnine. 

     run  perl bidirectional.pl accession control reference
     e.g. perl bidirectional.pl SRR8181712 default TAIR10

     This script is for stand alone machine.
     If you want to detect polymorphisms between target and control (e.g.SRR1581142),
     set SRR1581142 as the control.
     If you want to use the control, downloading control sequeces and making the sort_uniq data are required, before the run the kmer.pl.
     Before run the script, downloading target reads and reference genome data
     are required.
       For example,
       % perl download.pl SRR1581142
       % perl sort_uniq.pl SRR1581142

     To make the reference data,
     run  perl mkref.pl reference
     e.g. perl mkref.pl TAIR10

     Run without reference name, references in config file will be listed.
     You can select reference from the list.

     To download from SRA in NCBI,
     run  perl download.pl accession
     e.g. perl download.pl SRR8181712

     To analyze your data of short reads,
     mkdir mydata
     mkdir mydata/read
     and the copy your fastq data to mydata/read
     and then run perl bidirectional.pl mydata reference

Author: Akio Miyao <miyao@affrc.go.jp>

';

if ($ARGV[1] eq ""){
    print $usage;
    exit;
}

$target  = $ARGV[0];
$control = $ARGV[1];
$ref     = $ARGV[2];

$control = $ref if $control eq "default";

if (! -e "$target/$target.sort_uniq"){
    report("Making $target.sort_uniq.");
    system("perl sort_uniq.pl $target");
}

&sortWait("$target/$target.sort_uniq");

if (! -e "$control/$control.sort_uniq"){
    report("Making $control.sort_uniq.");
    system("perl sort_uniq.pl $control");
}

&sortWait("$control/$control.sort_uniq");

report("Aligning of $target sequence to $ref genome. margin = 5");
system("perl align.pl $target $ref 5");

## If you want speed up rather than sensitivity, comment out with '#'
## system command for perl align.pl.

report("Aligning of $target sequence to $ref genome. margin = 0");
system("perl align.pl $target $ref 0");
report("Aligning of $target sequence to $ref genome. margin = 10");
system("perl align.pl $target $ref 10");
report("Aligning of $target sequence to $ref genome. margin = 15");
system("perl align.pl $target $ref 15");

report("Splitting of Indel alignment.");
system("perl split_indel.pl $target $ref");
report("Splitting of SNP alignment.");
system("perl split_snp.pl $target $ref");

opendir(DIR, $target);
foreach(sort readdir(DIR)){
    if (/$target.indel.[0-9][0-9]$/){
	$number = (split('\.', $_))[2];
	report("Verifying indel $number");
	system("perl verify_indel.pl $target $control $ref $number bi");
    }elsif(/$target.snp.[0-9][0-9]$/){
	$number = (split('\.', $_))[2];
	report("Verifying SNPs $number");
	system("perl verify_snp.pl $target $control $ref $number bi");
    }
}
closedir(DIR);

report("Making vcf file of SNP");
system("cat $target/$target.indel.verify.* > $target/$target.indel && rm $target/$target.indel.verify.* $target/$target.indel.??");
system("cat $target/$target.snp.verify.* > $target/$target.bi.snp && rm $target/$target.snp $target/$target.snp.verify.* $target/$target.snp.??");
system("perl snp2vcf.pl $target");
report("bidirectional.pl complete.");

sub report{
    my $message = shift;
    my $now = `date`;
    chomp($now);
    print "$now $message\n";
}

sub sortWait{
    my $file = shift;
    while(1){
	$mtime = (stat($file))[9];
	if (time > $mtime + 5){
	    return;
	}
	sleep 1;
    }
}
