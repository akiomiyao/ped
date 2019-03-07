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

     run  perl bidirectional.pl accession reference
     e.g. perl bidirectional.pl SRR8181712 TAIR10

     This script is for stand alone machine.

     Before run the script, downloading target reads and reference genome data
     are required.

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

$target = $ARGV[0];
$ref    = $ARGV[1];

if (! -e "$target/$target.sort_uniq"){
    report("Making $target.sort_uniq.");
    system("perl sort_uniq.pl $target");
}

report("Aligning of $target sequence to $ref genome.");
system("perl align.pl $target $ref 5");

## If you want more sensitivity, remove #s below.
#system("perl align.pl $target $ref 0");
#system("perl align.pl $target $ref 10");
#system("perl align.pl $target $ref 15");

report("Splitting alignment of Indel.");
system("perl split_indel.pl $target $ref");
report("Splitting alignment of SNP.");
system("perl split_snp.pl $target $ref");

opendir(DIR, $target);
foreach(sort readdir(DIR)){
    if (/$target.indel.[0-9][0-9]$/){
	$number = (split('\.', $_))[2];
	report("Verifying indel $number");
	system("perl verify_indel.pl $target default $ref $number bi");
    }elsif(/$target.snp.[0-9][0-9]$/){
	$number = (split('\.', $_))[2];
	report("Verifying SNPs $number");
	system("perl verify_snp.pl $target default $ref $number bi");
    }
}
closedir(DIR);

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
