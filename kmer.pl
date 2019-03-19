#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     kmer.pl - pipeline of kmer method for standalone machine. 

     run  perl kmer.pl accession control reference
     e.g. perl kmer.pl SRR8181712 default TAIR10

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
     and then run perl kmer.pl mydata reference

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
    system("perl sort_uniq.pl $target");
}

if (! -e "$control/$control.sort_uniq"){
    report("Making $control.sort_uniq.");
    system("perl sort_uniq.pl $control");
}

if (! -e "$control/$control.lbc.AAA.gz"){
    report("Making kmer count of $control.");
    system("perl count.pl $control");
    system("gzip $control/$control.lbc.*");
}

report("Making kmer count of $target.");
system("perl count.pl $target");
report("Making snp");
system("perl snp.pl $target $control");
if ($ref eq ""){
     system("cat $target.snp.* >$target.kmer");
     exit;
}
report("Making map");
system("perl map.pl $target $ref");
report("Verifying");
system("perl verify_snp.pl $target $control $ref AAA kmer");

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
system("rm $target/$target.snp.* $target/$target.map.* $target/$target.kmer.verify.* $target/$target.kmer_chr.* $target/$target.lbc.* ");
system("perl snp2vcf.pl $target kmer");
report("kmer.pl done.");

sub report{
    my $message = shift;
    my $now = `date`;
    chomp($now);
    print "$now $message\n";
}
