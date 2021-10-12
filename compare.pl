#!/usr/bin/perl
#
# This file is a script for Polymorphic Edge Detection.
#
# Copyright (C) 2021 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#
# Author: Akio Miyao <miyao@affrc.go.jp>
#

$usage = "
 $0 - output joined genotype data from two selected targets with primer sequence

 e.g. perl $0 target=SRR14477896:SRR14477897:SRR14477898
      perl $0 target=SRR14477896:SRR14477897:SRR14477898,type=sv
      Space after comma is not allowed.

      Result will be output to the standard output.
      Genotypes following primer data will be appear with specified order of
      target separated by ':'.
      Number of targets can be changed.

 Author: Akio Miyao <miyao\@affrc.go.jp>
 Scripts: https://github.com/akiomiyao/ped
 Web page: https://akiomiyao.github.io/ped
 Docker: https://hub.docker.com/r/akiomiyao/ped  

";

if ($ARGV[0] eq ""){
    print $usage;
    exit;
}elsif($ARGV[0] ne ""){
    foreach (split(",", $ARGV[0])){
	($name, $value) = split("=", $_);
	$$name = $value;
    }
}

@target = split(':', $target);

if($type eq "sv"){
    print "Chr\tPosition\tChr\tPosition\tType\tSize\tLeft_primer\tRight_primer\tSize";
}else{
    print "Chr\tPosition\tRef\tAlt\tLeft_primer\tRight_primer\tSize\tSequence around the mutation";
}


foreach $name (@target){
    print "\t$name";
}
print "\n";

foreach $name (@target){
    if($type eq "sv"){
	open(IN, "$name/$name.sv.primer");
    }else{
	open(IN, "$name/$name.bi.primer");
    }
    while(<IN>){
	chomp;
	@row = split;
	$row[1] = "00000000000$row[1]";
	$row[1] = substr($row[1], length($row[1])-11, 11);
	if ($type eq "sv"){
	    $dat = "$row[0] $row[1] $row[2] $row[3] $row[5] $row[6]";
	    $t{$dat} = "$row[13]\t$row[14]\t$row[15]";
	    $$name{$dat} = $row[11];
	}else{
	    $dat = "$row[0] $row[1] $row[2] $row[3]";
	    $t{$dat} = "$row[9]\t$row[10]\t$row[11]\t$row[12]";
	    $$name{$dat} = $row[8];
	}
    }
}

foreach (sort keys %t){
    @row = split;
    $row[1] += 0;
    if ($type eq "sv"){
	print "$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t";
    }else{
	print "$row[0]\t$row[1]\t$row[2]\t$row[3]\t";
    }
    print "$t{$_}";
    foreach $name (@target){
	$gt = $$name{$_};
	$gt = "W" if $gt eq "";
	print "\t$gt";
    }
    print "\n";
}
