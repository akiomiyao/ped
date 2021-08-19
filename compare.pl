#!/usr/bin/perl
#
# This file is a script for Polymorphic Edge Detection.
#
# Copyright (C) 2021 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#
# Author: MIYAO Akio <miyao@affrc.go.jp>
#

$usage = "
 $0 - output joined genotype data from two selected targets with primer sequence

 e.g. perl $0 a=SRR14477896,b=SRR14477897
      perl $0 a=SRR14477896,b=SRR14477897,type=sv
      Space after comma is not allowed.

      Alternatively,
      perl $0 SRR14477896 SRR14477897
      perl $0 SRR14477896 SRR14477897 sv

      Result will be output to the stanard output.

 Author: MIYAO Akio <miyao\@affrc.go.jp>
 Scripts: https://github.com/akiomiyao/ped
 Web page: https://akiomiyao.github.io/ped
 Docker: https://hub.docker.com/r/akiomiyao/ped  
";

if ($ARGV[0] eq ""){
    print $usage;
    exit;
}elsif($ARGV[0] =~ /,/){
    foreach (split(",", $ARGV[0])){
	($name, $value) = split("=", $_);
	$$name = $value;
    }
}else{
    $a = $ARGV[0];
    $b = $ARGV[1];
    $type = $ARGV[2];
}

if($type eq "sv"){
    open(IN, "$a/$a.sv.primer");
    print "Chr\tPosition\tChr\tPosition\tType\tSize\t$a\t$b\tLeft_primer\tRight_primer\tSize\n";
}else{
    open(IN, "$a/$a.bi.primer");
    print "Chr\tPosition\tRef\tAlt\t$a\t$b\tLeft_primer\tRight_primer\tSize\tSeqeunce aroud the mutation\n";
}
while(<IN>){
    chomp;
    @row = split;
    $row[1] = "00000000000$row[1]";
    $row[1] = substr($row[1], length($row[1])-11, 11);
    if ($type eq "sv"){
	$dat = "$row[0] $row[1] $row[2] $row[3] $row[5] $row[6]";
	$t{$dat} = "$row[13]\t$row[14]\t$row[15]";
	$a{$dat} = $row[11];
    }else{
	$dat = "$row[0] $row[1] $row[2] $row[3]";
	$t{$dat} = "$row[9]\t$row[10]\t$row[11]\t$row[12]";
	$a{$dat} = $row[8];
    }
}

if($type eq "sv"){
    open(IN, "$b/$b.sv.primer");
}else{
    open(IN, "$b/$b.bi.primer");
}
while(<IN>){
    chomp;
    @row = split;
    $row[1] = "00000000000$row[1]";
    $row[1] = substr($row[1], length($row[1])-11, 11);
    if ($type eq "sv"){
	$dat = "$row[0] $row[1] $row[2] $row[3] $row[5] $row[6]";
	$t{$dat} = "$row[13]\t$row[14]\t$row[15]";
	$b{$dat} = $row[11];
    }else{
	$dat = "$row[0] $row[1] $row[2] $row[3]";
	$t{$dat} = "$row[9]\t$row[10]\t$row[11]\t$row[12]";
	$b{$dat} = $row[8];
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
    $ga = $a{$_};
    $ga = "W" if $ga eq "";
    $gb = $b{$_};
    $gb = "W" if $gb eq "";
    print $ga;
    print "\t$gb";
    print "\t$t{$_}\n";
}
