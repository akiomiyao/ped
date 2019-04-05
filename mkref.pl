#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = "
     mkref.pl - script for making reference data. 

e.g. perl mkref.pl target
     perl mkref.pl TAIR10

For computer cluster,
     qsub -v target=target mkref.pl
     qsub -v target=hg38 mkref.pl

Author: Akio Miyao <miyao\@affrc.go.jp>

Currently, target listed below are available.

Customized reference can be made by adding or
modifiing the config file.

Name\tDescription
";

$cwd = `pwd`;
chomp($cwd);

$uname = `uname`;
chomp($uname);
if ($uname eq "FreeBSD"){
    $sort_opt = "-S 100M";
    $wget = "/usr/local/bin/wget";
}else{
    $wget = "/usr/bin/wget";
}

if ($ARGV[0] ne ""){
    $target = $ARGV[0];
}elsif ($ENV{target} ne ""){
    $target    = $ENV{target};
    $cwd       = $ENV{PBS_O_WORKDIR};
}

@nuc = ('A', 'C', 'G', 'T');

open(IN, "$cwd/config");
while(<IN>){
    chomp;
    @row = split('\t', $_);
    if ($row[1] eq "description"){
	$desc{$row[0]} = $row[2];
    }elsif($row[1] eq "wget"){
	$wget{$row[0]} = $row[2];
    }elsif($row[0] eq $target && $row[1] eq "chromosome"){
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

if ($target eq ""){
    print $usage;
    foreach (sort keys %desc){
	print "$_\t$desc{$_}\n"
    }
    exit;
}

if (! -d "$cwd/$target"){
    system("mkdir $cwd/$target");
}

chdir "$cwd/$target";

@row = split('/', $wget{$target});
$remote_file = $row[$#row];
if (! -e $remote_file){
    system("$wget -o wget-log $wget{$target}");
}

&mkChr;
&mk20;
&mkControlRead;

sub mkChr{
    my @file = split('/', $wget{$target});
    my $file = $file[$#file];
    my $i = 0;
    print "Making chromosome file\n";
    open(IN, "zcat $file|");
    while(<IN>){
	chomp;
	if (/^>/){
	    close(OUT);
	    $ref = "chr" . $chr[$i];
	    last if $chr[$i] eq "";
	    $i++;
	    open(OUT, "> $ref");
	}else{
	    y/a-z/A-Z/;
	    print OUT;
	}
    }
    close(IN);
    close(OUT);
}

sub mk20{
    print "Making 20mer position file.\n";
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		open($tag, "> ref20.$tag");
	    }
	}
    }
    
    foreach $i (@chr){
	print "Processing Chr$i\n";
	&mk20mer($i);
    }
    
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		close($tag);
	    }
	}
    }
    
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		&mkUniq($tag);
	    }
	}
    }    
}

sub mkControlRead{
    print "Making Control Read.\n";
    open(OUT, "|sort -T . $sort_opt |uniq > $target.sort_uniq");
    for(@chr){
	print "Processing Chr$_\n";
	open(IN, "chr$_");
	$data = <IN>;
	close(IN);
	$i = 0;
	while(1){
	    $read = substr($data, $i, 100);
	    last if length($read) != 100;
	    if ($read !~ /[MRWSYKVHDBN]/){
		print OUT "$read\n";
		print OUT &complement($read) . "\n";
	    }
	    $i += 2;
	}
    }
}

sub mkUniq{
    my $tag = shift;
    print "Making ref20_uniq.$tag.gz\n";
    system("sort -T . $sort_opt ref20.$tag > ref20_sort.$tag");
    open(OUT, "|gzip -f > ref20_uniq.$tag.gz");
    open(IN, "ref20_sort.$tag");
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
    sleep 1;
    system("rm ref20.$tag ref20_sort.$tag");
}

sub mk20mer{
    my $i;
    my $chr = shift;
    $file = "chr$chr";
    open(IN, $file);
    while(<IN>){
	$data = $_;
    }
    $length = length($data);
    for ($i = 0; $i <= $length - 20; $i++){
	$forward = substr($data, $i, 20);
	if ($forward !~ /[MRWSYKVHDBN]/){
	    $reverse = complement($forward);
	    $fpos = $i + 1;
	    $rpos = $fpos + 20 -1;
	    $tag = substr($forward, 0, 3);
	    print $tag "$forward\t$chr\t$fpos\tf\n";
	    $tag = substr($reverse, 0, 3);
	    print $tag "$reverse\t$chr\t$rpos\tr\n";
	}
    }
}

sub complement{
    my $seq = shift;
    my @seq = split('', $seq);
    my $length = length($seq);
    my $i = 0;
    my $out = "";
    for ($i = $length -1 ; $i >= 0; $i--){
        if ($seq[$i] eq "A"){
            $out .= "T";
        }elsif($seq[$i] eq "C"){
            $out .= "G";
        }elsif($seq[$i] eq "G"){
            $out .= "C";
        }elsif($seq[$i] eq "T"){
            $out .= "A";
        }
    }
    return $out;
}
