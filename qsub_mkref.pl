#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

if ($ENV{PBS_O_WORKDIR} ne ""){
    $cwd = $ENV{PBS_O_WORKDIR};
    chdir $cwd;
    require "$cwd/common.pl";
}elsif($ENV{SGE_O_WORKDIR} ne ""){
    $cwd = $ENV{SGE_O_WORKDIR};
    chdir $cwd;
    require "$cwd/common.pl";
}else{
    require './common.pl';
}

$usage = "
     qsub_mkref.pl - script for making reference data using computer cluster. 

     perl qsub_mkref.pl target
     e.g. perl qsub_mkref.pl TAIR10
     The qsub_mkref.pl runs on master server. Subprocess will run on nodes.

     If you want to run qsub_mkref.pl by a node,
     qsub -v target=target qsub_mkref.pl
     e.g. qsub -v target=TAIR10 qsub_mkref.pl

     Currently, target listed below are available.

     Target Name    Description
";

if ($ARGV[0] ne ""){
    $target = $ARGV[0];
    $file   = $ARGV[1];
}elsif ($ENV{target} ne ""){
    $target    = $ENV{target};
    $file      = $ENV{file};
    $tag       = $ENV{tag};
    $chr       = $ENV{chr};
    $execution = $ENV{execution};
}

&report("Job executed: qsub_mkref.pl");

open(IN, "$cwd/config");
while(<IN>){
    chomp;
    @row = split('\t', $_);
    if ($row[1] eq "description"){
	$desc{$row[0]} = $row[2];
    }elsif($row[1] eq "wget"){
	$wget{$row[0]} = $row[2];
    }elsif($row[1] eq "file"){
	$file{$row[0]} = $row[2];
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
	$name = $_ . "            ";
	$name = substr($name, 0, 14);
	print "      " . $name . "$desc{$_}\n";
    }
    print "
     Customized reference can be made by adding or modifiing the config file.

     If you want to make reference data which is not listed in config file,
     perl qsub_mkref.pl taget file_name

     Before run the script,
     mkdir ./target
     cp fasta_file ./target

     For example,
     mkdir ./RefBeet12
     cp somewhere/RefBeet-1.2.fna ./RefBeet12
     perl qsub_mkref.pl RefBeet12 RefBeet-1.2.fna

     Author: Akio Miyao <miyao\@affrc.go.jp>

";
    exit;
}

if ($execution eq ""){
    if (! -d "$cwd/$target"){
	system("mkdir $cwd/$target");
    }
    
    chdir "$cwd/$target";
    
    if ($file eq ""){
	@row = split('/', $wget{$target});
	$remote_file = $row[$#row];
	if (! -e $remote_file and $wget{$target} ne ""){
	    &report("Downloading $wget{$target}");
	    system("$wget -o wget-log $wget{$target}");
	}
	&mkChr;
    }else{
	&mkChrFromFile($file);
    }

    chdir $cwd;

    &report("Making 20mer position file for chromosome");
    foreach $chr (@chr){
	$qsub = "-v target=$target,chr=$chr,execution=mk20 $cwd/qsub_mkref.pl";
	&doQsub($qsub);
	$chr_count ++;
    }
    while(1){
	$file_count = 0;
	opendir(DIR, "$cwd/$target");
	foreach (readdir(DIR)){
	    if (/^ref20./){
		$mtime = (stat("$cwd/$target/$_"))[9];
		$ti = time - $mtime;
		if (time > $mtime + 10){
		    $file_count++;
		}
	    }
	}
	closedir(DIR);
	if ($file_count == $chr_count * 64){
	    last;
	}
	sleep 10;
    }
 
    &report("Selecting uniqe position");
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		$qsub = "-v target=$target,tag=$tag,execution=mkUniq $cwd/qsub_mkref.pl";
		&doQsub($qsub);

	    }
	}
    }
 
    &report("Making control read");
    foreach $chr (@chr){
	$qsub = "-v target=$target,tag=$tag,chr=$chr,execution=mkControlRead $cwd/qsub_mkref.pl";
	&doQsub($qsub);
    }
    while(1){
	$file_count = 0;
	opendir(DIR, "$cwd/$target/sort_uniq");
	foreach (readdir(DIR)){
	    if (/^[ACGT][ACGT][ACGT].tmp/){
		$mtime = (stat("$cwd/$target/sort_uniq/$_"))[9];
		if (time > $mtime + 10){
		    $file_count++;
		}
	    }
	}
	closedir(DIR);
	if ($file_count == $chr_count * 64){
	    last;
	}
	sleep 10;
    }

    &report("Sorting control read");
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		$qsub = "-v target=$target,tag=$tag,execution=sortControlRead $cwd/qsub_mkref.pl";
		&doQsub($qsub);

	    }
	}
    }

    while(1){
	$file_count = 0;
	opendir(DIR, "$cwd/$target/sort_uniq");
	foreach (readdir(DIR)){
	    if (/[ACGT][ACGT][ACGT].gz$/){
		$mtime = (stat("$cwd/$target/sort_uniq/$_"))[9];
		if (time > $mtime + 10){
		    $file_count++;
		}
	    }
	}
	closedir(DIR);
	if ($file_count ==  64){
	    last;
	}
	sleep 10;
    }

    &report("Job complete: qsub_mkref.pl");
}elsif($execution eq "mk20"){
    &mk20;
}elsif($execution eq "mkUniq"){
    &mkUniq;
}elsif($execution eq "mkControlRead"){
    &mkControlRead;
}elsif($execution eq "sortControlRead"){
    &sortControlRead;
}

sub mkChr{
    chdir "$cwd/$target";
    my @file = split('/', $wget{$target});
    my $file = $file[$#file];
    $file =~ s/\ +$//g;
    my $i = 0;
    if ($target eq "hg38"){
	open(IN, "$zcat $file|");
	while(<IN>){
	    chomp;
	    if (/^>/){
		close(OUT);
		$ref = (split)[0];
		$ref =~ s/>//;
		if ($ref =~ /\_/){
		    $flag = 1;
		}else{
		    $flag = 0;
		    $ref =~ s/chr//i;
		    $ref += 0 if $ref =~ /^[0-9]*$/;
		    $ref = "chr$ref";
		    open(OUT, "> $ref");
		}
	    }elsif(! $flag){
		y/a-z/A-Z/;
		print OUT;
	    }
	}
	close(IN);
	close(OUT);
    }else{
	if ($target eq "IWGSC1.0"){
	    system("unzip iwgsc_refseqv1.0_all_chromosomes.zip && mv iwgsc_refseqv1.0_all_chromosomes/$file{$target} .");
	}elsif($target eq "SL3"){
	    system("tar xvfz $file && rm $file");
	    $file =~ s/\.tar\.gz//;
 	}elsif($target =~ /^B73/){
	    for ($i = 1; $i <= 10; $i++){
		$file = $wget{$target} . "chr$i.fna.gz";
		&report("Downloading $file");
		system("$wget -o wget-log $file");
		$file = "chr$i.fna.gz";
		open(IN, "$zcat $file|");
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
	    open(IN, "$zcat $file|");
	}elsif ($file =~ /bz2$/){
	    open(IN, "bzcat $file|");
	}elsif ($file =~ /xz$/){
	    open(IN, "xzcat $file|");
	}else{
	    $file = $file{$target} if $file{$target} ne "";
	    open(IN, $file);
	}
	while(<IN>){
	    chomp;
	    if (/^>/){
		close(OUT);
		$ref = "chr" . $chr[$i];
		last if $chr[$i] eq "";
		open(OUT, "> $ref");
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

sub mkChrFromFile{
    chdir "$cwd/$target";
    my $file = shift;
    if ($file eq ""){
	my @file = split('/', $wget{$target});
	$file = $file[$#file];
    }
    if ($file =~ /gz$/){
	open(IN, "$zcat $file|");
    }elsif($file =~ /bz2$/){
	open(IN, "bzcat $file|");
    }else{
	open(IN, $file);
    }
    while(<IN>){
	chomp;
	if (/^>/){
	    close(OUT);
	    $ref = (split)[0];
	    $ref =~ s/>//;
	    $ref =~ s/^chr//i;
	    $ref += 0 if $ref =~ /^[0-9]*$/;
	    push(@chr, $ref) if ! $chr_flag;
	    $ref = "chr$ref";
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
    chdir "$cwd/$target";
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		open($tag, "> ref20.$tag.$chr");
	    }
	}
    }
    
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
    &closeTag;
}

sub mkUniq{
    system("cat $cwd/$target/ref20.$tag.* |sort -T $cwd/$target $sort_opt  > $cwd/$target/ref20_sort.$tag");
    &waitFile("$cwd/$target/ref20_sort.$tag");
    open(OUT, "|gzip -f > $cwd/$target/ref20_uniq.$tag.gz");
    open(IN, "$cwd/$target/ref20_sort.$tag");
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
    &waitFile("$cwd/$target/ref20_uniq.$tag.gz");
    system("rm $cwd/$target/ref20.$tag.* $cwd/$target/ref20_sort.$tag");
 }

sub mkControlRead{
    if (! -e "$target/sort_uniq"){
	system("mkdir $target/sort_uniq");
    }
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		open($tag, "> $cwd/$target/sort_uniq/$tag.tmp.$chr")
	    }
	}
    }
    open(IN, "$cwd/$target/chr$chr");
    $data = <IN>;
    close(IN);
    $i = 0;
    while(1){
	$read = substr($data, $i, 100);
	last if length($read) != 100;
	if ($read !~ /[MRWSYKVHDBN]/){
	    $tag = substr($read, 0, 3);
	    print $tag "$read\n";
	    $read = &complement($read);
	    $tag = substr($read, 0, 3);
	    print $tag "$read\n";
	}
	$i += 2;
    }
    &closeTag;
}

sub sortControlRead{
    system("cat $cwd/$target/sort_uniq/$tag.tmp.* | sort -T $cwd/$target $sort_opt |uniq |gzip > $cwd/$target/sort_uniq/$target.$tag.gz");
    &waitFile("$cwd/$target/sort_uniq/$target.$tag.gz");
    system("rm $cwd/$target/sort_uniq/$tag.tmp.*");
}