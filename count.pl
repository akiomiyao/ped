#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     count.pl - create count of kmers from sort_uniq data. 

e.g. qsub -v target=ERR194147,number=0001,tmpdir=/mnt/ssd count.pl

     For stand alone machine
e.g. perl count.pl SRR8181712
     perl count.pl DRR054198

     target is name of target.
     number is number of split target file (e.g. 0001).
     tmpdir can be ommited.

Author: Akio Miyao <miyao@affrc.go.jp>

';

$uname = `uname`;
chomp($uname);
if ($uname eq "FreeBSD"){
    $sort_opt = "-S 100M";
}

if($ENV{target} ne ""){
    $target    = $ENV{target};
    $number    = $ENV{number};
    $cwd       = $ENV{PBS_O_WORKDIR};
    $tmpdir    = $ENV{tmpdir};
    $workdir = "$cwd/$target";
    &cluster;
}elsif($ARGV[0] ne ""){
    &standalone;
}else{
    print $usage;
}

sub cluster{
    if ($tmpdir eq ""){
	$tmpdir = $workdir;
	$input_file = "$workdir/tmp/$target.sort_uniq.$number";
    }else{
	if (-e "$tmpdir/$target"){
	    system("rm -r $tmpdir/$target");
	}
	system("mkdir $tmpdir/$target");
	$tmpdir = "$tmpdir/$target";
	system("cp $workdir/tmp/$target.sort_uniq.$number $tmpdir");
	$input_file = "$target.sort_uniq.$number";
    }

    chdir $tmpdir;
    
    open(IN, "$input_file");
    open(OUT, "|sort $sort_opt -T $tmpdir | uniq -c | awk '{print \$2 \"\t\" \$1}' | /usr/bin/perl $cwd/split_count.pl");
    while(<IN>){
	chomp;
	$length = length($_);
	for ($i = 0; $i <= $length - 20; $i++){
	    $seq = substr($_, $i, 20);
	    if ($seq !~ /N/){
		print OUT "$seq\n";
	    }
	}
    }
    close(IN);
    close(OUT);
}

sub standalone{
    $target = $ARGV[0];
    chdir $target;
    $input_file = $target . ".sort_uniq";
    if (-e "tmp"){
	system("rm -rf tmp");
    }
    system("mkdir tmp");
    $count = 0;
    $file_count = 1;
    open(IN, $input_file);
    open(OUT, "> tmp/$input_file.$file_count");
    while(<IN>){
	if ($count == 10000000){
	    $file_count ++;
	    open(OUT, "> tmp/$target.sort_uniq.$file_count");
	    $count = 0;
	}
	$count++;
	print OUT;
    }
    close(IN);
    close(OUT);
    for ($j = 1; $j <= $file_count; $j++){
	open(IN, "tmp/$input_file.$j");
	$number = "000$j";
	$number = substr($number, length($number) -4, 4);
	&report("$input_file.$j / $file_count is processing.");
	$reads = 0;
	open(OUT, "|sort $sort_opt -T . | uniq -c | awk '{print \$2 \"\t\" \$1}' > $target.count.$number");
	while(<IN>){
	    $reads++;
	    chomp;
	    $length = length($_);
	    for ($i = 0; $i <= $length - 20; $i++){
		$seq = substr($_, $i, 20);
		if ($seq !~ /N/){
		    print OUT "$seq\n";
		}
	    }
	}
	close(IN);
	close(OUT);
	&split_count($number);
    }
    &merge;
}

sub split_count{
    my $number = shift;
    open(IN, "$target.count.$number");
    while(<IN>){
	$head = substr($_, 0, 3);
	$head{$head} = 1;
	if ($prev ne $head){
	    $output_file = "$target.count.$head.$number";
	    open(OUT, "> $output_file");
	}
	print OUT "$_";
	$prev = $head;
    }
    close(IN);
    close(OUT);
    system("rm $target.count.$number");
}

sub merge{
    foreach $head (sort keys %head){
	while(1){
	    @file = ();
	    opendir(DIR, ".");
	    @file = sort(grep(/count\.$head\.[0-9]/, readdir(DIR)));
	    closedir(DIR);
	    $total = @file;
	    if ($total == 1){
		$tag = "";
		$prev = "";
		$final = "$target.lbc.$head";
		open(IN, $output);
		open(OUT, "> $final");
		while(<IN>){
		    chomp;
		    ($seq, $count) = split;
		    $tag = substr($seq, 0, 19);
		    $nuc = substr($seq, 19, 1);
		    
		    if ($tag ne $prev and $prev ne ""){
			print OUT "$prev\t$A\t$C\t$G\t$T\n";
			$A = 0;
			$C = 0;
			$G = 0;
			$T = 0;
		    }
		    $$nuc = $count;
		    $prev = $tag;
		}
		close(IN);
		print OUT "$prev\t$A\t$C\t$G\t$T\n";
		close(OUT);
		system("rm $target.count.$head.*");
		last;
	    }
	    @last = split('\.', $file[$#file]);
	    $last = $last[$#last];
	    $last ++;
	    $last = "000" . $last;
	    $last = substr($last, length($last) - 4, 4);
	    $last[$#last] = "";
	    
	    $output = join('.', @last) . $last;
	    $cmd = "join -a 1 -a 2 $file[0] $file[1] | awk '{print \$1 \"\t\" \$2 + \$3}' > $output";
	    system($cmd);
	    system("rm $file[0] $file[1]");
	}
    }
}

sub report{
    my $message = shift;
    my $now = `date`;
    chomp($now);
    print "$now $message\n";
}
