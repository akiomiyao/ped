#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2019 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

# $qsub_opt = "-S /usr/bin/perl -jc M.c8"; # example of additional option.
$sort_opt = "-S 1M";

$sort_opt = "-S 1M"; # This value is very small. But if the buffer size is specified to 100M or more, Gnu sort sometimes freezes.

if ($cwd eq ""){
    $cwd = `pwd`;
    chomp($cwd);
}

$uname = `uname`;
chomp($uname);
if ($uname eq "FreeBSD"){
    $wget = "/usr/local/bin/wget";
    $rsync = "/usr/local/bin/rsync";
}else{
    $wget = "/usr/bin/wget";
    $rsync = "/usr/bin/rsync";
}

@nuc = ('A', 'C', 'G', 'T');

sub mkSortUniq{
    my $target = shift;
    my ($job, @row, $flag, $count_gz, $count_seq);
    if (! -e "$target/sort_uniq"){
	system ("mkdir $target/sort_uniq");
    }
    opendir(DIR, "$target/sort_uniq");
    foreach (readdir(DIR)){
	if (/gz/){
	    $count_gz++;
	}elsif(/seq/){
	    $count_seq++;
	}
	
    }
    if ($count_gz == 64){
	return;
    }elsif($count_seq == 64){
	foreach $nuca (@nuc){
	    foreach $nucb (@nuc){
		foreach $nucc (@nuc){
		    $tag = $nuca . $nucb . $nucc;
		    if ($tmpdir eq ""){
			$qsub = "-v target=$target,tag=$tag sort_uniq_sub.pl";
		    }else{
			$qsub = "-v target=$target,tag=$tag,tmpdir=$tmpdir sort_uniq_sub.pl";
		    }
		    &doQsub($qsub);
		}
	    }
	}
    }else{
	&report("Making $target.sort_uniq.");
	$qsub = "-v target=$target qsub_sort_uniq.pl";
	&doQsub($qsub);
    }
}

sub report{
    my $message = shift;
    my $now = `date`;
    chomp($now);
    print "$now $message\n";
}

sub doQsub{
    my $qsub = shift;
    my $job;
    open(IN, "qsub $qsub_opt $qsub |");
    $job = <IN>;
    close(IN);
    chomp($job);
    $job = (split(/\s+/, $job))[2] if $job =~ /Your job/;
    push(@job, $job);
    sleep 1;
}

sub cleanupLog{
    my ($job, @dat, $num);
    chdir $cwd;
    opendir(DIR, '.');
    foreach (sort readdir(DIR)){
	if(/\.pl\.e[0-9]+$/){
	    if (-s $_ == 0){
		system("rm $_");
	    }
	}elsif(/\.pl\.o[0-9]+$/){
	    system("rm $_");
	}
    }
}

sub holdUntilJobEnd{
    my (@row,  %stat, $job, $flag);
    while(1){
	%stat = ();
	$flag = 1;
	open(IN, "qstat |");
	while(<IN>){
	    chomp;
	    @row = split;
	    $stat{$row[0]} = $_;
	}
	close(IN);
	foreach $job(@job){
	    if($stat{$job} ne ""){
		@row = split(/\s+/, $stat{$job});
		foreach (@row){
		    if (/^[qre]$|^qw$/i){
			$flag = 0;
		    }
		}
	    }
	}
	if ($flag){
	    last;
	}
	sleep 10;
    }
}

sub waitFile{
    my $file = shift;
    while(1){
	last if -e $file;
	sleep 10;
    }
    while(1){
	$mtime = (stat($file))[9];
	if (time > $mtime + 5){
	    return;
	}
	sleep 1;
    }
}

sub closeTag{
    foreach $nuca (@nuc){
	foreach $nucb (@nuc){
	    foreach $nucc (@nuc){
		$tag = $nuca . $nucb . $nucc;
		close($tag);
	    }
	}
    }
}

sub complement{
    my $seq = shift;
    my @seq = split('', $seq);
    my $i = length($seq);
    my $out = "";
    while($i > 0){
        $i--;
        if ($seq[$i] eq "A"){
            $out .= "T";
        }elsif($seq[$i] eq "C"){
            $out .= "G";
        }elsif($seq[$i] eq "G"){
            $out .= "C";
        }elsif($seq[$i] eq "T"){
            $out .= "A";
        }else{
            $out .= "N";
        }
    }
    return $out;
}

sub bynumber{
    $a <=> $b;
}

1;
