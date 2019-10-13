#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2019 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

if ($cwd eq ""){
    $cwd = `pwd`;
    chomp($cwd);
}

$uname = `uname`;
chomp($uname);
if ($uname eq "FreeBSD"){
    $wget = "/usr/local/bin/wget";
    $rsync = "/usr/local/bin/rsync";
    $sort_opt = "-S 100M";
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
		    $qsub = "-v target=$target,tag=$tag $cwd/sort_uniq_sub.pl";
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
    open(IN, "qsub $qsub |");
    $job = <IN>;
    close(IN);
    chomp($job);
    push(@job, $job);
    sleep 1;
}

sub checkQsub{
    my ($job, @dat, $num);
    chdir $cwd;
    opendir(DIR, '.');
    foreach (sort readdir(DIR)){
	if(/.e[0-9]+$/){
	    @dat = split('\.', $_);
	    ($num = $dat[$#dat]) =~ y/0-9//cd;
	    if (-s $_ > 0){
		print "Error JobId $num.s2\n";
		system("cat $_");
	    }else{
		system("rm *.e$num *.o$num");
	    }
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
		    if (/^[qre]$/i){
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
