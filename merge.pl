#!/usr/bin/perl
#
# This file is a part of PED; a script for Polymorphic Edge Detection.
#
# Copyright (C) 2017 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/ped
#

$usage = '
     merge.pl - merge count.pl data. 

e.g. qsub -v target=ERR194147,tag=AAA,tmpdir=/mnt/ssd merge.pl

     target is name of target.
     tmpdir can be ommited.

Author: Akio Miyao <miyao@affrc.go.jp>

';

if($ENV{target} ne ""){
    $target    = $ENV{target};
    $tag       = $ENV{tag};
    $cwd       = $ENV{PBS_O_WORKDIR};
    $tmpdir    = $ENV{tmpdir};
    $workdir = "$cwd/$target";
}else{
    print $usage;
    exit;
}

if ($tmpdir eq ""){
    $tmpdir = $workdir;
}else{
    system("mkdir $tmpdir/$target");
    $tmpdir = "$tmpdir/$target";
}

chdir $tmpdir;

system("cp $workdir/*.$tag.*.gz . && gzip -d *.gz");

while(1){
    opendir(DIR, ".");
    @file = sort(grep(/$tag/, readdir(DIR)));
    closedir(DIR);
    $total = @file;
    if ($total == 1){
	last;
    }
    @last = split('\.', $file[$#file]);
    $last = $last[$#last];
    $last ++;
    $last = "000" . $last;
    $last = substr($last, length($last) - 4, 4);
    $output = "$last[0].$last[1].$last[2].$last";
    $cmd = "join -a 1 -a 2 $file[0] $file[1] | awk '{print \$1 \"\t\" \$2 + \$3}' > $output && rm $file[0] $file[1]";
    system($cmd);
}

$final = "$last[0].$last[1].$last[2]";
system("mv $output $final && gzip $final && mv $final.gz $workdir && cd .. && rm -r $target && rm $workdir/*.$tag.*.gz");
