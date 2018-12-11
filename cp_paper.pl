
$target = $ARGV[0];

system("cat $target/$target.indel.verify.* > ~/paper/trio/$target.indel");
system("cp $target/$target.indel.sort  ~/paper/trio/");
system("cat $target/$target.snp.verify.* > ~/paper/trio/$target.snp");
