# PED : Polymorphic Edge Detection

Polymorphic Edge Detection (PED) is the analysis flow for DNA polymorphism detection from short reads of next generation sequencer (NGS). We developed two methods to detect polymorphisms based on detection of the polymorphic edge. One is based on bidirectional alignment and the other is based on comparison of *k*-mers.

### Polymorphic Edge

DNA polymorphism is any difference of DNA sequence between individuals. These differences are single nucleotide polymorphism (SNP), insertion, deletion, inversion, translocation and copy number variation. On the non-polymorphic region, sequences between two individuals are completely same. At the position of SNP, or at the beginning of other polymorphisms, the nucleotide must be different between individuals.

### Bidirectional alignment method

                                                                    Chr11 80443004
                                                                    |
    TTTTTAATTGAAAAGGCATTAAGCTGGGTCTATGCAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTAGATAGGTAGAAAAAAAAAACCACTATCAGCAACA Reference sequence matching from 5'-end
    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  | |     | ||||||||   |  |       |  
    TTTTTAATTGAAAAGGCATTAAGCTGGGTCTATGCAGTGTGTGTGTGTGTGTGTGTGTGTGTGTAGATAGGTAGAAAAAAAAAACCACTATCAGCAACAGT Short read sequence
    |||       ||             |  | |     |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    TTTAATTGAAAAGGCATTAAGCTGGGTCTATGCAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTAGATAGGTAGAAAAAAAAAACCACTATCAGCAACAGT Reference sequence matching from 3'-end
                                       |
                                       Chr11 80442977

Short read sequence is aligned with reference sequence from both 5'- and 3'-ends. Positions indicated over and bellow of the alignment are first mismatched base, *i.e.,* polymorphic edge. The bidirectional alignment clearly indicates two bases (GT) deletion in the short read. The bidirectional can detect not only deletion but also SNP, insertion, inversion and translocation.

### *K*-mer method

    Individual_A AAATGGTACATTTATATTAT
    Individual_B AAATGGTACATTTATATTAC
          
All short reads from Individual_A and Individual_B are sliced to *k*-mer (*e.g. k* = 20) in each position. For example, the Individual_A has the *k*-mer sequence of AAATGGTACATTTATATTAT but does not have AAATGGTACATTTATATTAC. On the other hand, the Individual_B has the AAATGGTACATTTATATTAC but does not have AAATGGTACATTTATATTAT. The last base of *k*-mer of Individual_A is T, and Individual_B is C. The last base of *k*-mers must be SNP or edge of insertion, deletion, inversion, translocation or copy number variation. The *k*-mer method detects edges of polymorphism by difference of last base of *k*-mers. This method enables to detect polymorphisms by direct comparison of NGS data.



## Installation

- Programs run on Unix platforms (Linux, MacOS, FreeBSD). 
- Download zip file of PED from https://github.com/akiomiyao/ped and extract.  
or  
% git clone https://github.com/akiomiyao/ped.git  

- To download sequence data, fastq-dump from NCBI is required.  
    Tool kit can be download from  
    https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/ 
- To download reference date, wget is required.  
    If your machine do not have wget program, install wget from package.
    
## Simple instruction for bidirectional method  
- % perl mkref.pl reference  
  For example,  
  % perl mkref.pl WBcel235  
  Directory WBcel235 for reference of Caenorhabditis elegans WBcel235 will be created.  
- % perl download.pl accession  
  For example,  
  % perl download.pl ERR3063486  
  % perl download.pl ERR3063487   
  Data directory of ERR3063486 and ERR3063487 will be created.  
  Fastq data will be downloaded in read subdirectory.  
  ERR3063486 is read data of Caenorhabditis elegans mutant.  
  ERR3063487 is read data of Caenorhabditis elegans wild-type.  
- If you want to analyze local file,  
  % mkdir mydata1  
  % mkdir mydata1/read  
  % cp somewhere/mydata1.fastq mydata1/read  
  % perl bidirectional.pl mydata1 control reference  
- % perl bidirectional.pl target control reference  
  For example,  
  % perl bidirectional.pl ERR3063487 default WBcel235  
   After two hours, you will find results in ERR3063487 directory.  
   ERR3063487.indel is list of structural variation.  
   ERR3063487.snp is list of SNPs.  
   ERR3063487.snp.vcf is the vcf file for SNPs.  
    
- Our verify process counts reads containing polymorphic region.  
  Basically, counts are from target and reference.  
  If control is specified, counts from target and control will be listed.  
  For example,  
  % perl bidirectional.pl ERR3063487 ERR3063486 WBcel235  
  retruns 'M' or 'H' marked SNPs and indels of ERR3063487 which are absent in ERR3063486.  
- To confirm the alignment for detected polymorphisms,  
  % perl search.pl target chr position  
  e.g. % perl search.pl ERR3063487 II 948033  
  Alignments will be selected by the search script.  
- if you want to run with computer cluster,  
  % perl qsub_bidirectional.pl ERR3063487 default WBcel235  
- Run without arguments, help for script will be shown.  
- I searched small size short reads suitable for demonstration of PED from SRA in NCBI.  
  I found data set of Caenorhabditis elegans.
  But, I could not contact to uploaded scientist because of no email or contact address was descirbed.  
  I express special thanks for the release of short reads to SRA by ENS Lyon.  

## Simple instruction for kmer method
- % perl kmer.pl target control reference  
- % perl kmer.pl ERR3063487 ERR3063486 WBcel235  
- if you want to run with computer cluster,  
  % perl qsub_kmer.pl ERR3063487 ERR3063486 WBcel235
  ERR3063487 specific SNPs will be detected.  
- if you want to detect SNPs against reference genome,  
  % perl kmer.pl ERR3063487 default WBcel235

## Making reference data sets
- % perl mkref.pl  
  Without specify reference, reference name and description are listed.  
  If you want to make new reference, add the information about the reference to 'config' file.  
  For example,  
  % perl mkref.pl hg38  
    Data set of human genome hg38 will be made. It takes about two days.  
  % perl mkref.pl GRCm38  
    Data set of mouse genome GRCm38 will be made. It takes about two days.  
  % perl mkref.pl dmel626  
    Data set of Drosophila melanogaster r6.26 will be made. It takes one hour.  
  % perl mkref.pl IRGSP1.0  
    Data set of rice (Olyza sativa L. cv. Nipponbare) will be made.  
  % perl mkref.pl TAIR10  
    Data set of Arabidopsis thaliana will be made.  
- To run by computer cluster,  
  % qsub -v target=target mkref.pl  
  For example,    
  % qsub -v target=hg38 mkref.pl  

## Examples of result
A part of SNP list of bidirectional method is  
```
1       994949  C       T       11      50      0       21      17      H
1       994962  G       T       1       50      0       21      0       
1       994997  C       T       11      50      0       21      23      H
1       995371  C       G       5       50      0       15      24      H
1       995450  A       G       1       50      0       22      0       
1       995512  T       C       6       50      0       0       22      M
1       995543  A       G       6       50      0       0       22      M

Column 1: Chromosome number
Column 2: Position of SNP
Column 3: Reference base at the SNP position
Column 4: Alternative base
Column 5: Number of detected reads with alternative base
Column 6: Number of reads in the control sort_uniq file with control type polymorphism
Column 7: Number of reads in the control sort_uniq file with target type polymorphism
Column 8: Number of reads in the target sort_uniq file with control type polymorphism
Column 9: Number of reads in the target sort_uniq file with target type polymorphism
Column 10: Genotype (M: homozygous, H: heterozygous)
```

- A part of indel result is
```
1       923312  1       923312  f       deletion        -1      50      0       0       31      M
1       931147  1       931131  f       insertion       4       50      0       5       19              CCCTCCCTCCC
1       932618  1       932617  f       deletion        -4      50      0       0       25      M       TTTC
1       933555  1       933548  f       insertion       1       50      0       0       0               GGGGG
1       933750  1       933741  f       insertion       1       50      0       4       20              GGGGGGG
1       937406  1       937398  f       deletion        -2      50      0       6       11      H       CCCCCCCCC
1       939490  1       939445  f       insertion       36      50      0       32      1               ATCTCCCC
1       939575  1       939570  f       insertion       12      50      0       23      15      H       GGAGGACC

Column 1: Chromosome number of junction detected by 5' to 3' matching
Column 2: Position of junction detected by 5' to 3' matching
Column 3: Chromosome number of junction detected by 3' to 5' matching
Column 4: Position of junction detected by 3' to 5' matching
Column 5: Direction
Column 6: Type of polymorphism (insertion, deletion, inversion and translocation)
Column 7: Size of insertion or deletion
Column 8: Number of reads in the control sort_uniq file with control type polymorphism
Column 9: Number of reads in the control sort_uniq file with target type polymorphism
Column 10: Number of reads in the target sort_uniq file with control type polymorphism
Column 11: Number of reads in the target sort_uniq file with target type polymorphism
Column 12: Genotype (M: homozygous, H: heterozygous)
Column 13: Sequence between junctions
```

- A part of SNP result by kmer method is
```
X       54009891        AAAAAAAAAAGTGGCTCTT     T       T       GT      f       0       0       0       40      1       1       17      21      50      0       12      0       
6       112904084       AAACGACACTTTTTTTTTT     C       C       AC      r       0       0       40      0       0       1       31      22      50      0       30      0       
3       125369417       AAAAAAAAAAGTGTGCCTC     T       T       AT      f       0       0       0       40      19      1       0       17      50      0       12      0       
X       154687353       AAAAAAAAAAGTGTTAGGC     C       C       CT      f       0       40      0       0       1       21      1       21      50      0       1       0       
3       139458183       CGCCAACACTTTTTTTTTT     T       T       CT      r       41      0       0       0       46      1       12      0       50      0       20      0       
6       147492241       AAAAAAAAAAGTTATGATC     G       G       C       f       0       0       40      0       7       28      0       1       50      0       0       34      M
8       10651787        AAAAAAAAAAGTTATGCTC     A       A       AT      f       40      0       0       0       35      0       0       17      50      0       12      0       
1       155410117       AAAAAAAAAAGTTCGTATT     A       A       AT      f       40      0       0       0       44      1       0       10      50      0       33      0       
14      100919244       TAAACGAACTTTTTTTTTT     A       A       AT      r       0       0       0       40      11      1       0       44      50      0       12      0       
5       59826499        ATCGAGAACTTTTTTTTTT     A       A       AC      r       0       0       0       41      0       0       27      47      50      0       25      0       
5       164623412       AAAAAAAAAAGTTGCCACT     G       G       GT      f       0       0       41      0       1       0       15      46      50      0       4       0       

Column 1: Chromosome number
Column 2: Position of SNP
Column 3: (k-1)-mer (k = 20)
Column 4: Base of reference at the position of SNP
Column 5: Base of control
Column 6: Base of target
Column 7: Direction of k-mer sequence on the reference
Column 8: Number of k-mer with A at the end of k-mer in the control
Column 9: Number of k-mer with C at the end of k-mer in the control
Column 10: Number of k-mer with G at the end of k-mer in the control
Column 11: Number of k-mer with T at the end of k-mer in the control
Column 12: Number of k-mer with A at the end of k-mer in the target
Column 13: Number of k-mer with C at the end of k-mer in the target
Column 14: Number of k-mer with G at the end of k-mer in the target
Column 15: Number of k-mer with T at the end of k-mer in the target
Column 16: Number of reads in the control sort_uniq file with control type base
Column 17: Number of reads in the control sort_uniq file with target type base
Column 18: Number of reads in the target sort_uniq file with control type base
Column 19: Number of reads in the target sort_uniq file with target type base
Column 20: Genotype (M: homozygous, H: heterozygous)
```

## Detection of polymorphisms between target and control

  
SNPs between ERR3063487 and ERR3063486  
```
I       27950   A       T       9       31      0       0       17      M
I       892680  C       G       1       8       0       5       3       H
I       1196268 T       A       1       9       0       8       4       H
I       1380502 A       T       9       22      0       0       13      M
I       1728826 T       G       1       8       0       6       3       H
I       3203676 T       G       3       16      1       11      5       H
I       3407954 C       A       16      28      0       0       12      M
I       3656814 A       C       3       8       1       8       4       H
I       5001132 T       G       1       8       1       6       3       H
I       6324213 G       T       9       24      0       0       21      M
I       7249395 T       G       3       19      1       7       4       H
I       7263091 T       G       4       14      1       9       6       H
I       9136539 T       A       11      20      0       0       16      M
I       10137843        T       G       3       6       0       9       4       H
I       14097862        A       T       6       17      1       9       13      H
II      86891   A       G       2       7       0       8       4       H
II      271122  C       A       1       16      0       6       3       H
II      1179320 C       A       9       11      0       0       17      M
II      2500482 C       A       13      15      0       0       17      M
II      3552886 A       C       3       15      1       5       3       H
II      3624396 A       C       2       12      1       9       5       H
II      3648966 C       T       17      9       0       0       17      M
II      3771956 G       C       11      19      0       0       16      M
II      3935824 A       C       1       7       0       6       3       H
II      8284226 C       A       16      19      0       0       20      M
II      8553707 T       A       6       21      0       0       6       M
II      9410187 C       T       9       23      0       0       18      M
II      9937543 T       G       11      21      0       0       18      M
II      10629519        A       C       6       17      1       12      8       H
II      10685303        T       A       13      21      0       0       16      M
II      12732768        A       C       2       13      1       6       3       H
II      14096056        T       G       2       6       0       5       4       H
III     198231  T       A       9       16      0       0       18      M
III     3091906 T       G       3       19      0       9       4       H
III     3643248 G       C       1       6       0       6       3       H
III     4824486 C       G       15      33      0       0       23      M
III     7532164 T       G       1       8       1       6       3       H
III     9723566 G       A       17      25      0       0       23      M
III     11532166        C       T       6       15      0       1       18      M
III     13092063        C       T       5       13      0       6       6       H
III     13273147        A       T       1       15      1       6       4       H
III     13350315        A       C       3       14      0       7       4       H
IV      554740  A       C       5       11      1       9       4       H
IV      876681  A       G       1       20      0       7       4       H
IV      1159003 A       C       6       13      0       7       6       H
IV      2417073 G       A       17      32      0       0       21      M
IV      2858610 C       T       7       21      0       0       14      M
IV      3931877 G       A       8       13      0       0       20      M
IV      4298928 A       C       3       11      1       9       4       H
IV      5200355 G       T       13      27      0       0       23      M
IV      6481216 C       G       2       17      0       0       8       M
IV      6796145 C       T       11      16      0       0       19      M
IV      6967218 G       A       14      25      0       0       25      M
IV      8053747 T       A       14      23      0       0       19      M
IV      8951120 T       G       9       10      1       13      6       H
IV      9709645 C       A       17      22      0       0       27      M
IV      10195034        T       G       1       14      1       8       4       H
IV      13672183        C       T       10      18      0       0       18      M
IV      14217812        A       C       1       9       0       6       3       H
IV      14760376        T       G       3       9       0       5       3       H
IV      16870604        A       C       4       12      1       8       4       H
V       7011    T       G       3       10      1       11      5       H
V       974819  A       C       2       6       0       9       4       H
V       3052728 A       C       1       9       1       10      5       H
V       3277678 G       A       2       22      0       1       19      M
V       3786240 T       A       9       23      0       0       17      M
V       4261370 T       G       4       14      0       9       4       H
V       5134172 C       T       6       19      0       17      10      H
V       7816318 A       C       5       20      1       8       4       H
V       9771516 A       T       7       25      0       18      12      H
V       10513400        A       C       3       11      0       9       4       H
V       11310419        G       T       9       15      0       0       18      M
V       15880652        T       C       1       20      1       6       3       H
V       19657843        T       A       10      25      0       0       24      M
V       19718778        G       A       2       11      0       0       11      M
V       19733914        A       C       1       6       0       6       3       H
V       20413571        T       G       4       13      1       6       3       H
X       2843588 A       C       4       13      0       9       4       H
X       4194330 T       C       15      23      0       0       22      M
X       4247242 T       G       2       7       1       6       3       H
X       5446454 C       T       4       21      0       0       7       M
X       6994561 A       T       20      29      0       0       29      M
X       7299494 A       G       2       7       0       6       3       H
X       7586849 T       G       1       8       1       8       4       H
X       10486673        A       T       12      31      0       0       23      M
X       14544549        A       G       21      13      0       0       30      M
X       14632040        T       A       13      21      0       0       21      M
X       15395882        A       G       3       17      0       13      6       H
X       15815152        T       G       3       13      1       5       3       H
X       16861146        C       T       6       9       0       0       11      M
X       16964164        T       G       3       14      0       5       3       H
```
indels between ERR3063487 and ERR3063486  
```
I       834776  I       834746  f       deletion        -12     6       0       7       4       H       TCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT
I       1251597 I       1251573 f       deletion        -2      5       0       8       4       H       TGTGTGTGTGTGTGTGTGTGTGTGT
I       1412142 I       1411952 f       insertion       102     12      1       6       3       H       CCCCCCGCTGACCCCAAACCAATATCCCGTCAAAAAACGAAAATTCATATTTTTCTTAATCTACAGTAATCCTACAGTGCCCCTACA
I       1560682 I       1560674 f       deletion        -1      17      0       0       19      M       TTTTTTTT
I       1560682 I       1561254 r       inversion               13      0       0       6       M       TTAAAGGTGGTGTGGTCGAATTTTTTTT
I       2160919 I       2160898 f       deletion        -1      9       0       8       4       H       TTTTTTTTCAAAAAAAAAAAA
I       2384261 I       2384244 f       insertion       1       14      1       8       4       H       CAAAAAAAAAAAAAA
I       2715230 III     12981863        r       translocation           5       1       5       3       H       CGTATTGCACAGCACATTTGACGCGCAAAAT
I       5081495 I       5081478 f       deletion        -2      7       1       6       4       H       TTTTTTTTTTTTTTTTTT
I       8028077 I       8028066 f       deletion        -1      6       0       6       3       H       TTTTTTTTTTT
I       8028078 I       8028065 f       insertion       1       6       0       6       3       H       TTTTTTTTTTT
I       9825014 I       9825007 f       insertion       1       23      0       8       13      H       AAAAA
I       10622455        IV      9192034 f       translocation           5       1       6       3       H       GTTCAAATAAAAATATTTTTTT
I       10887954        I       10887958        f       deletion        -6      11      0       1       10      M       A
I       11005509        I       11005489        f       deletion        -1      5       1       0       5       M       AAATTTTTTTTTTTTTTTTT
I       12856424        I       12856415        f       deletion        -1      14      0       8       5       H       TTTTTTTTT
I       13734806        I       13734788        f       deletion        -1      5       0       6       4       H       TTTTTTTTTTTTTTTTTT
II      191614  II      191601  f       insertion       1       11      0       9       4       H       AAAAAAAAAAA
II      948033  II      948099  f       deletion        -69     20      0       1       15      M       AT
II      1777672 II      1770125 f       insertion       7506    8       1       9       4       H       ATGGTGAGTAGCCGGTAATTTCATAGTTATTGAAATTTGA
II      2361947 II      2361931 f       deletion        -1      10      1       0       15      M       TTTTTCTTTTTTTTTT
II      4463297 V       1000739 r       translocation           6       0       6       4       H       TTTCGATTTTCCAGAAAATCAAAAAAAAA
II      4895056 V       18778607        r       translocation           5       0       5       3       H       TTCTACGTTTTGCAATGTGTTTTTT
II      5473562 II      3675168 r       inversion               8       1       5       3       H       TTTTACTCAGTTATGTTTTTTTT
II      12746320        III     4520340 f       translocation           6       1       5       3       H       TGTAAAATTGTTTTTTTTT
II      13112130        V       20740040        f       translocation           6       0       7       4       H       AAAAAAAAACGCATGCATTTTTCG
II      13327324        II      13327697        r       inversion               6       1       6       5       H       TTTTGACACTTTTTAGTAATAAATGCAAAAAAAATCAACAAAAATAGACTAAACATTGTAAAAACTGTAAA
AACTAAGAGAAAAAAT
III     1663181 III     1663357 f       deletion        -209    6       1       7       5       H       TTTTTTCCAGAAATTAATATTTCTAGAAAAAT
III     2300741 X       15839886        r       translocation           19      0       11      6       H       TTAAAGGTGGAGTAGCGCCAGTGGGAAAATTGCTTTAAAACATGCCTATGGTACCACAATGACCAAATATCAT
III     2520715 X       97782   f       translocation           7       1       6       5       H       TATTTTTTCGCCATTTTTTTT
III     2985966 IV      876821  f       translocation           5       1       6       3       H       AAAAAAATTTTTTTTTT
III     3566956 III     3566962 f       deletion        -48     8       1       6       4       H       TCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT
III     6231151 III     6231140 f       deletion        -1      14      0       6       3       H       AAAAAAAAAAA
III     10402397        III     10402387        f       deletion        -2      5       1       7       4       H       TTTTTTTTTTT
III     11280814        X       2272016 r       translocation           12      0       6       3       H       TTTTTTTCAAAAAAAAAAAAA
III     12543896        III     12543907        f       deletion        -62     8       1       9       4       H       AAATTTCCGGAAAACATGCAAATTGCCAGAATTGAAAATTTCCGGCAAAT
III     13419473        III     13419443        f       insertion       6       5       1       10      5       H       TGTGTGTGTGTGTGTGTGTGTGT
IV      1786345 IV      1786337 f       deletion        -1      15      0       0       12      M       AAAAAAAA
IV      2112309 IV      2112287 f       deletion        -1      7       1       5       3       H       TTTTTTGTTTTTTTTTTTTTTT
IV      2314588 IV      2314573 f       deletion        -1      10      1       8       4       H       TTTTTTTTTTTTTTT
IV      2445289 IV      2445276 f       deletion        -1      15      0       5       3       H       TTTTTTTTTTTTT
IV      3192017 IV      3192001 f       deletion        -1      10      1       8       5       H       AAAAAAAAAAAAAAAA
IV      3192018 IV      3192000 f       insertion       1       10      0       7       4       H       AAAAAAAAAAAAAAAA
IV      3336297 IV      3336328 f       deletion        -31     21      1       0       23      M       A
IV      3867139 IV      3867116 f       deletion        -2      11      1       6       3       H       ATATATATATATATATATATATAT
IV      4399486 IV      4399441 f       insertion       2       8       0       10      6       H       GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA
IV      4489131 IV      4489122 f       deletion        -1      16      0       1       17      M       TTTTTTTTT
V       1027708 V       1027700 f       deletion        -39     9       1       10      5       H       GCCTATGGCCTACGCCTATGGCCTACGCCTATGGCCTACGCCTATG
V       1989616 V       1989642 f       deletion        -3      6       1       6       3       H       GATAAAAACTACTTGGATAAATGA
V       9026732 V       9026728 f       deletion        -3      6       1       7       4       H       ATTATT
V       11884926        II      3151004 r       translocation           5       1       6       5       H       GTGCCGAGTGCCGATCGGCACAATGTG
V       18672436        IV      2681633 r       translocation           7       1       5       3       H       GGGAAAATTGCTTTAAAACATGCCTATGGTACTACAA
V       18890234        V       18892184        r       inversion               6       1       6       3       H       TTTTTATTGAAAACTAGTATAAAAATATA
V       20123760        V       20123893        f       deletion        -150    12      0       7       4       H       GGGGTTCGAACCCCGG
X       99547   X       99523   f       deletion        -1      8       1       9       5       H       AAAAAAAATTTTTTTTTTTTTTTT
X       438457  X       438446  f       insertion       1       12      0       5       5       H       TTTTTTTTT
X       562164  X       562155  f       insertion       1       22      0       0       22      M       AAAAAAA
X       1522918 X       1522902 f       deletion        -1      6       1       5       3       H       AAAAAAAAAAAAAAAA
X       2498483 X       2498466 f       deletion        -1      6       1       5       3       H       TTTTTTTTTTTTTTTTT
X       3823840 X       3823828 f       deletion        -1      10      0       8       4       H       AAAAAAAAAAAA
X       5885390 X       5885377 f       deletion        -1      11      0       10      5       H       TTTTTTTTTTTTT
X       7312229 X       11435923        r       inversion               6       0       5       4       H       TATTCACCCCGTTCGACTGTGCAATGGGTTTAATCTATTCACTTTGTAAATCAAAGAATCGACGACCGCCTCCTGAA
X       10023790        III     7850798 r       translocation           5       0       6       3       H       ATATCAAAATTTCATTTTTTTT
X       14258766        III     303520  f       translocation           6       0       9       5       H       TCACAAAATTCTTTGGCCGCCCCAAGTGTCCTAACTCGAAG
```

## Author
Akio Miyao, Ph.D. miyao@affrc.go.jp  
Institute of Crop Science / National Agriculture and Food Research Organization  
2-1-2, Kannondai, Tsukuba, Ibaraki 305-8518, Japan  

## Version
Version 1.0

## Citing PED
Cite this article as: TBA

## License
NARO NON-COMMERCIAL LICENSE AGREEMENT Version 1.0

This license is for 'Non-Commercial' use of software for polymorphic edge detection (PED)

1. Scientific use of PED is permitted free of charge.
1. Modification of PED is only permitted to the person of downloaded and his colleagues.
1. The National Agriculture and Food Research Organization (hereinafter referred to as NARO) does not guarantee that defects, errors or malfunction will not occur with respect to PED.
1. NARO shall not be responsible or liable for any damage or loss caused or be alleged to be caused, directly or indirectly, by the download and use of PED.
1. NARO shall not be obligated to correct or repair the program regardless of the extent, even if there are any defects of malfunctions in PED.
1. The copyright and all other rights of PED belong to NARO.
1. Selling, renting, re-use of license, or use for business purposes etc. of PED shall not be allowed. For commercial use, license of commercial use is required. Inquiries for such commercial license are directed to ped_request@ml.affrc.go.jp.
1. The PED may be changed, or the distribution maybe canceled without advance notification.
1. In case the result obtained using PED in used for publication in academic journals *etc.,* please refer the publication of PED and/or acknowledge the use of PED in the publication. 

Copyright (C) 2017 [National Agriculture and Food Research Organization](https://www.naro.affrc.go.jp/english/index.html). All rights reserved. 
