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
  e.g. perl mkref.pl dmel626  
  Directory dmel626 for reference of Drosophila melanogaster release 6.26
  will be created.  
- % perl download.pl accession  
  e.g. perl download.pl SRR5989890  
  Data directory of SRR5989890 will be created.  
  Fastq data will be downloaded in SRR5989890/read subdirectory.
- % perl bidirectional.pl target control reference  
  e.g. perl bidirectional.pl SRR5989890 default dmel626  
    After two hours, you will find results in SRR5989890 directory.  
    SRR5989890.indel is list of structural variation.  
    SRR5989890.snp is list of SNPs.  
    SRR5989890.snp.vcf is the vcf file for SNPs.  
    If you want to detect polymorphisms between target and control (e.g.SRR8182352),  
    set SRR8182352 as the control.  
    At first, bidirectional method detects polymorphisms between target and reference.  
    The detected polymorphisms are verified by counting target and control reads with/without polymorphisms.  
    If you want to use the control, downloading control sequeces and making the sort_uniq data are required, before the run the bidirectional.pl.  
    For example,  
    % perl download.pl SRR8182352  
    % perl sort_uniq.pl SRR8182352  
- To confirm the alignment for detected polymorphisms,  
  % perl search.pl target chr position  
  e.g. % perl search.pl SRR5989890 2L 15920731  
  Alignments will be selected by the search script.  
- if you want to run with computer cluster,  
  % perl qsub_bidirectional.pl SRR5989890 default dmel626  
- Run without arguments, help for script will be shown.  

## Simple instruction for kmer method
- % perl kmer.pl target control reference  
- % perl kmer.pl SRR8181712 default TAIR10  
- if you want to run with computer cluster,  
  % perl qsub_kmer.pl SRR8181712 default TAIR10

## Making reference data sets
- % perl mkref.pl  
    Without specify reference, reference name and description are listed.  
    If you want to make new reference, add the information about the reference to 'config' file.  
- To run by computer cluster,  
  % qsub -v target=target mkref.pl
- For example,  
  % perl mkref.pl hg38  
    Data set of human genome hg38 will be made. It takes about two days.  
  % perl mkref.pl dmel626  
    Data set of Drosophila melanogaster r6.26 will be made. It takes one hour.  
  % qsub -v target=dmel626 mkref.pl

## Set up data directory
- If you want to download short read sequence from NCBI SRA,  
  % perl download.pl accession  
  For example,  
  % perl download.pl SRR5989890  
    Directory SRR5989890 will be made, fastq files will be downloaded to  
    SRR5989890/read directory.
- If you want to analyze local file,  
  % mkdir mydata1  
  % mkdir mydata1/read  
  % cp somewhere/mydata1.fastq mydata1/read  
  % perl bidirectional.pl mydata1 reference

  target : target name, *e.g.* ERR194147  
  control : *e.g.* ERR194146 or 'default', if you want to use reference data for control  
  referene : reference genome, *e.g.* hg38  
  tmpdir : specify temporary directofy on local disk (can be omitted).
           tmpdir is specified as fourth argument of scrpit.

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
Our verify process counts reads containing polymorphic region.  
Basically, counts are from target and reference.  
If control is specified, counts from target and control will be listed.  
For example,  
ERR3063488 is read data of Caenorhabditis elegans mutant.  
ERR3063489 is read data of Caenorhabditis elegans wild-type.  
% perl bidirectional.pl ERR3063488 ERR3063489 WBcel235  
retruns 'M' or 'H' marked SNPs and indels of ERR3063488 which are absent in ERR3063489.  

I appreciate short read release to SRA by ENS Lyon.  
  
SNPs between ERR3063488 and ERR3063489  
```
I       332222  T       G       1       11      0       5       3       H
I       1448800 G       A       2       5       1       9       4       H
I       6283523 A       C       1       12      1       6       3       H
I       6415699 A       C       5       17      1       13      7       H
I       8167808 T       A       10      22      0       0       13      M
I       8702729 T       C       21      21      0       0       28      M
I       10509827        G       A       3       21      0       0       18      M
I       10885278        A       C       3       7       0       6       3       H
I       11533386        A       C       5       18      1       10      5       H
I       11600291        A       C       2       12      1       5       3       H
I       12219825        T       G       2       11      0       7       5       H
I       12595119        T       G       2       14      0       9       4       H
I       12902496        A       C       2       9       1       6       3       H
I       14207159        A       C       6       12      1       11      7       H
I       14208132        A       C       4       7       0       8       5       H
I       14403748        G       A       5       15      0       0       9       M
I       14457330        A       C       1       14      1       7       4       H
I       14668351        T       G       1       9       1       5       4       H
II      3941726 A       C       2       11      0       9       4       H
II      5996649 A       T       10      23      0       0       18      M
II      6900518 T       G       2       9       1       7       4       H
II      7801318 A       C       4       9       0       9       4       H
II      7880040 A       G       13      10      0       0       17      M
II      10529763        T       G       1       13      1       9       4       H
II      10836830        C       A       2       10      0       5       3       H
II      11490799        C       T       10      19      0       0       16      M
II      11664180        T       C       2       8       1       5       3       H
II      13072968        A       C       1       5       0       10      5       H
II      13315213        A       T       6       14      1       6       4       H
II      13359250        A       T       12      14      0       0       21      M
II      14004296        A       C       3       19      0       6       3       H
II      15240134        T       C       1       11      1       5       3       H
III     2056209 T       C       8       20      0       13      11      H
III     3065696 A       C       1       7       0       6       3       H
III     4176111 T       C       12      25      0       0       21      M
III     6896648 C       A       1       12      0       5       3       H
III     6932472 T       G       16      19      0       0       27      M
III     6972389 T       G       12      20      0       0       18      M
III     8190132 T       G       2       10      1       9       4       H
III     8843993 A       T       12      7       0       0       17      M
III     9052039 A       C       4       10      1       6       3       H
III     10961250        T       C       3       10      0       8       4       H
III     11051821        C       T       17      17      0       0       21      M
III     11557788        T       G       3       9       0       9       6       H
III     12099100        T       G       2       22      0       9       4       H
III     12508403        C       A       18      14      0       0       27      M
IV      448629  T       G       16      22      0       0       24      M
IV      1001401 A       T       12      13      0       0       21      M
IV      1372875 T       G       2       16      1       9       4       H
IV      1476073 A       C       3       12      0       6       4       H
IV      1781639 C       T       7       17      0       0       11      M
IV      1845669 C       G       1       16      0       0       18      M
IV      1991053 G       T       12      11      0       0       16      M
IV      3626736 T       A       9       14      0       0       17      M
IV      3660762 T       G       4       17      0       7       4       H
IV      3990444 A       C       4       17      0       8       4       H
IV      6946911 C       G       14      18      0       0       20      M
IV      7957975 T       G       3       7       0       5       3       H
IV      9594931 G       A       9       15      0       0       14      M
IV      10545314        G       A       14      19      0       0       22      M
IV      12625239        A       C       11      11      0       0       12      M
IV      12804123        T       C       1       13      0       8       4       H
IV      13058123        C       T       12      16      0       0       14      M
IV      15914692        C       A       17      17      0       0       25      M
IV      16388039        A       C       4       11      0       9       4       H
IV      16667859        A       C       1       7       0       7       4       H
IV      17009515        G       A       9       27      0       0       22      M
V       8959    A       C       1       12      0       5       4       H
V       1317773 A       C       2       14      0       9       4       H
V       2028980 T       G       14      11      0       0       19      M
V       4646273 T       C       1       11      0       7       4       H
V       6100769 T       G       2       11      0       11      5       H
V       8259099 C       T       14      27      0       0       19      M
V       10100458        A       C       3       14      1       11      5       H
V       11871371        A       T       10      24      0       0       17      M
V       13536894        T       C       15      19      0       1       23      M
V       14563775        T       C       4       14      0       6       3       H
V       15197495        A       C       3       18      0       11      5       H
V       15790839        C       T       1       20      0       0       12      M
V       16365541        A       G       9       12      0       0       17      M
V       16726243        T       G       3       15      1       10      5       H
V       16967163        T       G       2       12      1       7       4       H
V       18144888        C       A       11      19      0       0       22      M
V       18999888        C       A       14      26      0       0       26      M
V       19044283        T       G       3       12      1       6       3       H
X       1521358 T       G       7       15      0       0       15      M
X       1911051 T       G       4       11      1       5       3       H
X       2526842 T       G       1       17      0       7       4       H
X       3455431 T       G       3       16      0       5       3       H
X       4554851 C       T       13      25      0       0       19      M
X       5458200 A       C       1       12      1       6       3       H
X       6988370 C       T       18      14      0       0       23      M
X       7536774 C       A       5       17      0       0       13      M
X       11154570        G       A       17      29      0       0       17      M
X       11576895        A       G       1       9       1       6       4       H
X       11577936        T       G       7       16      0       8       4       H
X       11732798        C       A       7       8       0       0       9       M
X       12180689        A       C       3       13      1       6       4       H
X       13358034        C       A       14      15      0       0       22      M
X       13672823        A       C       1       10      0       5       3       H
X       14693495        A       G       16      19      0       0       23      M
X       15361048        A       C       2       14      0       8       4       H
X       16244315        C       G       10      14      0       0       9       M
X       17007533        G       A       9       18      0       0       19      M
X       17030601        T       G       2       5       1       6       3       H
X       17487527        G       A       19      28      0       0       29      M
```
indels between ERR3063488 and ERR3063489  
```
I       1411098 I       1411130 f       deletion        -33     18      0       0       13      M       
I       1946171 I       1944320 f       insertion       1858    5       0       5       8       H       TAAAATTC
I       1968845 X       15237407        f       deletion        -13268601       7       1       6       4       H       GCCTACTTTCTGGCGCGAAAATAGCGGCAACAGAGAGA
I       10376611        I       10376593        f       deletion        -2      5       1       7       4       H       ATATATATATATATATATA
I       11506461        I       11506451        f       deletion        -1      13      0       8       4       H       TTTTTTTTTT
I       13101741        I       13101724        f       deletion        -1      8       1       5       3       H       TTTTTTTTTTTTTTTTT
I       13110447        I       13110435        f       deletion        -1      8       1       6       3       H       AAAAAAAAAAAA
II      2567493 II      2568128 f       deletion        -679    11      1       11      5       H       AAATGCAATTTTTAACGAAAATTTGTCAATTTTTCGATTAAAA
II      4780838 II      4780830 f       deletion        -3      6       0       6       4       H       TTTTTTTTTT
II      4895056 V       18778607        r       inversion               7       1       5       3       H       TTCTACGTTTTGCAATGTGTTTTTT
II      5403699 IV      16418046        r       inversion               5       0       5       3       H       AAAAAAACGAAAAAAAAA
II      10031763        II      10031748        f       deletion        -2      6       1       6       4       H       TTTTTTTTTTTTTTTT
II      13002251        V       1414357 r       inversion               8       1       8       8       H       AATTTTTTTTGTTCGACTTCCAAAA
II      13223262        V       3742448 r       inversion               9       1       9       4       H       ATTTGCCCATTTGCCGAAAAAAAAA
II      15200675        II      15200789        f       deletion        -109    22      0       0       21      M       CCCCCC
III     1753031 III     1753017 f       insertion       1       5       1       6       4       H       TTTTTTTTTTTT
III     3115212 III     3115193 f       insertion       2       6       1       5       3       H       AAAAAAAAAAAAAAAA
III     3145001 III     3144991 f       deletion        -2      6       0       5       3       H       TTTTTTTTTTT
III     7157630 V       20500233        f       deletion        -13342624       11      0       5       3       H       ATATATTTCTCAAAAAAAAA
III     8500712 III     8500689 f       insertion       2       5       0       6       3       H       GAGAGAGAGAGAGAGAGAGA
III     9066341 III     9066324 f       deletion        -2      7       1       7       4       H       ATATATATATATATATAT
III     11530173        IV      13391387        r       inversion               10      1       5       3       H       TTTTCAAAAATCTAATTTTGTTC
III     11839745        III     11839727        f       deletion        -2      10      1       8       4       H       TATATATATATATATATAT
III     11943870        III     11943862        f       insertion       2       9       1       7       8       H       AAAAA
III     12473776        II      1240280 f       insertion       11233463        8       1       6       3       H       CTAGGCCTAAGAATAAGCCTAAGCCTAAGCCT
IV      197758  IV      197740  f       insertion       1       23      0       5       10      H       AAAAAAATTTTTTTTT
IV      1015325 IV      1015310 f       deletion        -1      5       1       10      5       H       AAAAAAAAAAAAAAA
IV      3170651 IV      3170636 f       insertion       1       7       0       8       4       H       AAAAAAAAAAAAA
IV      3182386 IV      3182367 f       deletion        -19     17      0       10      11      H       ATAGCTAGGTGCCTATCTCATACCTAGGTGCCTACCT
IV      3225586 IV      3225570 f       deletion        -1      5       1       5       7       H       AAAAAAAAAAAAAAAA
IV      12045733        IV      12045717        f       insertion       1       11      1       9       4       H       AAAAAAAAAAAAAA
IV      14749269        IV      14749252        f       deletion        -2      13      0       12      6       H       ATATATATATATATATAT
IV      16927209        IV      16927194        f       insertion       1       11      0       8       4       H       TTTTTTTTTTTTT
V       1348255 V       20002829        r       inversion               7       0       8       4       H       TAAGCCTAAGCCTAAGCCTGAGCC
V       1646502 X       2286086 f       deletion        -639610 6       0       9       4       H       AATTGTCTGAAAACATCGAATTTCA
V       5181873 V       5181871 f       insertion       9       29      0       24      21      H       TTTTTCGT
V       5956837 V       5957448 f       deletion        -614    22      0       0       23      M       AA
V       6499662 IV      2007203 f       insertion       4492438 5       1       5       4       H       TTATTGTAAAATATCTAAAA
V       7592908 V       7592888 f       deletion        -2      7       1       6       3       H       TATATATATATATATATATAT
V       14946423        V       14946612        f       deletion        -191    20      0       0       20      M       T
V       15842476        V       15855565        f       deletion        -13178  10      0       6       4       H       TTTTTTTAATCTATAAACCACACATTTTGAGCAATCAATTTTGCGTCTTTTTGATCAGGAAGATTCACAAA
TGCAGGAGCAGGGAACG
V       17966942        V       17966926        f       deletion        -1      6       0       5       3       H       AAAAAAAAAAAAAAAA
V       19879865        V       19792905        r       inversion               5       1       9       4       H       TGTCAAAGATTCCTGATGTCAAAATGAATGCAAGTGAGAAATCCAAAGTGATTGAGGAG
X       1058707 IV      5697577 r       inversion               8       1       5       3       H       TTTATCAAAAACTTTTTTTTTT
X       2335393 X       2335378 f       deletion        -1      5       1       8       8       H       TTTTTTTTTTTTTTT
X       7069115 IV      6306975 r       inversion               10      0       6       6       H       CTCTATGCAAAAAGCGATCAATGTCGTCACTGATTGGTGTGCAAAGTG
X       7738629 X       1734175 r       inversion               8       1       6       4       H       GTTCAAATCAGATTTTAAGAACACGACTTGCGGTCCGC
X       9034324 IV      199585  f       insertion       8834720 5       0       7       4       H       TAAAATACGTTTTTTTTT
X       14381618        IV      14919917        r       inversion               5       0       7       4       H       AGCTCATCTACACTGTGAGCAAATTTGCATTGCTCCC
X       14719283        V       18041043        r       inversion               8       0       6       3       H       GGAAAATTGCCTTTTTTCCGGCAACTTCGGC
X       14934806        I       710368  f       insertion       14224394        9       0       6       3       H       ATTAGACTCAAAATTGTCTGAAAACACCAAATTTCATAATGAA
X       16825915        X       16825903        f       deletion        -2      5       0       5       3       H       AAAAAAAAAAAAA
X       16872379        X       16872368        f       insertion       1       8       0       1       16      M       TTTTTTTTT
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
