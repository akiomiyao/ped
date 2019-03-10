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
- % perl mkref.pl dmel626  
- % perl download.pl SRR5989890
- % perl bidirectional.pl SRR5989890 default dmel626  
    After two hours, you will find results in SRR5989890 directory.  
    SRR5989890.indel is list of structural variation.  
    SRR5989890.snp is list of SNPs.  
    SRR5989890.snp.vcf is the vcf file for SNPs.  
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
- For example,  
  % perl mkref.pl hg38  
    Data set of human genome hg38 will be made. It takes about two days.  
  % perl mkref.pl dmel626  
    Data set of Drosophila melanogaster r6.26 will be made. It takes one hour.

## Set up data directory
- If you want to download short read sequence from NCBI SRA,  
  % perl download.pl Accession  
  For example,  
  % perl download.pl SRR5989890  
    Directory SRR5989890 will be made, fastq files will be downloaded to  
    SRR5989890/read directory.
- If you want to analyze local file,  
  % mkdir mydata1  
  % mkdir mydata1/read  
  % cp somewhere/mydata1.fastq mydata1/read  
  % perl bidirectional.pl mydata1 Reference

  target : target name, *e.g.* ERR194147  
  control : *e.g.* ERR194146 or 'default', if you want to use reference data for control  
  referene : reference genome, *e.g.* hg38  
  tmpdir : specify temporary directofy on local disk (can be omitted).
           tmpdir is specified as fourth argument of scrpit.

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
