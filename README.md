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

### Installation

- Programs run on Unix platforms (Linux, MacOS, FreeBSD). At least 2 TB disk space is required. For analysis by *k*-mer method, computer cluster and scheduler are required. Programs for bidirectional method can be run on a computer, but using computer cluster is recommended.
- Download zip file of PED from https://github.com/akiomiyao/ped and extract.  
or  
% git clone https://github.com/akiomiyao/ped.git
-  To make reference data of *Arabidopsis thaliana,*  
    % cd TAIR10  
    % sh setup.sh  
-  To make reference data of rice,  
    % cd IRGSP1.0  
    % sh setup.sh  
-  To make reference data of hg38,  
    % cd hg38  
    % sh setup.sh  
- For soybean and mouse, execute shetup.sh in Gmax275v2.0 and GRCm38 directory, respectively.
- 'wget' is required. If your machine do not have wget program, install wget from package before run setup.sh. Making reference data of human genome takes two days. Making the data sets of reference genome is only required at the first usage of PED. 

## Demonstration of bidirectional alignment method

- To download sequence data, fastq-dump from NCBI is required.
    Tool kit can be download from
    https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/ 

    Because data of human is big, it takes long time for downloading from SRA. To check the performance, data of *Arabidopsis* is recommended.
    
### *Arabidopsis thaliana*
- Short read data SRR8181712 is used for demonstration.
- % perl download.pl SRR8181712  
    SRR8181712_1.fastq will be downloaded in SRR8181712/read directory.
-  To make sorted unique sequence,  
    % perl sort_uniq.pl SRR8181712
-  Because the bidirectional alignment method can detect not only SNP but also insertion, deletion, inversion and translocation, we recommend the bidirectional alignment method for usual analysis.
    For stand alone computer,  
    % perl align.pl target_name reference margin tmpdir  
    'tmpdir' is optional.    
    *e.g.*  
    % perl align.pl SRR8181712 TAIR10 0  
    % perl align.pl SRR8181712 TAIR10 5  
    % perl align.pl SRR8181712 TAIR10 10  
    % perl align.pl SRR8181712 TAIR10 15  
    At the first step of the method, *k*-mer (*k* = 20) sequences from both ends of short read will be mapped to the reference genome. If the *k*-mer sequence is repetitive to the genome sequence, the alignment will be failed. If the margin is given, *k*-mers begin inside of the margin. This will increase the coverage of polymorphisms.  


    For the computer cluster,  
    % qsub -v target=target_name,ref=reference,margin=0,tmpdir=path_of_tmpdir align.pl  
    *e.g.*  
    % qsub -v target=SRR8181712,ref=TAIR10,margin=0 align.pl  
    % qsub -v target=SRR8181712,ref=TAIR10,margin=5 align.pl  
    % qsub -v target=SRR8181712,ref=TAIR10,margin=10 align.pl  
    % qsub -v target=SRR8181712,ref=TAIR10,margin=15 align.pl  
-  From result of align.pl, snp and indel data are gathered and then split for verification.  
    % perl split_snp.pl SRR8181712  
    % perl split_indel.pl SRR8181712  
    or  
    % qsub -v target=SRR8181712 split_snp.pl  
    % qsub -v target=SRR8181712 split_indel.pl  

- To verify SNPs,  
    % perl verify_snp.pl target control reference number type tmpdir  
    *e.g.*  
    % perl verify_snp.pl SRR8181712 default TAIR10 01 bi  
    or  
    % qsub -v target=SRR8181712,control=default,ref=TAIR10,number=01,type=bi verify_snp.pl  
    
    target : target name, *e.g.* SRR8181712  
    control : *e.g.* SRR7686247 or 'default', if you want to use reference data for control  
    ref : reference genome, *e.g.* TAIR10  
    number : specify number of splited subfile  
    type : specify type of input data (bi, kmer or vcf)  
    tmpdir : specify temporary directofy on local disk (can be omitted)  
- To verify indel,  
    % perl verify_indel.pl target control reference number type tmpdir  
    *e.g.*  
    % perl verify_indel.pl SRR8181712 default TAIR10 01 bi
    or  
    % qsub -v target=SRR8181712,control=default,ref=TAIR10,number=01,type=bi verify_indel.pl 

### Human
- Short read data ERR194147 is used for demonstration. ERR194147 is a short reads of whole genome sequence from a member of the [Coriell CEPH/UTAH 1463 family](https://www.coriell.org/0/Sections/Collections/NIGMS/CEPHFamiliesDetail.aspx?PgId=441&fam=1463&) for a "[platinum](https://www.illumina.com/platinumgenomes.html)" standard comprehensive set for variant calling improvement.
- % perl download.pl ERR194147  
    ERR194147_1.fastq and ERR194147_2.fastq will be downloaded in ERR194147/read directory.
- To make sorted unique sequence,  
    % perl sort_uniq.pl ERR194147
- Because the bidirectional alignment method can detect not only SNP but also insertion, deletion, inversion and translocation, we recommend the bidirectional alignment method for usual analysis.  
    For stand alone computer,  
    % perl align.pl target_name reference margin tmpdir  
    *e.g.*  
    % perl align.pl ERR194147 hg38 0 /mnt/ssd  
    % perl align.pl ERR194147 hg38 5 /mnt/ssd  
    % perl align.pl ERR194147 hg38 10 /mnt/ssd  
    % perl align.pl ERR194147 hg38 15 /mnt/ssd  
    At the first step of the method, *k*-mer (*k* = 20) sequences from both ends of short read will be mapped to the reference genome. If the *k*-mer sequence is repetitive to the genome sequence, the alignment will be failed. If the margin is given, *k*-mers begin inside of the margin. This will increase the coverage of polymorphisms.  
    /mnt/ssd (tmpdir, specify to a directory on the fast local disk) is optional.  
    
    For the computer cluster,  
    % qsub -v target=target_name,ref=reference,margin=0,tmpdir=path_of_tmpdir align.pl  
    *e.g.*  
    % qsub -v target=ERR194147,ref=hg38,margin=0,tmpdir=/mnt/ssd align.pl  
    % qsub -v target=ERR194147,ref=hg38,margin=5,tmpdir=/mnt/ssd align.pl  
    % qsub -v target=ERR194147,ref=hg38,margin=10,tmpdir=/mnt/ssd align.pl  
    % qsub -v target=ERR194147,ref=hg38,margin=15,tmpdir=/mnt/ssd align.pl  
    'tmpdir' on local disk for qsub makes speed up the processing. More than 1 TB SSD of local disk for tmpdir is recommended.
- From result of align.pl, snp and indel data are gathered and then split for verification.  
    % perl split_snp.pl ERR194147  
    % perl split_indel.pl ERR194147  
    or  
    % qsub -v target=ERR194147,tmpdir=/mnt/ssd split_snp.pl  
    % qsub -v target=ERR194147 split_indel.pl  
    tmpdir for qsub run of split_snp.pl is optional.  
- To verify SNPs,  
    % perl verify_snp.pl target control reference number type tmpdir  
    *e.g.*  
    % perl verify_snp.pl ERR194147 default hg38 01 bi /mnt/ssd  
    or  
    % qsub -v target=ERR194147,control=default,ref=hg38,number=01,type=bi,tmpdir=/mnt/ssd verify_snp.pl  

  target : target name, *e.g.* ERR194147  
  control : *e.g.* ERR194146 or 'default', if you want to use reference data for control  
  ref : reference genome, *e.g.* hg38  
  number : specify number of splited subfile  
  type : specify type of input data (bi, kmer or vcf)  
  tmpdir : specify temporary directofy on local disk (can be omitted)  

A part of verify_snp.pl result is  
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
- To verify indel,  
    % perl verify_indel.pl target control reference number type tmpdir  
    *e.g.*  
    % perl verify_indel.pl ERR194147 default hg38 01 bi /mnt/ssd  
    or  
    % qsub -v target=ERR194147,control=default,ref=hg38,number=01,type=bi,tmpdir=/mnt/ssd verify_indel.pl
- A part of verify_indel.pl result is
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
## Demonstration of *k*-mer method
- Grid engine, *e.g.* Torque, is required for *k*-mer method.
- At first, all *k*-mers (*k* = 20) on each position of short read sequence from the sort_uniq data. Because the sort_uniq data is too big to handling, the sort_uniq data is split to sub files. Usually, the control is reference sequence data. To create *k*-mer data of control,  
    *e.g.*  
    % perl split_sort_uniq.pl hg38  
    or  
    % qsub -v target=hg38 split_sort_uniq.pl  
- To create count data of *k*-mer,  
    *e.g.*  
    % perl count_qsub.pl hg38 /mnt/ssd  
    count_qsub.pl is a launcher of count.pl for all sort_uniq subfiles.  
    The second argument is temporally directory of local disk on the node. This can be omitted, but encouraged.  
- Output of count.pl is separated by subfiles. merge.pl merges separated data.  
    *e.g.*  
    % perl merge_qsub.pl hg38 /mnt/ssd  
    merge_qsub.pl is a launcher of merge.pl for all count subfiles.  
- merged files are converted to last-base-count data.  
    *e.g.*  
    % perl lbc_qsub.pl hg38 /mnt/ssd  
    lbc_qsub.pl is a launcher of last_base_count.pl for all merged subfiles.  
- To create lbc files of target (e.g ERR194147),  
    % perl split_sort_uniq.pl ERR194147  
    % perl count_qsub.pl ERR194147 /mnt/ssd  
    % perl merge_qsub.pl ERR194147 /mnt/ssd  
    % perl lbc_qsub.pl ERR194147 /mnt/ssd  
- SNP detection between target(ERR194147) and control(hg38).  
    % perl snp_qsub.pl ERR194147 hg38 /mnt/ssd 10  
    snp_qsub.pl is a launcher of snp.pl for detection of SNP at the last base of *k*-mer between target and control.
- To map SNPs,  
    % perl map_qsub.pl ERR194147 hg38 /mnt/ssd  
- To verify mapped SNPs,  
    % perl verify_snp_qsub.pl ERR194147 default hg38 kmer /mnt/ssd
- A part of verify_snp.pl result is
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
