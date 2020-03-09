# PED : Polymorphic Edge Detection

Polymorphic Edge Detection (PED) is the analysis flow for DNA polymorphism detection from short reads of next generation sequencer (NGS). We developed two methods to detect polymorphisms based on detection of the polymorphic edge. One is based on bidirectional alignment and the other is based on comparison of *k*-mers. Examples of PED result and useful information are shown in the [wiki](https://github.com/akiomiyao/ped/wiki) and [Web page (English)](https://akiomiyao.github.io/ped/) [(Japanese)](https://akiomiyao.github.io/ped/index_ja.html).

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

## Simplified instruction
- The ped.pl is a multithreaded script, suitable for the multi-core CPU like as 16 or 8 cores.  
  Of course, the ped.pl can run with the 2 or single core machine, but slow.  
  The ped.pl runs on Linux (or FreeBSD) machine and Mac with at least 4 GB RAM and 1 TB hard disk (or SSD).  

  Following is a demonstration of spontaneous SNPs and SVs detection from a *Caenorhabditis elegans* with 250-times repeated generations.  
```
cd ped
perl download.pl accession=ERR3063486
perl download.pl accession=ERR3063486
perl ped.pl target=ERR3063487,control=ERR3063486,ref=WBcel235
```
  Installation of fastq-dump and ped scripts is described below.  
  The docker container for Linux includes fastq-dump and ped scripts.  
```
docker run -v `pwd`:/work -w /ped akiomiyao/ped perl download.pl accession=ERR3063486,wd=/work
docker run -v `pwd`:/work -w /ped akiomiyao/ped perl download.pl accession=ERR3063487,wd=/work
docker run -v `pwd`:/work -w /ped akiomiyao/ped perl ped.pl target=ERR3063487,control=ERR3063486,ref=WBcel235,wd=/work
```
-  ERR3063487 sequence is after 250 generations of the nematoda (ERR3063486).  
  Downloading fastq files may take several hours, because connection of fastq-dump to NCBI-SRA is slow.  
  Sometimes, download.pl returns the timeout of network connection. In the case, network will be reconnected and redumed the download.  
  Fastq files will be saved in ERR3063486/read and ERR3063487/read.  
  SNPs and SVs in ERR3063487 against ERR3063486, *i.e.* spontaneous mutations during 250 generations, will be saved in ERR3063487 directory.  
  If control is ommitted, polymorphisms against reference genome will be saved in target directory.  
  If script runs without arguments, description of how to use the script will be shown.  
- Options,  
  thread=8 : specify the max thread number. Default is the number of logical core or 14.  
  tmpdir=/mnt/ssd : specify the temporally directory to /mnt/ssd. Default is target directory.  

## Installation
- If you do not want to use the docker container, downloading of programs is required.  
- Programs run on Unix platforms (FreeBSD or Linux) and Mac.  
- Download zip file of PED from https://github.com/akiomiyao/ped and extract.  
or  
```
git clone https://github.com/akiomiyao/ped.git  
```
If you got scripts by clone command of git, update to newest version is very easy using pull command of git.  
```
git pull  
```
- To download sequence data, fastq-dump from NCBI is required.  
    Tool kit can be download from  
    https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/ 
- To download reference data, curl is required.  
    If your machine do not have curl program, install curl from package.  
```
sudo apt-get install curl (Ubontu)
sudo yum install curl (CentOS)
sudo pkg install curl (FreeBSD)
```

## Setup of Docker
In the case of docker, zombie processes due to execution of sub process will be incleased.  
When the ped analysis is finished, zombie processes will be removed.   
On the run of docker container, the --init option is effective to kill zombie processes.  
But premature termination of sort command is observed with --init options.  
If the premature termination is observed, direct run of ped script from the github instead of docker is recommended.

https://docs.docker.com/install/linux/docker-ce/ubuntu/
```
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
```
To get or update the container,
```
sudo docker pull akiomiyao/ped
```
To check running containers,
```
sudo docker stats
```
To kill running container,
```
sudo docker kill Container_ID
```
If you want run the docker container without sudo or su,
```
sudo usermod -a -G docker your_username
```
After the new login, docker commands can be execute with your account.  

## Instruction for bidirectional method  
```
perl ped.pl target=target,control=control,ref=reference
```
- Run without argument, ped.pl returns the usage and options.  
If reference data is absent, ped.pl downloads the reference sequence and makes reference dataset.  
```
perl ped.pl target=ERR3063487,control=ERR3063486,ref=WBcel235  
```
- ERR3063487 sequence is after 250 generations of the nematoda (ERR3063486).  
  SNPs and SVs in ERR3063487 against ERR3063486, *i.e.* spontaneous mutations during 250 generations, will be saved in ERR3063487 directory.  
  If control is ommitted, polymorphisms against reference genome will be saved in target directory.  
  If script runs without arguments, description of how to use the script will be shown.  
- If you want to make reference data separately,  
```
perl mkref.pl reference  
```
For example,  
```
perl mkref.pl WBcel235  
```
or  
```
perl qsub_mkref.pl WBcel235 (For computer cluster)  
```
- Directory WBcel235 for reference of *Caenorhabditis elegans* WBcel235 will be created.  
If run without argument, help and suported reference will be listed.  
If you want to new reference, add the reference information to the config file.    
Format is described in the comment in config file.  
To make reference of human,  
```
perl mkref.pl hg38  
```
- Currently, supporting reference genomes are  
```
  Reference Name  Description  
  B73v4           Corn (Zea mays B73) RefGen v4
  Bomo            Silkworm (Bombyx mori) Genome assembly (Nov.2016)
  GRCm38          Mouse (Mus musculus) Genome Reference Consortium Mouse Build 38
  Gmax275v2.0     Soybean (Glycine max) genome project assemble version 2
  IBSC2           Barley (Hordeum vulgare L. cv. Molex) Release 44
  IRGSP1.0        Rice (Olyza sativa L. cv. Nipponbare) version 1.0
  IWGSC1.0        Wheat (Triticum aestivum L. cv. Chinese Spring) Version 1.0
  LJ3             Lotus japonicus MG20 v3.0
  SL3             Tomato (Solanum lycopersicum cv. Heinz 1706) Build 3.0
  SScrofa11.1     Pig (Sus scrofa) Release-97
  TAIR10          Arabidopsis thaliana version TAIR10
  UMD3.1          Cow (Bos taurus L1 Dominette 01449) UMD 3.1
  WBcel235        Caenorhabditis elegans WBcel235
  danRer11        Zebrafish (Danio rerio) Genome Reference Consortium Zebrafish Build 11
  dmel626         Drosophila melanogaster
  hg19            Human (Homo sapiens) Genome Reference Consortium Human Build 19
  hg38            Human (Homo sapiens) Genome Reference Consortium Human Build 38
```
  If fetch the fasta file is failed by the script, fetch the file separately and save in the reference directory and run the script.  
- Otherwise,  
```
perl mkref.pl reference fasta_file_name  
```
For example,  
```
mkdir IRGSP1.0  
cp somewhere/IRGSP-1.0_genome.fasta.gz IRGSP1.0  
perl mkref.pl IRGSP1.0  IRGSP-1.0_genome.fasta.gz  
```
or  
```
qsub -v target=IRGSP1.0,file=IRGSP-1.0_genome.fasta.gz mkref.pl  
```
- If you want to analyze public data in SRA.  
```
perl download.pl accession=accession  
```
For example,  
```
perl download.pl accession=ERR3063486  
perl download.pl accession=ERR3063487   
```
- Data directory of ERR3063486 and ERR3063487 will be created.  
  Fastq data from SRA in NCBI will be downloaded in read subdirectory.  
  ERR3063486 is the read data of *Caenorhabditis elegans* wild-type.  
  ERR3063487 is the read data of *Caenorhabditis elegans* mutant.  
  BioPoject https://www.ncbi.nlm.nih.gov/bioproject/PRJEB30822  
- If you want to analyze your private file,  
```
mkdir mydata1  
mkdir mydata1/read  
cp somewhere/mydata1.fastq mydata1/read  
perl bidirectional.pl mydata1 control reference  
```
or  
```
  perl ped.pl target=mydata1,control=control,ref=reference
```
- perl bidirectional.pl target control reference  
  For example,  
```
perl ped.pl target=ERR3063487,ref=WBcel235  
```
or  
```  
perl bidirectional.pl ERR3063487 default WBcel235  
```
or  
```
qsub -v target=ERR3063487,control=default,ref=WBcel235 bidirectional.pl 
```
-  After four hours, you will find results in ERR3063487 directory.  
   ERR3063487.sv is the list of structural variation.  
   ERR3063487.bi.snp is the list of SNPs.  
   ERR3063487.bi.vcf is the vcf file for SNPs.  
   ERR3063487.sv.primer is the list of primers to detect stuructural varilations.  
   ERR3063487.bi.primer is the list of primers to detect SNPs.  
   Primer files are experimental.  
   The algorithm of detection primer sequences has been developed by my experience of PCR experiment.
   
- Our verify process counts reads containing polymorphic region.  
  Basically, counts are from target and reference.  
  If control is specified, counts from target and control will be listed.  
  For example,  
```
perl ped.pl target=ERR3063487,control=ERR3063486,ref=WBcel235  
```
or  
```
perl bidirectional.pl ERR3063487 ERR3063486 WBcel235 
```
  returns 'M' or 'H' marked SNPs and structural variations of ERR3063487 which are absent in ERR3063486.  
- To confirm the alignment for detected polymorphisms,  
  perl search.pl target chr position  
  For example,  
```
perl search.pl target=ERR3063487,chr=II,pos=948033  
```
Alignments will be selected by the search script.  
- if you want to run with computer cluster, 
```
perl qsub_bidirectional.pl ERR3063487 default WBcel235  
```
- qsub_bidirectional.pl will work with Torque installed with default settings.  
  If output format of qsub and qstat are customized, modification of the
  script will be required.  
  Run without arguments, help for script will be shown.  
- ERR3063487.sv is the list of structural variations.  
  ERR3063487.bi.snp is the list of SNPs.  
  ERR3063487.bi.vcf is the vcf file for SNPs.  
  Quality score in vcf file is fixed to 1000.  
  Because our system does not use aligner program, *e.g.* bwa, output of quality score is difficult.  
  Please check quarity of polymorphism with depth (DP) in vcf file.  
- I searched small size short reads suitable for demonstration of PED from SRA in NCBI.  
  I found data set of *Caenorhabditis elegans.*  
  https://www.ncbi.nlm.nih.gov/bioproject/PRJEB30822  
  I express special thanks for the release of short reads to SRA by ENS Lyon.  

## Instruction for *k*-mer method
- For example,  
```
perl ped.pl target=ERR3063487,control=ERR3063486,ref=WBcel235,method=kmer  
```
or  
```
perl kmer.pl ERR3063487 ERR3063486 WBcel235  
```
- if you want to run with computer cluster,  
```
perl qsub_kmer.pl ERR3063487 ERR3063486 WBcel235   
```
  ERR3063487 specific SNPs against ERR3063486 will be detected.  
- if you want to detect SNPs against reference genome,  
```
perl ped.pl target=ERR3063487,ref=WBcel235,method=kmer  
```
or  
```
perl kmer.pl ERR3063487 default WBcel235  
```
- ERR3063487.kmer.snp is the list of SNPs.  
  ERR3063487.kmer.vcf is the vcf file for SNPs.   
  ERR3063487.kmer.primer is the list of primer sequence to detect SNPs.   
  Primer files are experimental.  
  The algorithm of detection primer sequences has been developed by my experience of PCR experiment.  
  Quality score in vcf file is fixed to 1000.  
  Because our system does not use aligner program, *e.g.* bwa, output of quality score is difficult.  
  Please check quarity of polymorphism with depth (DP) in vcf file.  
- The *k*-mer method is able to detect polymorphisms by the direct comparison between two short read data.  
  If you want to SNP detection between target and control without reference data,  
  run script without reference specification.  
  For example,  
```
perl ped.pl target=ERR3063487,control=ERR3063486,method=kmer  
```
or  
```
perl kmer.pl ERR3063487 ERR3063486  
```
or  
```
perl qsub_kmer.pl ERR3063487 ERR3063486  
```
- ERR3063487.kmer is the list of polymorphic edge.  
  SNPs tagged with first 19-mer will be used for the genetic analysis,  
  such as segregation analysis.  
  The 19-mer can be used as the identifire (*i.e.* name) for analysis.  

## Examples of result
A part of SNP list of the bidirectional method is  
```
1       3189273 T       C       50      0       15      9       H
1       3189333 T       C       50      0       14      0       _
1       3189345 A       G       50      0       13      17      H
1       3189429 G       A       50      0       15      0       _
1       3189498 G       A       50      0       0       39      M
1       3189503 T       G       50      0       2       1       _
1       3189527 T       G       50      0       10      0       _
1       3189609 T       A       50      0       41      1       _
1       3189704 C       G       50      0       0       11      M
1       3189741 A       G       50      0       0       23      M
1       3189819 A       G       50      0       0       4       _
1       3189833 C       T       50      0       0       46      M

Column 1: Chromosome number
Column 2: Position of SNP
Column 3: Reference base at the SNP position
Column 4: Alternative base
Column 5: Number of detected reads with alternative base
Column 6: Number of reads in the control sort_uniq file with control type polymorphism
Column 7: Number of reads in the control sort_uniq file with target type polymorphism
Column 8: Number of reads in the target sort_uniq file with control type polymorphism
Column 9: Number of reads in the target sort_uniq file with target type polymorphism
Column 10: Genotype (M: homozygous, H: heterozygous, _: reference type)
Following columns will be appeared in the primer file.
Column 11: Left primer sequence 
Column 12: Right primer sequence
Column 13: Estimated size of the amplified fragment
Column 14: Upstream and downstream sequence around the SNP
Primer files are experimental.  
The algorithm of detection primer sequences has been developed by my experience of PCR experiment.  
```
- A part of the structural variation result is
```
1       902788  1       902774  f       deletion        -1      50      0       10      9       H       AAAAAAAAAAAAAA
1       907869  1       907835  f       insertion       32      50      0       19      11      H       C
1       911222  1       911221  f       deletion        -1      50      0       21      15      H       T
1       923312  1       923312  f       deletion        -1      50      0       0       34      M       _
1       931147  1       931131  f       insertion       4       50      0       7       19      _       CCCTCCCTCCC
1       932618  1       932617  f       deletion        -4      50      0       1       32      M       TTTC

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
Column 12: Genotype (M: homozygous, H: heterozygous, _: reference type)
Column 13: Sequence between junctions
Following columns will be appeared in the primer file.
Column 14: Left primer sequence 
Column 15: Right primer sequence
Column 16: Estimated size of the amplified fragment
Column 17: Upstream and downstream sequences around the junction
Primer files are experimental.  
The algorithm of detection primer sequences has been developed by my experience of PCR experiment.
```

- A part of SNP result by the *k*-mer method is
```
X       14544549        ACAGTTGTATTTTTCAATT     A       A       G       r       0       0       0       13      0       27      0       0       13      0       0       30      M
X       14544549        TACACCACTGTAAGTCAAC     A       A       G       f       13      0       0       0       0       0       27      0       13      0       0       30      M
X       14592687        AAAAATTTGGATTTTTGGA     G       G       GT      r       0       15      0       0       20      20      0       1       5       0       10      0       _
X       14592707        AAAAATTTGGATTTTTGGA     G       G       AG      f       0       0       15      0       20      0       20      0       5       0       10      0       _
X       14615799        ACCCCTATATATAGTGTTT     T       T       AT      r       25      0       0       0       16      0       0       12      26      0       18      0       _
X       14615819        ACCCCTATATATAGTGTTT     T       T       AT      f       0       0       0       25      12      0       0       16      26      0       23      0       _
X       14624654        GTTGTGCTTTATTTATTTG     A       A       AT      r       0       0       0       19      15      0       0       20      19      0       19      0       _
X       14624674        GTTGTGCTTTATTTATTTG     A       A       AG      f       18      0       0       0       18      0       16      0       19      0       17      0       _
X       14632040        ATTTTTTTCACAAACAAGG     T       T       A       r       26      0       0       0       0       0       0       23      21      0       0       21      M
X       14632040        CTTAAATCGGAGAACAAAT     T       T       A       f       0       0       0       25      22      0       0       0       21      0       0       21      M
X       14682371        TCAGTTCAATCACGATAAA     A       A       AT      r       0       0       0       19      10      0       0       19      24      0       20      0       _

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
Column 20: Genotype (M: homozygous, H: heterozygous, _: reference type)
Following columns will be appeared in the primer file.
Column 21: Left primer sequence 
Column 22: Right primer sequence
Column 23: Estimated size of the amplified fragment
Column 24: Upstream and downstream sequence around the SNP```
Primer files are experimental.  
The algorithm of detection primer sequences has been developed by my experience of PCR experiment.
```
## Detection of polymorphisms between control and target
- SNPs between ERR3063486 (wild-type) and ERR3063487 (mutant)  
```
I       27950   A       T       31      0       0       17      M
I       892680  C       G       8       1       5       3       H
I       1196268 T       A       9       0       8       4       H
I       1380502 A       T       22      0       0       13      M
I       1728826 T       G       8       0       6       3       H
I       3203676 T       G       16      1       11      5       H
I       3407954 C       A       28      0       0       23      M
I       3656814 A       C       8       1       8       4       H
I       5001132 T       G       8       1       6       3       H
I       6324213 G       T       24      0       0       21      M
I       7249395 T       G       19      1       7       4       H
I       7263091 T       G       14      1       9       6       H
I       9136539 T       A       20      0       0       16      M
I       10137843        T       G       6       0       9       4       H
I       11168832        A       C       5       1       5       3       H
I       14097862        A       T       17      1       9       13      H
II      86891   A       G       7       0       8       4       H
II      271122  C       A       16      0       6       3       H
II      1179320 C       A       11      0       0       17      M
II      2500482 C       A       15      0       0       24      M
II      3552886 A       C       15      1       5       3       H
II      3624396 A       C       12      1       9       5       H
II      3648966 C       T       9       0       0       27      M
II      3771956 G       C       19      0       0       16      M
II      3935824 A       C       7       0       6       3       H
II      3979685 T       G       9       0       5       3       H
II      8284226 C       A       19      0       0       20      M
II      8447001 T       G       10      0       9       4       H
II      8553707 T       A       21      0       0       6       M
II      9410187 C       T       23      0       0       18      M
II      9937543 T       G       21      0       0       18      M
II      10629519        A       C       17      1       12      8       H
II      10685303        T       A       21      0       0       16      M
II      12732768        A       C       13      1       6       3       H
II      14096056        T       G       6       0       5       4       H
III     198231  T       A       16      0       0       18      M
III     3091906 T       G       19      0       9       4       H
III     3643248 G       C       6       0       6       3       H
III     4824486 C       G       33      0       0       23      M
III     7532164 T       G       8       1       6       3       H
III     9723566 G       A       25      0       0       23      M
III     11532166        C       T       15      0       1       18      M
III     13092063        C       T       13      0       6       6       H
III     13273147        A       T       15      1       6       4       H
III     13350315        A       C       14      0       7       4       H
IV      554740  A       C       11      1       9       4       H
IV      876681  A       G       20      0       7       4       H
IV      1159003 A       C       13      0       7       6       H
IV      1327294 G       A       26      0       0       31      M
IV      2417073 G       A       32      0       0       21      M
IV      2858610 C       T       21      0       0       14      M
IV      3931877 G       A       13      0       0       20      M
IV      4298928 A       C       11      1       9       4       H
IV      5200355 G       T       27      0       0       23      M
IV      6481216 C       G       17      0       0       8       M
IV      6796145 C       T       16      0       0       19      M
IV      6967218 G       A       25      0       0       25      M
IV      8053747 T       A       23      0       0       19      M
IV      8951120 T       G       10      1       13      6       H
IV      9709645 C       A       22      0       0       27      M
IV      10195034        T       G       14      1       8       4       H
IV      13672183        C       T       18      0       0       18      M
IV      14217812        A       C       9       0       6       3       H
IV      14760376        T       G       9       0       5       3       H
IV      16870604        A       C       12      1       8       4       H
V       7011    T       G       10      1       11      5       H
V       974819  A       C       6       0       9       4       H
V       3052728 A       C       9       1       10      5       H
V       3277678 G       A       22      0       1       19      M
V       3786240 T       A       23      0       0       17      M
V       4261370 T       G       14      0       9       4       H
V       5134172 C       T       19      0       17      10      H
V       7816318 A       C       20      1       8       4       H
V       9771516 A       T       25      0       18      12      H
V       10513400        A       C       11      0       9       4       H
V       11310419        G       T       15      0       0       18      M
V       15880652        T       C       20      1       6       3       H
V       19657843        T       A       25      0       0       24      M
V       19718778        G       A       11      0       0       11      M
V       19733914        A       C       6       0       6       3       H
V       19742410        T       G       6       1       5       3       H
V       20202209        T       G       5       0       5       3       H
V       20413571        T       G       13      1       6       3       H
X       2843588 A       C       13      0       9       4       H
X       4194330 T       C       23      0       0       22      M
X       4247242 T       G       7       1       6       3       H
X       5446454 C       T       21      0       0       7       M
X       6994561 A       T       29      0       0       29      M
X       7299494 A       G       7       0       6       3       H
X       7586849 T       G       8       1       8       4       H
X       10486673        A       T       31      0       0       23      M
X       12500491        T       G       7       1       6       3       H
X       14544549        A       G       13      0       0       30      M
X       14632040        T       A       21      0       0       21      M
X       14735899        T       G       8       0       5       3       H
X       15395882        A       G       17      0       13      6       H
X       15815152        T       G       13      1       5       3       H
X       16861146        C       T       9       0       0       11      M
X       16964164        T       G       14      0       5       3       H
```
- Structural variations between ERR3063486 (wild-type) and ERR3063487 (mutant)  
```
I       834776  I       834746  f       deletion        -12     6       0       7       4       H       TCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT
I       1251597 I       1251573 f       deletion        -2      5       0       8       4       H       TGTGTGTGTGTGTGTGTGTGTGTGT
I       1412142 I       1411952 f       insertion       102     12      1       6       3       H       CCCCCCGCTGACCCCAAACCAATATCCCGTCAAAAAACGAAAATTCATATTTTTCTTAATCTACAGTAATCCTACAGTGCCCCTACA
I       1560682 I       1560674 f       deletion        -1      17      0       0       19      M       TTTTTTTT
I       1560682 I       1561254 r       inversion       _       13      0       0       6       M       TTAAAGGTGGTGTGGTCGAATTTTTTTT
I       2160919 I       2160898 f       deletion        -1      9       0       8       4       H       TTTTTTTTCAAAAAAAAAAAA
I       2384261 I       2384244 f       insertion       1       14      1       8       4       H       CAAAAAAAAAAAAAA
I       2715230 III     12981863        r       translocation   _       5       1       5       3       H       CGTATTGCACAGCACATTTGACGCGCAAAAT
I       5081495 I       5081478 f       deletion        -2      7       1       6       4       H       TTTTTTTTTTTTTTTTTT
I       8028077 I       8028066 f       deletion        -1      6       0       6       3       H       TTTTTTTTTTT
I       8028078 I       8028065 f       insertion       1       6       0       6       3       H       TTTTTTTTTTT
I       9825014 I       9825007 f       insertion       1       23      0       8       13      H       AAAAA
I       10622455        IV      9192034 f       translocation   _       5       1       6       3       H       GTTCAAATAAAAATATTTTTTT
I       10887954        I       10887958        f       deletion        -6      11      0       1       10      M       A
I       11005509        I       11005489        f       deletion        -1      5       1       0       5       M       AAATTTTTTTTTTTTTTTTT
I       12856424        I       12856415        f       deletion        -1      14      0       8       5       H       TTTTTTTTT
I       13734806        I       13734788        f       deletion        -1      5       0       6       4       H       TTTTTTTTTTTTTTTTTT
II      191614  II      191601  f       insertion       1       11      0       9       4       H       AAAAAAAAAAA
II      948033  II      948099  f       deletion        -69     20      0       1       15      M       AT
II      1777672 II      1770125 f       insertion       7506    8       1       9       4       H       ATGGTGAGTAGCCGGTAATTTCATAGTTATTGAAATTTGA
II      2361947 II      2361931 f       deletion        -1      10      1       0       15      M       TTTTTCTTTTTTTTTT
II      4463297 V       1000739 r       translocation   _       6       0       6       4       H       TTTCGATTTTCCAGAAAATCAAAAAAAAA
II      4895056 V       18778607        r       translocation   _       5       0       5       3       H       TTCTACGTTTTGCAATGTGTTTTTT
II      5473562 II      3675168 r       inversion       _       8       1       5       3       H       TTTTACTCAGTTATGTTTTTTTT
II      12746320        III     4520340 f       translocation   _       6       1       5       3       H       TGTAAAATTGTTTTTTTTT
II      13112130        V       20740040        f       translocation   _       6       0       7       4       H       AAAAAAAAACGCATGCATTTTTCG
II      13327324        II      13327697        r       inversion       _       6       1       6       5       H       TTTTGACACTTTTTAGTAATAAATGCAAAAAAAATCAACAAAAATAGACTAAACATTGTAAAAACTGTAAAAACTAAGAGAAAAAAT
III     1663181 III     1663357 f       deletion        -209    6       1       7       5       H       TTTTTTCCAGAAATTAATATTTCTAGAAAAAT
III     2300741 X       15839886        r       translocation   _       19      0       11      6       H       TTAAAGGTGGAGTAGCGCCAGTGGGAAAATTGCTTTAAAACATGCCTATGGTACCACAATGACCAAATATCAT
III     2520715 X       97782   f       translocation   _       7       1       6       5       H       TATTTTTTCGCCATTTTTTTT
III     2985966 IV      876821  f       translocation   _       5       1       6       3       H       AAAAAAATTTTTTTTTT
III     3566956 III     3566962 f       deletion        -48     8       1       6       4       H       TCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT
III     6231151 III     6231140 f       deletion        -1      14      0       6       3       H       AAAAAAAAAAA
III     10402397        III     10402387        f       deletion        -2      5       1       7       4       H       TTTTTTTTTTT
III     11280814        X       2272016 r       translocation   _       12      0       6       3       H       TTTTTTTCAAAAAAAAAAAAA
III     12543896        III     12543907        f       deletion        -62     8       1       9       4       H       AAATTTCCGGAAAACATGCAAATTGCCAGAATTGAAAATTTCCGGCAAAT
III     13419473        III     13419443        f       insertion       6       5       1       10      5       H       TGTGTGTGTGTGTGTGTGTGTGT
IV      1786345 IV      1786337 f       deletion        -1      15      0       0       12      M       AAAAAAAA
IV      2112309 IV      2112287 f       deletion        -1      7       1       5       3       H       TTTTTTGTTTTTTTTTTTTTTT
IV      2314588 IV      2314573 f       deletion        -1      10      1       8       4       H       TTTTTTTTTTTTTTT
IV      2445289 IV      2445276 f       deletion        -1      15      0       5       3       H       TTTTTTTTTTTTT
IV      3192017 IV      3192001 f       deletion        -1      10      1       8       5       H       AAAAAAAAAAAAAAAA
IV      3192018 IV      3192000 f       insertion       1       10      0       7       4       H       AAAAAAAAAAAAAAAA
IV      3336297 IV      3336328 f       deletion        -31     21      0       0       23      M       A
IV      3867139 IV      3867116 f       deletion        -2      11      1       6       3       H       ATATATATATATATATATATATAT
IV      4399486 IV      4399441 f       insertion       2       8       0       10      6       H       GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA
IV      4489131 IV      4489122 f       deletion        -1      16      0       1       17      M       TTTTTTTTT
V       1027708 V       1027700 f       deletion        -39     9       1       10      5       H       GCCTATGGCCTACGCCTATGGCCTACGCCTATGGCCTACGCCTATG
V       1989616 V       1989642 f       deletion        -3      6       1       6       3       H       GATAAAAACTACTTGGATAAATGA
V       9026732 V       9026728 f       deletion        -3      6       1       7       4       H       ATTATT
V       11884926        II      3151004 r       translocation   _       5       1       6       5       H       GTGCCGAGTGCCGATCGGCACAATGTG
V       18672436        IV      2681633 r       translocation   _       7       1       5       3       H       GGGAAAATTGCTTTAAAACATGCCTATGGTACTACAA
V       18890234        V       18892184        r       inversion       _       6       1       6       3       H       TTTTTATTGAAAACTAGTATAAAAATATA
V       20123760        V       20123893        f       deletion        -150    12      0       7       4       H       GGGGTTCGAACCCCGG
X       99547   X       99523   f       deletion        -1      8       1       9       5       H       AAAAAAAATTTTTTTTTTTTTTTT
X       438457  X       438446  f       insertion       1       12      0       5       5       H       TTTTTTTTT
X       562164  X       562155  f       insertion       1       22      0       0       22      M       AAAAAAA
X       1522918 X       1522902 f       deletion        -1      6       1       5       3       H       AAAAAAAAAAAAAAAA
X       2498483 X       2498466 f       deletion        -1      6       1       5       3       H       TTTTTTTTTTTTTTTTT
X       3823840 X       3823828 f       deletion        -1      10      0       8       4       H       AAAAAAAAAAAA
X       5885390 X       5885377 f       deletion        -1      11      0       10      5       H       TTTTTTTTTTTTT
X       7312229 X       11435923        r       inversion       _       6       0       5       4       H       TATTCACCCCGTTCGACTGTGCAATGGGTTTAATCTATTCACTTTGTAAATCAAAGAATCGACGACCGCCTCCTGAA
X       10023790        III     7850798 r       translocation   _       5       0       6       3       H       ATATCAAAATTTCATTTTTTTT
X       14258766        III     303520  f       translocation   _       6       0       9       5       H       TCACAAAATTCTTTGGCCGCCCCAAGTGTCCTAACTCGAAG
```

## Author
MIYAO Akio, Ph.D. miyao@affrc.go.jp  
Institute of Crop Science / National Agriculture and Food Research Organization  
2-1-2, Kannondai, Tsukuba, Ibaraki 305-8518, Japan  

## Version
Version 1.3 Update for search.pl for comfirmation of alignment. Improuvement of making sort_uniq data.  
Version 1.2 sort_uniq files are divided to 64 subfiles by first three nucleotide sequence. Remake of reference data is required.   
Version 1.1 sort_uniq files are compressed by gzip. Requirement of disk space is reduced, but requires more CPU time.  
Version 1.0 Original version for PED paper.  

## Citing PED
Cite this article as: Polymorphic edge detection (PED): two efficient methods of polymorphism detection from next-generation sequencing data  
Akio Miyao, Jianyu Song Kiyomiya, Keiko Iida, Koji Doi, Hiroshi Yasue  
BMC Bioinformatics. 2019 20(1):362.  
URL https://doi.org/10.1186/s12859-019-2955-6  
PDF https://rdcu.be/bH7e8  

## License
NARO NON-COMMERCIAL LICENSE AGREEMENT Version 1.0  

This license is for 'Non-Commercial' use of software for polymorphic edge detection (PED)

1. Scientific use of PED is permitted free of charge.
1. Modification of PED is only permitted to the person of downloaded and his/her colleagues.
1. The National Agriculture and Food Research Organization (hereinafter referred to as NARO) does not guarantee that defects, errors or malfunction will not occur with respect to PED.
1. NARO shall not be responsible or liable for any damage or loss caused or be alleged to be caused, directly or indirectly, by the download and use of PED.
1. NARO shall not be obligated to correct or repair the program regardless of the extent, even if there are any defects of malfunctions in PED.
1. The copyright and all other rights of PED belong to NARO.
1. Selling, renting, re-use of license, or use for business purposes etc. of PED shall not be allowed. For commercial use, license of commercial use is required. Inquiries for such commercial license are directed to ped_request@ml.affrc.go.jp.
1. The PED may be changed, or the distribution maybe canceled without advance notification.
1. In case the result obtained using PED in used for publication in academic journals *etc.,* please refer the publication of PED and/or acknowledge the use of PED in the publication. 

Copyright (C) 2017 [National Agriculture and Food Research Organization](https://www.naro.affrc.go.jp/english/index.html). All rights reserved.  
