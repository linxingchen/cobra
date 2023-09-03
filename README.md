# COBRA
COBRA (Contig Overlap Based Re-Assembly) is a bioinformatics tool to get higher quality viral genomes assembled from metagenomes of short paired-end reads. COBRA was written in Python. COBRA has so far only been tested on assembled contigs from metaSPAdes, IDBA_UD, and MEGAHIT.

## Introduction

**1. Why metagenomic contigs are fragmented?**

The genomes assembled from short paired-end reads based metagenomes are usually fragmented due to (1) intra-genome repeats, (2) inter-genome shared region, and (3) within-population variations, as the widely utilized assemblers based on de Bruijn graphs, e.g., metaSPAdes, IDBA_UD and MEGAHIT, tend to have a breaking point when multiple paths are available instead of making risky extension (see example in **Figure 1**). 

![image](https://user-images.githubusercontent.com/46725273/111676563-8a21f180-87db-11eb-9b8c-4c63fb993936.png)

Figure 1. Example of how assemblers break in assembly when within-population occurs.

**2. Contigs may be joined with expected end overlap.**

According to the principles of the abovementioned assemblers, the broken contigs have an end overlap with determined length, that is the max-kmer (maxK hereafter) used in de nono assembly for metaSPAdes and MEGAHIT, and the maxK-1 for IDBA_UD, which we termed as "expected overlap length" (Figures 1 and 2). 

* Note: as COBRA will use the information provided by paired-end reads, thus only those samples sequenced by paired-end technology should work.

![image](https://user-images.githubusercontent.com/46725273/111677281-4c719880-87dc-11eb-85a9-a62906f4e10b.png)

Figure 2. The "expected overlap length" has been documented in manual genome curation, see [Chen et al. 2020](https://genome.cshlp.org/content/30/3/315.short) for details.

##
## How COBRA works

COBRA determines the "expected overlap length" (both the forward direction and reverse complement direction) for all the contigs from an assembly, then looks for the valid joining path for each query that users provide (should be a fraction of the whole assembly) based on a list of features including contig coverage, contig overlap relationships, and contig continuity (based on paired-end reads mapping) (Figure 3).

Note that scaffolds (for example, metaSPAdes assembly) could be used as input for COBRA extension, however, we suggest not using scaffolds from IDBA_UD as the potential errors in the scaffolding step (see [Chen et al. 2020](https://genome.cshlp.org/content/30/3/315.short) for details). Thus, for IDBA_UD and MEGAHIT assembly, the contigs should be used. 

Given that COBRA has only tested for contigs/scaffolds from IDBA_UD, metaSPAdes and MEGAHIT, it will be risky to use it on contigs/scaffolds from any other assemblers.

![Figure 1](https://github.com/linxingchen/cobra.github.io/assets/46725273/b4dfe26a-794a-4384-9c15-12762b5e7e10)

Figure 3. The workflow of COBRA.

##
## Installation
COBRA is a Python script (tested for version 3.7 or higher) that uses a list of frequently used Python packages including:
```
Bio
Bio.Seq
collections
argparse
math
pysam
time
```

Once these packages are available, the only third-party software that the user should have is [BLASTn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

##
## Input files
(1) COBRA needs four files as inputs, i.e., 

* ```-f/--fasta```: A fasta format file containing all the contigs from a single assembly, note that IDBA_UD and MEGAHIT usually save contigs with a minimum length of 200 bp.

* ```-c/--coverage```: a two columns (separated by tab) file of the sequencing coverage of all contigs in the ```-f/--fasta``` file, example below:

```contig-140_0    25.552
contig-140_1    42.1388
contig-140_2    14.6023
contig-140_3    15.4817
contig-140_4    41.2746
...
```

* ```-q/--query```: the query contigs that the user wants COBRA to extend, could be provided in a fasta format file, or a one-column text file with the names of the query contigs. Please make sure the names are exactly the same format as in the ```-f/--fasta``` file, otherwise, COBRA may have problems extending them.
  
* ```-m/--mapping```: the paired-end reads mapping file of all contigs in the ```-f/--fasta``` file, could be sam or bam file.

(2) and three parameters
* ```-a/--assembler```: the name of the de novo assembler used, currently only 'idba' (for IDBA_UD), 'metaspades' (for metaSPAdes), and 'megahit' (for MEGAHIT).
* ```-maxk/--maxk```: the largest kmer used in de novo assembly.
* ```-mink/--mink```: the smallest kmer used in de novo assembly.


(3) Optional flags
* ```-lm/--linkage_mismatch```: the number of read mapping mismatches allowed when determining if two contigs were spanned by paired reads.
* ```-o/--output```: the name of the output folder, otherwise it will be "{```-q/--query```}.COBRA" if not provided.
* ```-t/--threads```: the number of threads used for BLASTn search.

##
## How to obtain the mapping file

The mapping file could be obtained with tools like Bowtie2 and BBMap, please refer to the manual descriptions for details of the tools. Below is the general way to get the sorted sam/bam file, you thus need to be available to samtools (which could be get here - https://github.com/samtools/samtools).

For example, 

* ```contig file = "contigs.fasta"```

* ```first read file = "R1.fastq.gz"```

* ```second read file = "R2.fastq.gz"```

(1) with Bowtie2 (https://github.com/BenLangmead/bowtie2)

* ```bowtie2-build contigs.fasta contigs.fasta```

* ```bowtie2 -p 16 -x contigs.fasta -1 R1.fastq.gz -2 R2.fastq.gz -S output.sam && samtools view -bS output.sam | samtools sort -o sorted_output.bam -```


(2) with BBMap (https://github.com/BioInfoTools/BBMap)

* ```bbmap.sh ref=contigs.fasta in1=R1.fastq.gz in2=R2.fastq.gz threads=16 out=output.sam``` (good)
* ```samtools view -bS output.sam > output.bam```
* ```samtools sort -o sorted_output.bam output.bam```


##  How to obtain the coverage file

(1) with jgi_summarize_bam_contig_depths

Once the sorted sam or bam file is ready, the tool of "jgi_summarize_bam_contig_depths" from MetaBAT (https://bitbucket.org/berkeleylab/metabat/src/master/), or could be used to obtain the coverage file, the resulting profile should be transferred to get a two-column file divided by tab.

* ```jgi_summarize_bam_contig_depths --outputDepth original.coverage.txt *sam``` 

* ```jgi_summarize_bam_contig_depths --outputDepth original.coverage.txt *bam```

The output file from jgi_summarize_bam_contig_depths could be converted to a two-column file divided by tab using the script provided in this study (coverage.transfer.py).

* ```python coverage.transfer.py -i original.coverage.txt -o coverage.txt```


(2) CoverM 

CoverM is a fast DNA read coverage and relative abundance calculator focused on metagenomics applications. Usage could be found here (https://github.com/wwood/CoverM).


(3) pyCoverM 

pyCoverM is a simple Python interface to CoverM's fast coverage estimation functions, which could be found here (https://github.com/apcamargo/pycoverm).


##
## How to run
(1) The users can only specify the required parameters:

```
python cobra.py -f input.fasta -q query.fasta -c coverage.txt -m mapping.sam -a idba -mink 20 -maxk 140
```

(2) The users could also include the optional parameters like output name (-o), mismatch of mapped reads for linkage identification (-mm)

```
python cobra.py -f all.contigs.fasta -q query.fasta -o query.fasta.COBRA.out -c coverage.txt -m mapping.sam -a idba -mink 20 -maxk 140 -mm 2
```

```
python cobra.py -f all.contigs.fasta -q query.fasta -o query.fasta.COBRA.out -c coverage.txt -m mapping.sam -a metaspades -mink 21 -maxk 127 -mm 2
```

```
python cobra.py -f all.contigs.fasta -q query.fasta -o query.fasta.COBRA.out -c coverage.txt -m mapping.sam -a megahit -mink 21 -maxk 141 -mm 2
```

##
## Output files
Below is a general list of output files in the output folder:

```
COBRA_category_i_self_circular_queries_trimmed.fasta
COBRA_category_i_self_circular_queries_trimmed.fasta.summary.txt
COBRA_category_ii_extended_circular_unique (folder)
COBRA_category_ii_extended_circular_unique.fasta
COBRA_category_ii_extended_circular_unique.fasta.summary.txt
COBRA_category_ii_extended_circular_unique_joining_details.txt
COBRA_category_ii_extended_failed.fasta
COBRA_category_ii_extended_failed.fasta.summary.txt
COBRA_category_ii_extended_partial_unique (folder)
COBRA_category_ii_extended_partial_unique.fasta
COBRA_category_ii_extended_partial_unique.fasta.summary.txt
COBRA_category_ii_extended_partial_unique_joining_details.txt
COBRA_category_iii_orphan_end.fasta
COBRA_category_iii_orphan_end.fasta.summary.txt
COBRA_joining_status.txt
COBRA_joining_summary.txt
intermediate.files (folder)
log
debug.txt
contig.new.fa
```

For all the queries, COBRA assigns them to different categories based on their joining status (detailed in the ```COBRA_joining_status.txt``` file), i.e.,

* "self_circular" - the query contig itself is a circular genome.
* "extended_circular" - the query contig was joined and extended into a circular genome.
* "extended_partial" - the query contig was joined and extended but not to circular.
* "extended_failed" - the query contig was not able to be extended due to COBRA rules. 
* "orphan end" - neither end of a given contig shares "expected overlap length" with others.

For the joined and extended queries in each category, only the unique ones (```*.fasta```) will be saved for users' following analyses, and the sequence information (e.g., length, coverage, GC, num of Ns) is summarized in the ```*fasta.summary.txt``` files. For categories of "extended_circular", and "extended_partial", the joining details of each query are included in the corresponding folder and ```*joining_details.txt``` file, and summarized in the ```COBRA_joining_summary.txt``` file, an example shown below:

```
QuerySeqID      QuerySeqLen     TotRetSeqs      TotRetLen       AssembledLen    ExtendedLen     Status
contig-140_100  47501   3       50379   49962   2461    Extended_circular
contig-140_112  45060   3       62549   62132   17072   Extended_circular
contig-140_114  44829   2       45342   45064   235     Extended_circular
contig-140_160  40329   2       41018   40740   411     Extended_circular
contig-140_188  38386   5       48986   48291   9905    Extended_circular
...
```


* **log file:** The ```log``` file includes the content of each processing step, an example shown below:

```
1. INPUT INFORMATION
# Assembler: IDBA_UD
# Min-kmer: 20
# Max-kmer: 140
# Overlap length: 139 bp
# Read mapping max mismatches for contig linkage: 2
# Query contigs: file-path
# Whole contig set: file-path
# Mapping file: file-path
# Coverage file: file-path
# Output folder: file-path

2. PROCESSING STEPS
[01/22] [2023/05/04 12:48:15] Reading contigs and getting the contig end sequences. A total of 311739 contigs were imported.
[02/22] [2023/05/04 12:48:33] Getting shared contig ends.
[03/22] [2023/05/04 12:48:40] Writing contig end joining pairs.
[04/22] [2023/05/04 12:48:41] Getting contig coverage information.
[05/22] [2023/05/04 12:48:42] Getting query contig list. A total of 2304 query contigs were imported.
[06/22] [2023/05/04 12:48:43] Getting contig linkage based on sam/bam. Be patient, this may take long.
[07/22] [2023/05/04 13:00:01] Detecting self_circular contigs.
[08/22] [2023/05/04 13:00:37] Detecting joins of contigs. 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, 100% finished.
[09/22] [2023/05/04 14:29:20] Saving potential joining paths.
[10/22] [2023/05/04 14:29:23] Checking for invalid joining: sharing queries.
[11/22] [2023/05/04 14:29:28] Getting initial joining status of each query contig.
[12/22] [2023/05/04 14:29:28] Getting final joining status of each query contig.
[13/22] [2023/05/04 14:29:28] Getting the joining order of contigs.
[14/22] [2023/05/04 14:29:28] Getting retrieved contigs.
[15/22] [2023/05/04 14:29:32] Saving joined seqeuences.
[16/22] [2023/05/04 14:29:35] Checking for invalid joining using BLASTn: close strains.
[17/22] [2023/05/04 14:30:11] Saving unique sequences of "Extended_circular" and "Extended_partial" for joining checking.
[18/22] [2023/05/04 14:30:12] Getting the joining details of unique "Extended_circular" and "Extended_partial" query contigs.
[19/22] [2023/05/04 14:30:12] Saving joining summary of "Extended_circular" and "Extended_partial" query contigs.
[20/22] [2023/05/04 14:30:15] Saving joining status of all query contigs.
[21/22] [2023/05/04 14:30:15] Saving self_circular contigs.
[22/22] [2023/05/04 14:30:16] Saving the new fasta file.

======================================================================================================================================================
3. RESULTS SUMMARY
# Total queries: 2304
# Category i   - Self_circular: 74
# Category ii  - Extended_circular: 120 (Unique: 82)
# Category ii  - Extended_partial: 1088 (Unique: 889)
# Category ii  - Extended_failed (due to COBRA rules): 245
# Category iii - Orphan end: 777
# Check "COBRA_joining_status.txt" for joining status of each query.
# Check "COBRA_joining_summary.txt" for joining details of "Extended_circular" and "Extended_partial" queries.
======================================================================================================================================================
```

##
## Citation
The preprint is available at bioRxiv (doi: https://doi.org/10.1101/2023.05.30.542503). Please cite if you find it helpful for your analyses.

##
## Questions
Please do not hesitate to contact me for any questions you have when using COBRA via my email (linkingchan@gmail.com), or Twitter (https://twitter.com/ChenLinxing).
