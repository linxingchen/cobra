# COBRA
COBRA (Contig Overlap Based Re-Assembly) is a bioinformatics tool to get higher quality virus genomes assembled from short-read metagenomes. Which was written in python.

## Introduction
* The genomes assembled from short-reads sequenced metagenomes are usually fragmented due to intra-genome repeats and within-population variations (or subpopulation diversity, or local diversity), as the widely used assemblers based on de Bruijn graphs, e.g., metaSPAdes, IDBA_UD and MEGAHIT, tend to have a breaking point when multiple paths are available instead of making risky extension (see example in **Figure 1**). 

![image](https://user-images.githubusercontent.com/46725273/111676563-8a21f180-87db-11eb-9b8c-4c63fb993936.png)

**Figure 1. Example of how assemblers break in assembly when within-population occurs.**

* According to the principles of the abovementioned assemblers, the broken contigs have an end overlap with determined length, that is the *max-kmer* used in de nono assembly for metaSPAdes and MEGAHIT, and the *max-kmer* - 1 for IDBA_UD, which we termed as ```expected overlap length (EOL)``` (**Figures 1 and 2**). 

*Note: as COBRA will use information provided by paired-end reads, thus only those samples sequenced by paired-end technology should work.*

![image](https://user-images.githubusercontent.com/46725273/111677281-4c719880-87dc-11eb-85a9-a62906f4e10b.png)

**Figure 2. The EOL have been documented in manual genome curation, see [Chen et al. 2020. Genome Research](https://genome.cshlp.org/content/30/3/315.short) for details.**

##
## How COBRA works
* COBRA determines the EOL (both the forward direction and reverse complement direction) for all the contigs from an assembly, then looks for the valid joining path for each query that users provide (should be a fraction of the whole assembly) based on a list of features including contig coverage, contig overlap relationships, and contig continuity (based on paired end reads mapping) (**Figrue 3**).

![Figure 1](https://github.com/linxingchen/cobra.github.io/assets/46725273/a5148ae5-50ee-4b3c-855a-acff10311a18)

**Figure 3. The workfolw of COBRA.**

##
## Installation
COBRA is a python script (tested for version 3.7 or higher) that uses a list of frequently used python packages including:
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
COBRA needs four files as inputs, i.e., 

* ```all.contigs.fasta``` - the whole contig set from the assembly, note that IDBA_UD and MEGAHIT usually save contigs with a minimun length of 200 bp.

* ```coverage.txt``` - a two columns (tab) file of the coverage of all contigs, example below:

```contig-140_0    25.552
contig-140_1    42.1388
contig-140_2    14.6023
contig-140_3    15.4817
contig-140_4    41.2746
...
```

* ```queries.fasta``` - the fasta file containing all the query contigs for joining path search.

* ```mapping.sam (or mapping.bam)``` - the paired-end reads mapping file of all contigs.


Two additional flags should be provided for COBRA to determine the EOL, i.e.,

* ```assembler``` - currently only 'idba' (for IDBA_UD), 'metaspades' (for metaSPAdes), and 'megahit' (for MEGAHIT).

* ```kmer``` - the largest kmer used in de novo assembly.


Optional flag
* ```mismatch``` - the number of read mapping mismatches allowed when determining if two contigs were spanned by paired reads.

##
## How to run

The output file will be ```{query}.COBRA``` if not specified via the ```-o``` flag.

```
COBRA.py -f all.contigs.fasta -q queries.fasta -o queries.fasta.COBRA.out -c coverage.txt -m mapping.sam -a idba -k 140 -mm 2
```

```
COBRA.py -f all.contigs.fasta -q queries.fasta -o queries.fasta.COBRA.out -c coverage.txt -m mapping.sam -a metaspades -k 127 -mm 2
```

```
COBRA.py -f all.contigs.fasta -q queries.fasta -o queries.fasta.COBRA.out -c coverage.txt -m mapping.sam -a megahit -k 141 -mm 2
```

##
## Output files
Below is a general list of output files in the ```queries.fasta.COBRA.out``` folder:

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
COBRA_category_iii_DNA_break.fasta
COBRA_category_iii_DNA_break.fasta.summary.txt
COBRA_joining_status.txt
COBRA_joining_summary.txt
intermediate.files (folder)
log
```

For all the quries, COBRA assigns them to different categories based on their joining status (detailed in the ```COBRA_joining_status.txt``` file), i.e.,

* "self_circular" - the query contig itself is a circular genome.
* "extended_circular" - the query contig was joined and exteneded into a circular genome.
* "extended_partial" - the query contig was joined and extended but not to circular.
* "extended_failed" - the query contig was not able to be extended due to COBRA rules. 
* "DNA break" - neither end of a given contig share EOL with others.

For the joined and extended queries in each category, only the unique ones (```*.fasta```) will be saved for users' following analyses, and the sequence information (e.g., length, coverage, GC, num of Ns) is summarized in the ```*fasta.summary.txt``` files. For categories of "extended_circular", and "extended_partial", the joining details of each query are included in the corresponding folder and ```*joining_details.txt``` file, and summarized in the ```COBRA_joining_summary.txt``` file, example shown below:

```
QuerySeqID      QuerySeqLen     TotRetSeqs      TotRetLen       AssembledLen    ExtendedLen     Status
contig-140_100  47501   3       50379   49962   2461    Extended_circular
contig-140_112  45060   3       62549   62132   17072   Extended_circular
contig-140_114  44829   2       45342   45064   235     Extended_circular
contig-140_160  40329   2       41018   40740   411     Extended_circular
contig-140_188  38386   5       48986   48291   9905    Extended_circular
...
```


* **log file:** The ```log``` file includes the content of each processing step, example shown below:

```
************************************* COBRA analyses for contig.with.97.identity.10k.alignment.to.polished.fasta *************************************

# Key parameters:
# Assembler: IDBA_UD
# Max-kmer: 140
# Overlap length: 139 bp
# Mismatch: 2

[01/20] [Sat Oct 16 13:52:29 2021] Reading contigs and getting contig ends ... A total of 311739 contigs were imported.
[02/20] [Sat Oct 16 13:52:37 2021] Getting shared contig ends ...
[03/20] [Sat Oct 16 13:52:43 2021] Writing contig end joining pairs ...
[04/20] [Sat Oct 16 13:52:43 2021] Getting contig coverage information ...
[05/20] [Sat Oct 16 13:52:43 2021] Getting query contig list ... 551 query contigs were imported ...
[06/20] [Sat Oct 16 13:52:43 2021] Getting contig linkage based on sam/bam ... Be patient, this may take long ...
[07/20] [Sat Oct 16 13:54:49 2021] Walking joins... 5%, 10%, 15%, 20%, 25%, 30%, 35%, 40%, 45%, 50%, 55%, 60%, 65%, 70%, 75%, 80%, 85%, 90%, 95%, 100% finished ...
[08/20] [Sat Oct 16 13:57:32 2021] Saving potential joining paths ...
[09/20] [Sat Oct 16 13:57:32 2021] Saving contig DNA break information ...
[10/20] [Sat Oct 16 13:57:32 2021] Checking for invalid joining - sharing queries ...
[11/20] [Sat Oct 16 13:57:32 2021] Getting the joining order of contigs ...
[12/20] [Sat Oct 16 13:57:32 2021] Getting retrieved contigs ...
[13/20] [Sat Oct 16 13:57:33 2021] Saving joined seqeuences ...
[14/20] [Sat Oct 16 13:57:33 2021] Getting initial joining status of each query contig ...
[15/20] [Sat Oct 16 13:57:33 2021] Getting final joining status of each query contig ...
[16/20] [Sat Oct 16 13:57:33 2021] Checking for invalid joining using BLASTn ...
[17/20] [Sat Oct 16 13:58:31 2021] Saving joining summary of retrieved contigs ...
[18/20] [Sat Oct 16 13:58:31 2021] Saving unique sequences of "Extended_circular" and "Extended_partial" ...
[19/20] [Sat Oct 16 13:59:42 2021] Saving joining details of unique "Extended_circular" and "Extended_partial" queries ...
[20/20] [Sat Oct 16 13:59:43 2021] Saving self_circular contigs ...
```

The ```log``` file also gives a summary of the joining status of all queries, example shown below:

```
======================================================================================================================================================
Final summary
Total queries: 551
Category i - Self_circular: 24
Category ii - Extended_circular: 53 (Unique: 37)
Category ii - Extended_partial: 289 (Unique: 239)
Category ii - Failed due to COBRA rules: 67
Category iii - Failed due to DNA break: 118
======================================================================================================================================================
```

##
## Citation
The manuscript is in preparation (Chen et al., in prep).
