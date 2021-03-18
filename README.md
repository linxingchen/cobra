# COBRA
COBRA (Contig Overlap Based Re-Assembly) is a bioinformatics tool to get higher quality virus genomes assembled from short-read metagenomes.

## Introduction
Virus genomes assembled from short-reads sequenced metagenomes are usually fragmented due to intra-genome repeats and within-population variations (or subpopulation diversity, or local diversity), as the widely used assemblers based on de Bruijn graphs, e.g., metaSPAdes, IDBA_UD and MEGAHIT, tend to have a breaking point when multiple paths are available instead of making risky extension. According to the principles of the abovementioned assemblers, the broken contigs have an end overlap with determined length, that is the *max-kmer* used in de nono assembly for metaSPAdes and MEGAHIT, and the *max-kmer* - 1 for IDBA_UD, which we termed as ```expected overlap length (EOL)``` (**Figure 1**). COBRA determines the EOL (both the forward direction and reverse complement direction) for all the contigs from an assembly, then looks for the valid joining path for each query that users provide (should be a fraction of the whole assembly) based on a list of features including contig coverage, contig overlap relationships, and contig continuity (based on paired end reads mapping) (**Figrue 2**).

*Note: as COBRA will use information provided by paired-end reads, thus only those samples sequenced by paired-end technology should work.*

<img src="https://user-images.githubusercontent.com/46725273/111421478-49768b00-86aa-11eb-8bea-9d4aa060a5e0.png">

**Figure 1. The EOL have been documented in manual genome curation, see [Chen et al. 2020. Genome Research](https://genome.cshlp.org/content/30/3/315.short) for details.**

## How COBRA works

![image](https://user-images.githubusercontent.com/46725273/111675243-25b26280-87da-11eb-9f28-60c8625e48c7.png)

**Figure 2. The workfolw of COBRA.**


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
*Note that metaSPAdes contigs have coverage information (not absolute coverage though) in their headers, which will be used by COBRA if no ```coverage.txt``` file is provided.*

* ```queries.fasta``` - the fasta file containing all the query contigs for joining path search.

* ```mapping.sam (or mapping.bam)``` - the paired-end reads mapping file of all contigs.

Two additional flags should be provided for COBRA to determine the EOL, i.e.,

* ```assembler``` - currently only 'idba' (for IDBA_UD), 'metaspades' (for metaSPAdes), and 'megahit' (for MEGAHIT).

* ```kmer``` - the largest kmer used in de novo assembly.


## How to run

The output file will be ```queries.fasta.COBRA``` if not specified via the ```-o``` flag.

```
COBRA.py -f all.contigs.fasta -c coverage.txt -q queries.fasta -m mapping.sam -a idba -k 100
```

If assembled with metaSPAdes, the coverage information in the contig headers could be used when no additional coverage file is provided.

```
COBRA.py -f all.contigs.fasta -q queries.fasta -m mapping.sam -a metaspades -k 127
```

## Output files
Below is a general list of output files in the ```queries.fasta.COBRA``` folder:

```
COBRA_extended_10k_fewer
COBRA_extended_10k_fewer.fasta
COBRA_extended_10k_fewer.fasta.summary.txt
COBRA_extended_10k_fewer_joining_details.txt
COBRA_extended_10k_or_more_unique
COBRA_extended_10k_or_more_unique.fasta
COBRA_extended_10k_or_more_unique.fasta.summary.txt
COBRA_extended_10k_or_more_unique_joining_details.txt
COBRA_extended_circular_dominant_unique
COBRA_extended_circular_dominant_unique.fasta
COBRA_extended_circular_dominant_unique.fasta.summary.txt
COBRA_extended_circular_dominant_unique_joining_details.txt
COBRA_extended_circular_rare_unique
COBRA_extended_circular_rare_unique.fasta
COBRA_extended_circular_rare_unique.fasta.summary.txt
COBRA_extended_failed.fasta
COBRA_extended_failed.fasta.summary.txt
COBRA_joining_status.txt
COBRA_joining_summary.txt
COBRA_self_circular_queries_trimmed.fasta
COBRA_self_circular_queries_trimmed.fasta.summary.txt
intermediate.files
log
```

For all the quries, COBRA assigns them to different categories based on their joining status (detailed in the ```COBRA_joining_status.txt``` file), i.e.,

* "self_circular" - the query contig itself is a circular genome.
* "extended_circular" - the query contig was joined and exteneded into a circular genome.
* "extended >= 10k" - the query contig was joined and extended at least 10 kbp.
* "extended < 10k" - the query contig was joined and extended fewer than 10 kbp.
* "extended_failed" - tthe query contig was not able to be extended due to different reasons (i.e., "COBRA rules", "DNA break or low abundance" and "Unexpected_assembly_break_or_short_piece_missing"). 

For the joined and extended queries in each category, only the unique ones (```*.fasta```) will be saved for users' following analyses, and the sequence information (e.g., length, coverage, GC, num of Ns) is summarized in the ```*fasta.summary.txt``` files. For categories of "extended_circular", "extended >= 10k" and "extended < 10k", the joining details of each query are included in the corresponding folder and ```*joining_details.txt``` file, and summarized in the ```COBRA_joining_summary.txt``` file, example shown below:

```
QuerySeqID      QuerySeqLen     TotRetSeqs      TotRetLen       AssembledLen    ExtendedLen     Status
contig-140_100  47501   3       50379   49962   2461    Extended_circular
contig-140_112  45060   3       62549   62132   17072   Extended_circular
contig-140_114  44829   2       45342   45064   235     Extended_circular
contig-140_160  40329   2       41018   40740   411     Extended_circular
contig-140_188  38386   5       48986   48291   9905    Extended_circular
...
```


For those queries resulting in the "extended_circular" category, if "equal_path" exists (**Figure 3**), COBRA will output both the ```dominant genome``` with the dominant path and the ```rare genome``` with the rare path (**Figure 4**).

![image](https://user-images.githubusercontent.com/46725273/111669292-e41eb900-87d3-11eb-8f23-b23b0b5cdb3c.png)

**Figure 3. The equal pathes originating from within-population diversity (or local diversity).**


![image](https://user-images.githubusercontent.com/46725273/111668676-390dff80-87d3-11eb-87e5-b16251f06b73.png)

**Figure 4. Information of COBRA "extended_circular" genomes (dominant and rare).**

* **log file:** The ```log``` file includes the content of each processing step, example shown below:

```
[Thu Mar 18 10:35:03 2021] [1/20] Reading contigs and getting contig ends ... A total of 311739 contigs were imported.
[Thu Mar 18 10:35:12 2021] [2/20] Getting shared contig ends ...
[Thu Mar 18 10:35:18 2021] [3/20] Writing contig end joining pairs ...
[Thu Mar 18 10:35:19 2021] [4/20] Getting contig coverage information ...
[Thu Mar 18 10:35:19 2021] [5/20] Getting query contig list ... 551 query contigs were imported ...
[Thu Mar 18 10:35:19 2021] [6/20] Getting contig linkage based on sam/bam ... Be patient, this may take long ... 
[Thu Mar 18 10:38:40 2021] [7/20] Walking joins. 1% finished ... 2% finished ... 3% finished ... 4% finished ... 5% finished ... 6% finished ... 7% finished ... 8% finished ... 9% finished ... 10% finished ... 11% finished ... 12% finished ... 13% finished ... 14% finished ... 15% finished ... 16% finished ... 17% finished ... 18% finished ... 19% finished ... 20% finished ... 21% finished ... 22% finished ... 23% finished ... 24% finished ... 25% finished ... 26% finished ... 27% finished ... 28% finished ... 29% finished ... 30% finished ... 31% finished ... 32% finished ... 33% finished ... 34% finished ... 35% finished ... 36% finished ... 37% finished ... 38% finished ... 39% finished ... 40% finished ... 41% finished ... 42% finished ... 43% finished ... 44% finished ... 45% finished ... 46% finished ... 47% finished ... 48% finished ... 49% finished ... 50% finished ... 51% finished ... 52% finished ... 53% finished ... 54% finished ... 55% finished ... 56% finished ... 57% finished ... 58% finished ... 59% finished ... 60% finished ... 61% finished ... 62% finished ... 63% finished ... 64% finished ... 65% finished ... 66% finished ... 67% finished ... 68% finished ... 69% finished ... 70% finished ... 71% finished ... 72% finished ... 73% finished ... 74% finished ... 75% finished ... 76% finished ... 77% finished ... 78% finished ... 79% finished ... 80% finished ... 81% finished ... 82% finished ... 83% finished ... 84% finished ... 85% finished ... 86% finished ... 87% finished ... 88% finished ... 89% finished ... 90% finished ... 91% finished ... 92% finished ... 93% finished ... 94% finished ... 95% finished ... 96% finished ... 97% finished ... 98% finished ... 99% finished ... 100% finished ...
[Thu Mar 18 10:42:25 2021] [8/20] Saving potential joining paths ...
[Thu Mar 18 10:42:25 2021] [9/20] Saving contig DNA break information ...
[Thu Mar 18 10:42:25 2021] [10/20] Getting retrieved contigs ...
[Thu Mar 18 10:42:25 2021] [11/20] Analyzing for valid joining paths ...
[Thu Mar 18 10:42:25 2021] [12/20] Saving valid joining paths ...
[Thu Mar 18 10:42:26 2021] [13/20] Getting initial joining status of each query contig ...
[Thu Mar 18 10:42:26 2021] [14/20] Getting final joining status of each query contig ...
[Thu Mar 18 10:43:02 2021] [15/20] Detecting invalid joins due to similar direct terminal repeats of different viruses ...
[Thu Mar 18 10:43:15 2021] [16/20] Saving joining summary of retrieved contigs ...
[Thu Mar 18 10:43:16 2021] [17/20] Saving unique sequences of "Extended_circular" and "Extended_10k" ...
[Thu Mar 18 10:44:02 2021] [18/20] Saving joining details of unique "Extended_circular" and "Extended_10k" queries ...
[Thu Mar 18 10:44:03 2021] [19/20] Saving self_circular contigs ...
[Thu Mar 18 10:44:03 2021] [20/20] Saving less dominant genomes of "joined_circular" ...
```

The ```log``` file also gives a summary of the joining status of all queries, example shown below:

```
Final summary
Total queries: 616
Self circular sequences: 27
Extended_circular: 70 (Unique: 48)
Extended >= 10k: 151 (Unique: 108)
Extended < 10k: 155
Failed due to COBRA rules: 93
Failed due to DNA break or low abundance: 120
Failed due to Unexpected_assembly_break_or_short_piece_missing: 0
```

## Citation
The manuscript is in preparation (Chen et al., in prep).
