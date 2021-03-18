# COBRA
COBRA (Contig Overlap Based Re-Assembly) is a bioinformatic tool to get higher quality virus genomes assembled from short-reads metagenomes.

## Introduction
Virus genomes assembled from short-reads sequenced metagenomes are usually fragmented due to intra-genome repeats and within-population variations (or subpopulation diversity, local diversity), as the widely used assemblers based on de Bruijn graphs, e.g., metaSPAdes, IDBA_UD and MEGAHIT, tend to have a breaking point when multiple paths are available insteal of making risky extension. According to the principles of the abovementioned assemblers, the broken contigs have an end overlap with determined length, that is the ```max-kmer``` used in de nono assembly for metaSPAdes and MEGAHIT, and the ```max-kmer - 1``` for IDBA_UD, which we termed as ```expected overlap length (EOL)```. COBRA determines the EOL (both the forward direction and reverse complement direction) for all the contigs from an assembly, then looks for the valid joining path for each query that users provide (should be a fraction of the whole assembly) based on a list of features including contig coverage, contig overlap relationships, and contig continuity (based on paired end reads mapping).

![image](https://user-images.githubusercontent.com/46725273/111421478-49768b00-86aa-11eb-8bea-9d4aa060a5e0.png)

Figure 1. The EOL have been documented in manual genome curation, see [Chen et al. 2020. Genome Research](https://genome.cshlp.org/content/30/3/315.short) for details.

## Installation options
COBRA is a python script (version 3.7 or higher) that uses a list of frequently used python packages including:
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

input 1. ```all.contigs.fasta``` - the whole contig set from the assembly.

input 2. ```coverage.txt``` - a two columns (tab) file of the coverage of all contigs, example below:

```contig-140_0    25.552
contig-140_1    42.1388
contig-140_2    14.6023
contig-140_3    15.4817
contig-140_4    41.2746
...
```
* Note that metaSPAdes contigs have coverage information (not absolute coverage though) in their headers, which will be used by COBRA if no ```coverage.txt``` file is provided.

input 3. ```queries.fasta``` - the fasta file containing all the query contigs for joining path search.

input 4. ```mapping.sam (or mapping.bam)``` - the paired-end reads mapping file of all contigs.

Two additional flags should be provided for COBRA to determine the EOL, i.e.,

```assembler``` - currently only 'idba' (for IDBA_UD), 'metaspades' (for metaSPAdes), and 'megahit' (for MEGAHIT).

```kmer``` - the largest kmer used in de novo assembly.


## How to run

The output file will be ```queries.fasta.COBRA``` if not specified via the ```-o``` flag.

```COBRA.py -f all.contigs.fasta -c coverage.txt -q queries.fasta -m mapping.sam -a idba -k 100```

If assembled with metaSPAdes, the coverage information in the contig headers could be used when no additional coverage file is provided.

```COBRA.py -f all.contigs.fasta -q queries.fasta -m mapping.sam -a metaspades -k 127```


