# COBRA
COBRA (Contig Overlap Based Re-Assembly) is a bioinformatic tool to get higher quality virus genomes assembled from short-reads metagenomes.

## Introduction
Virus genomes assembled from short-reads sequenced metagenomes are usually fragmented due to intra-genome repeats and within-population variations (or subpopulation diversity, local diversity), as the widely used assemblers based on de Bruijn graphs, e.g., metaSPAdes, IDBA_UD and MEGAHIT, tend to have a breaking point when multiple paths are available insteal of making risky extension. According to the principles of the abovementioned assemblers, the broken contigs have an end overlap with determined length, that is the ```max-kmer``` used in de nono assembly for metaSPAdes and MEGAHIT, and the ```max-kmer - 1``` for IDBA_UD, which we termed as ```expected overlap length (EOL)```. COBRA determines the EOL (both the forward direction and reverse complement direction) for all the contigs from an assembly, then looks for the valid joining path for each query that users provide (should be a fraction of the whole assembly) based on a list of features including contig coverage, contig overlap relationships, and contig continuity (based on paired end reads mapping).

Manual extension.png![image](https://user-images.githubusercontent.com/46725273/111420720-f6e89f00-86a8-11eb-98be-0838f20b06eb.png)
Figure 1. The EOL have been documented in manual genome curation, see Chen et al. 2020 (https://genome.cshlp.org/content/30/3/315.short) for details.

## Installation options
COBRA is a python script that uses a list of frequently used python packages including:
```
Bio
Bio.Seq
collections
argparse
math
pysam
time
```

The only third-party software that the user should have is ```BLASTn```, which is available at https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download.
