# COBRA
COBRA (Contig Overlap Based Re-Assembly) is a bioinformatic tool to get higher quality virus genomes assembled from short-reads metagenomes.

## Introduction
Virus genomes assembled from short-reads sequenced metagenomes are usually fragmented due to intra-genome repeats and within-population variations (or subpopulation diversity, or local diversity), as the widely used assemblers based on de Bruijn graphs, e.g., metaSPAdes, IDBA_UD and MEGAHIT, tend to have a breaking point when multiple paths are available instead of making risky extension. According to the principles of the abovementioned assemblers, the broken contigs have an end overlap with determined length, that is the *max-kmer* used in de nono assembly for metaSPAdes and MEGAHIT, and the *max-kmer* - 1 for IDBA_UD, which we termed as ```expected overlap length (EOL)```. COBRA determines the EOL (both the forward direction and reverse complement direction) for all the contigs from an assembly, then looks for the valid joining path for each query that users provide (should be a fraction of the whole assembly) based on a list of features including contig coverage, contig overlap relationships, and contig continuity (based on paired end reads mapping).

*Note: as COBRA will use information provided by paired-end reads, thus only those samples sequenced by paired-end technology should work.*

<img src="https://user-images.githubusercontent.com/46725273/111421478-49768b00-86aa-11eb-8bea-9d4aa060a5e0.png" width="519" height="186">

Figure 1. The EOL have been documented in manual genome curation, see [Chen et al. 2020. Genome Research](https://genome.cshlp.org/content/30/3/315.short) for details.

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

```

For those queries resulting in the "extended_circular" category, if "equal_path" exists (Figure 2), COBRA will output both the ```dominant genome``` with the dominant path and the ```rare genome``` with the rare path (Figure 3).

![image](https://user-images.githubusercontent.com/46725273/111669292-e41eb900-87d3-11eb-8f23-b23b0b5cdb3c.png)

Figure 2. The equal pathes originating from within-population diversity (or local diversity).


![image](https://user-images.githubusercontent.com/46725273/111668676-390dff80-87d3-11eb-87e5-b16251f06b73.png)

Figure 3. Information of COBRA "extended_circular" genomes (dominant and rare).

### COBRA_joining_summary.txt file
All the queries in the categories of 

### log file
The ```log``` file gives a summary of the joining status of all queries, example shown below:

```
======================================================================================================================================================
Final summary
Total queries: 616
Self circular sequences: 27
Extended_circular: 70 (Unique: 48)
Extended >= 10k: 151 (Unique: 108)
Extended < 10k: 155
Failed due to COBRA rules: 93
Failed due to DNA break or low abundance: 120
Failed due to Unexpected_assembly_break_or_short_piece_missing: 0
======================================================================================================================================================
```

## Citation
The manuscript is in preparation (Chen et al., in prep).
