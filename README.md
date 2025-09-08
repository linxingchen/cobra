![COBRA_logo](https://github.com/linxingchen/cobra/assets/46725273/f21d1198-991c-42b2-9ec0-0ae8a662ba35)

COBRA (Contig Overlap Based Re-Assembly) is a bioinformatics tool to get higher quality viral genomes assembled from metagenomes of short paired-end reads. COBRA was written in Python. COBRA has so far only been tested on assembled contigs and scaffolds from metaSPAdes, IDBA_UD, and MEGAHIT.

```
# Developed by Dr. LinXing Chen
# University of California, Berkeley
# The Banfield Lab
# Email: linkingchan@gmail.com
```

## Versions
1. v1.2.2 (released on 2023-09-03) - initial release

2. v1.2.3 (released on 2024-02-26)
    - The GC function issue due to the update of Biopython.
    - The abnormal exit in the middle of processing some samples.
    - If none of the queries was extended, the process will break. If your runs do not have the expected output files, see the log file.

3. v1.3.0 (released on 2025-06-26)
    - refactor+feature: format code and add trim_readno
    - fix the handle of "6" shape path


## Citation
The paper is out at Nature Microbiology (https://www.nature.com/articles/s41564-023-01598-2). Please cite as follows if you find COBRA is helpful for your analyses.
```Chen, L., Banfield, J.F. COBRA improves the completeness and contiguity of viral genomes assembled from metagenomes. Nat Microbiol (2024). https://doi.org/10.1038/s41564-023-01598-2```


## Introduction

**1. Why are metagenomic contigs fragmented?**

The genomes assembled from short paired-end reads based metagenomes are usually fragmented due to (1) intra-genome repeats, (2) inter-genome shared region, and (3) within-population variations, as the widely utilized assemblers based on de Bruijn graphs, e.g., metaSPAdes, IDBA_UD and MEGAHIT, tend to have a breaking point when multiple paths are available instead of making risky extension (see example in **Figure 1**).

![image](https://user-images.githubusercontent.com/46725273/111676563-8a21f180-87db-11eb-9b8c-4c63fb993936.png)

Figure 1. An example of how assemblers break in assembly when within-population occurs.

**2. Contigs may be joined with expected end overlap.**

According to the principles of the abovementioned assemblers, the broken contigs have an end overlap with a determined length, that is, the max-kmer (maxK hereafter) used in de nono assembly for metaSPAdes and MEGAHIT, and the maxK-1 for IDBA_UD, which we termed as "expected overlap length" (Figures 1 and 2).

* Note: as COBRA will use the information provided by paired-end reads, only those samples sequenced by paired-end technology should work.

![image](https://user-images.githubusercontent.com/46725273/111677281-4c719880-87dc-11eb-85a9-a62906f4e10b.png)

Figure 2. The "expected overlap length" has been documented in manual genome curation, see [Chen et al. 2020](https://genome.cshlp.org/content/30/3/315.short) for details.


## How COBRA works

COBRA determines the "expected overlap length" (both the forward direction and reverse complement direction) for all the contigs from an assembly, then looks for the valid joining path for each query that users provide (should be a fraction of the whole assembly) based on a list of features including contig coverage, contig overlap relationships, and contig continuity (based on paired-end reads mapping) (Figure 3).

Note that scaffolds (for example, metaSPAdes assembly) could be used as input for COBRA extension; however, we suggest not using scaffolds from IDBA_UD as the potential errors in the scaffolding step (see [Chen et al. 2020](https://genome.cshlp.org/content/30/3/315.short) for details). Thus, for IDBA_UD and MEGAHIT assembly, the contigs should be used.

Given that COBRA has only tested for contigs/scaffolds from IDBA_UD, metaSPAdes, and MEGAHIT, it will be risky to use it on contigs/scaffolds from any other assemblers.

![Figure 1](https://github.com/linxingchen/cobra/assets/46725273/77a8285b-cce8-45f3-8e4e-4c46c3f45354)

Figure 3. The workflow of COBRA.


## Dependencies
* COBRA is a Python script (tested for version 3.7 or higher) that uses a list of frequently used Python packages, including:
```
Bio
Bio.Seq
collections
argparse
math
pysam
time
```

* The only third-party software that COBRA will use is [BLASTn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

## Installation
COBRA could now be installed in different ways.

* (1) git

```git clone https://github.com/linxingchen/cobra.git```

```cd cobra```

```python cobra.py -h```


* (2) pip

```pip install cobra-meta```

To confirm the installment,

```cobra-meta -h```

Which shows something like this

```
usage: cobra-meta [-h] -q QUERY -f FASTA -a {idba,megahit,metaspades} -mink MINK -maxk MAXK -m MAPPING -c COVERAGE [-lm LINKAGE_MISMATCH] [-o OUTPUT] [-t THREADS] [-v]

...
```

* (3) conda

```conda create -n cobra python=3.8```

```conda activate cobra```

```conda install bioconda::cobra-meta``` or ```conda install linxingchen1987::cobra-meta```

To confirm the installment,

```cobra-meta -h```

Which shows something like this

```
usage: cobra.py [-h] -q QUERY [-i IGNORE] -f FASTA -a {idba,megahit,metaspades} -mink MINK -maxk MAXK -m MAPPING
                [--mapping-link-cache [MAPPING_LINK_CACHE ...]] -c COVERAGE [-lm LINKAGE_MISMATCH] [-tr {no,trim,auto}]
                [--skip_new_assembly] [-o OUTPUT] [-t THREADS] [-v]

...
```

## Update
* pip

```pip install --upgrade cobra-meta```

* conda

``` conda activate cobra``` (if cobra is the conda environment name)

``` conda update cobra-meta```


## Input files
(1) COBRA needs four files as inputs, i.e.,

* ```-f/--fasta```: A fasta format file containing all the contigs from a single assembly. Note that IDBA_UD and MEGAHIT usually save contigs with a minimum length of 200 bp.

* ```-c/--coverage```: a two-column (separated by tab) file of the sequencing coverage of all contigs in the ```-f/--fasta``` file, example below:

```
contig-140_0    25.552
contig-140_1    42.1388
contig-140_2    14.6023
contig-140_3    15.4817
contig-140_4    41.2746
...
```

* ```-q/--query```: the query contigs that the user wants COBRA to extend, could be provided in a fasta format file, or a one-column text file with the names of the query contigs. Please make sure the names are exactly the same format as in the ```-f/--fasta``` file; otherwise, COBRA may have problems extending them.

* ```-m/--mapping```: the paired-end reads mapping file of all contigs in the ```-f/--fasta``` file, could be sam or bam file.

(2) and three parameters
* ```-a/--assembler```: the name of the de novo assembler used, currently only 'idba' (for IDBA_UD), 'metaspades' (for metaSPAdes), and 'megahit' (for MEGAHIT).
* ```-maxk/--maxk```: the largest kmer used in de novo assembly.
* ```-mink/--mink```: the smallest kmer used in de novo assembly.


(3) Optional flags
* ```-lm/--linkage_mismatch```: the number of read mapping mismatches allowed when determining if paired reads spanned two contigs.
* ```-o/--output```: the name of the output folder, otherwise it will be "{```-q/--query```}.COBRA" if not provided.
* ```-t/--threads```: the number of threads used for BLASTn search.


## How to obtain the mapping file

The mapping file could be obtained with tools like Bowtie2 and BBMap. Please refer to the manual descriptions for details of the tools. Below is the general way to get the sorted sam/bam file; you thus need to be available to samtools (which could be downloaded here - https://github.com/samtools/samtools).

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


## How to run
(1) The users can only specify the required parameters:

```
cobra-meta -f input.fasta -q query.fasta -c coverage.txt -m mapping.sam -a idba -mink 20 -maxk 140
```

(2) The users could also include the optional parameters like output name (-o), mismatch of mapped reads for linkage identification (-lm)

```
cobra-meta -f all.contigs.fasta -q query.fasta -o query.fasta.COBRA.out -c coverage.txt -m mapping.sam -a idba -mink 20 -maxk 140 -lm 2
```

```
cobra-meta -f all.contigs.fasta -q query.fasta -o query.fasta.COBRA.out -c coverage.txt -m mapping.sam -a metaspades -mink 21 -maxk 127 -lm 2
```

```
cobra-meta -f all.contigs.fasta -q query.fasta -o query.fasta.COBRA.out -c coverage.txt -m mapping.sam -a megahit -mink 21 -maxk 141 -lm 2
```


## Output files
Below is a general list of output files in the output folder:

```
COBRA_end_joining_pairs.tsv
COBRA_potential_joining_paths.tsv
blast_pairs.tsv_temp/
blast_pairs.tsv
debug.txt
COBRA_grouping_gv.tsv
COBRA_check_assembly_reason.tsv
COBRA_joining_summary.tsv
COBRA_joining_failed_paths.tsv
COBRA_joining_status.tsv
COBRA_category_ii-a_extended_circular_unique.fa
COBRA_category_ii-a_extended_circular_unique.fa.summary.tsv
COBRA_category_ii-b_extended_partial_unique.fa
COBRA_category_ii-b_extended_partial_unique.fa.summary.tsv
COBRA_category_i_self_circular.fa
COBRA_category_i_self_circular.fa.summary.tsv
COBRA_category_ii-c_extended_failed.fa
COBRA_category_ii-c_extended_failed.fa.summary.tsv
COBRA_category_iii-b_complex_end_remain.fa
COBRA_category_iii-b_complex_end_remain.fa.summary.tsv
COBRA_category_iii-a_orphan_end.fa
COBRA_category_iii-a_orphan_end.fa.summary.tsv
AcMG_contigs.fasta.new.fa
log
```

For all the queries, COBRA assigns them to different categories based on their joining status (detailed in the ```COBRA_joining_status.txt``` file), i.e.,

* "self_circular" - the query contig itself is a circular genome.
* "extended_circular" - the query contig was joined and extended into a circular genome.
* "extended_partial" - the query contig was joined and extended, but not circular.
* "extended_failed" - the query contig was not able to be extended due to COBRA rules.
* "orphan_end" - neither end of a given contig shares "expected overlap length" with others.
* "complex_end" - the query contig has multiple possible joining paths, so COBRA could not resolve them.

For the joined and extended queries in each category, only the unique ones (```*.fa```) will be saved for users' following analyses, and the sequence information (e.g., length, coverage, GC, num of Ns) is summarized in the ```*fa.summary.tsv``` files. For categories of "extended_circular", and "extended_partial", the joining details of each query are included in the corresponding folder and summarized in the ```COBRA_joining_summary.tsv``` file, an example shown below:

```
StatusReason	RepQuery	Status	Item	Contig	FinalSeqID	GroupID	JoinLen	StartOnJoin	Direction	JoinedReason	ContigLen	Cov	GC	IsQuery
circular	AcMG_13342	extended_circular	AcMG_4689667_Lrc AcMG_619867_Lrc AcMG_24500_R AcMG_13342	AcMG_4689667	AcMG_13342_extended_circular	165	6482	0	reverse	other_end_is_extendable	203	67.34694	0.6699507389162561	False
circular	AcMG_13342	extended_circular	AcMG_4689667_Lrc AcMG_619867_Lrc AcMG_24500_R AcMG_13342	AcMG_619867	AcMG_13342_extended_circular	165	6482	148	reverse	the_better_one	456	95.69927	0.6600877192982456	False
circular	AcMG_13342	extended_circular	AcMG_4689667_Lrc AcMG_619867_Lrc AcMG_24500_R AcMG_13342	AcMG_24500	AcMG_13342_extended_circular	165	6482	549	forward	the_better_one	2517	71.29798	0.5387365911799762	False
circular	AcMG_13342	extended_circular	AcMG_4689667_Lrc AcMG_619867_Lrc AcMG_24500_R AcMG_13342	AcMG_13342	AcMG_13342_extended_circular	165	6482	3011	forward	query	3526	116.24112	0.6154282473057289	True
circular	AcMG_14922	extended_circular	AcMG_18416_R AcMG_14922	AcMG_18416	AcMG_14922_extended_circular	238	6218	0	forward	the_better_one	2955	17.084356	0.3549915397631134	True
circular	AcMG_14922	extended_circular	AcMG_18416_R AcMG_14922	AcMG_14922	AcMG_14922_extended_circular	238	6218	2900	forward	query	3318	15.571528	0.3794454490657022	True
...
```


* **log file:** The ```log``` file includes the content of each processing step, an example shown below:

```
1. INPUT INFORMATION
# Max-kmer: 55
# Min-kmer: 21
# Overlap length: 55 bp
# Read mapping max mismatches for contig linkage: 2
# Query contigs: AcMG/AcMG_virus_gt2500_clean.fasta
# Whole contig set: AcMG/AcMG_contigs.fasta
# Mapping file: AcMG/AcMG.sorted.bam
# Coverage file: AcMG/AcMG
# Output folder: AcMG_virus_gt2500_clean.fasta_COBRA


2. PROCESSING STEPS
2.1. Loading assembly and mapping data
[0/5] [2024/11/03 11:51:54] Reading contigs and getting the contig end pairs. A total of 4775825 contigs were imported.
[1/5] [2024/11/03 11:53:14] Joining contigs by contig end (maxK). Among 570278 link pairs, found 407474 one path end, 158542 two paths end.
[2/5] [2024/11/03 11:53:18] Getting linkage based on sam/bam. Be patient, this may take long.
[3/5] [2024/11/03 11:55:27] Reading contig coverage information.
[4/5] [2024/11/03 11:55:30] Getting query contig list. A total of 3315 query contigs were imported, with 1495 query with unique end (orphan).
2.2. Analyzing assembled paths and solving conflicts
[0/8] [2024/11/03 11:55:30] Detecting self_circular contigs, independent of mapping linkage.
[1/8] [2024/11/03 11:55:31] Detecting joins of contigs.
[2/8] [2024/11/03 11:55:32] Saving potential joining paths.
[3/8] [2024/11/03 11:55:32] Getting the joining paths of contigs.
[4/8] [2024/11/03 11:55:32] Getting joined seqeuences.
[5/8] [2024/11/03 11:55:32] Checking for invalid joining using BLASTn: close strains.
[6/8] [2024/11/03 11:55:33] Grouping paths by sharing queries to check for invalid queries.
[7/8] [2024/11/03 11:55:33] Filtering paths according to COBRA rules.
3315																query_set
>	1495															orphan_end_query
>	0	486														complex_end_query
>	0	0	13													self_circular [overlap_maxk]
0	0	0	0	0												self_circular [overlap_flex]
>	0	0	0	0	1005											assembly_rep
>	0	73	0	0	132	305										failed_join_queries
>	0	68	0	0	70	>	149									failed_join [complex]
>	0	6	0	0	62	>	5	155								failed_join [conflict]
>	0	1	0	0	2	>	0	0	6							failed_join [circular_6]
>	0	0	1	0	1	1	0	0	1	2						failed_blast_half
0	0	0	0	0	0	0	0	0	0	0	0					redundant_circular_8
>	0	0	0	0	>	0	0	0	0	0	0	873				checked_rep
>	0	0	0	0	866	0	0	0	0	0	0	866	1078			checked_partial_queries
>	0	0	0	0	7	0	0	0	0	0	0	7	0	11		checked_circular_queries
>	0	0	0	0	>	0	0	0	0	0	0	>	0	>	7	path_circular_rep
2.3. Output extended paths and circulated paths
[0/5] [2024/11/03 11:55:33] Getting the joining details of extended query contigs.
[1/5] [2024/11/03 11:55:36] Getting the joining details of failed query contigs.
[2/5] [2024/11/03 11:55:36] Saving joining status of all query contigs.
[3/5] [2024/11/03 11:55:37] Saving identified and modified contigs.
[4/5] [2024/11/03 11:55:37] Summarise all query contigs and save the new assembled contigs.


3. RESULTS SUMMARY
# Total queries: 3315
# Category i   - self_circular: 13
# Category ii  - extended_circular: 11 (Unique: 7)
# Category ii  - extended_partial: 1078 (Unique: 866)
# Category ii  - extended_failed (due to COBRA rules): 305
# Category iii - orphan end: 1495
# Category iii - too complex to resolve: 413
# Check "COBRA_query_status.tsv" for joining status of each query.
# Check "COBRA_joining_summary.tsv" for joining details of "extended_circular" and "extended_partial" queries.
```
