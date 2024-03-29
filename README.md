		        __        __   ______   _________   _
		 _ __  |  \      /  | |  ____| |___   ___| | |
		| ^__| |   \    /   | | |___       | |     | |
		| |    | |\ \  / /| | |  ___|      | |     | |
		| |    | | \ \/ / | | | |____      | |     | |____
		|_|    |_|  \__/  |_| |______|     |_|     |______|
     
	rMETL - realignment-based Mobile Element insertion detection Tool for Long read
 **NOTE: The community users give the newest installation approach after 2023, which is referred to** [here](https://github.com/tjiangHIT/rMETL/issues/8).

[![PyPI version](https://badge.fury.io/py/rMETL.svg)](https://badge.fury.io/py/rMETL)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/rmetl/badges/version.svg)](https://anaconda.org/bioconda/rmetl)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/rmetl/badges/license.svg)](https://anaconda.org/bioconda/rmetl)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/rmetl/badges/platforms.svg)](https://anaconda.org/bioconda/rmetl)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/rmetl/badges/latest_release_date.svg)](https://anaconda.org/bioconda/rmetl)

---	
### Introduction
Mobile element insertion (MEI) is a significant category of structure variations (SVs). The rapid development of long-read sequencing technologies provides the opportunity to detect MEIs sensitively. However, the signals of MEI implied by noisy long reads are highly complex due to the repetitiveness of mobile elements and the high sequencing error rates. Herein, we propose the Realignment-based Mobile Element insertion detection Tool for Long read (rMETL). Benchmarking results of simulated and real datasets demonstrate that rMETL has the ability to discover MEIs sensitively as well as prevent false positives. It is suited to produce high-quality MEI callsets in many genomics studies.

---
### Simulated datasets

The simulated datasets used for benchmarking are available at [Google Drive](https://drive.google.com/open?id=1ujV2C8e1PNAVhSkh9vKtjWLdG_OHcH-k)

---
### Memory usage

The memory usage of rMETL can fit the configurations of most modern servers and workstations.
Its peak memory footprint is about 7.05 Gigabytes (default setting), on a server with Intel Xeon CPU at 2.00 GHz, 1 Terabytes RAM running Linux Ubuntu 14.04. These reads were aligned to human reference genome hs37d5.

---
### Dependences
	
	1. pysam
	2. Biopython
	3. ngmlr
	4. samtools
	5. cigar

	Python version 2.7

---
### Installation

	#install via pip
	$ pip install rMETL
	
	#install via conda
	$ conda install -c bioconda rmetl

	#install from GitHub
	$ git clone https://github.com/tjiangHIT/rMETL.git (git clone https://github.com/hitbc/rMETL.git)
	$ cd rMETL/
	$ pip install .

The current version of rMETL has been tested on a 64-bit Linux operating system.

**NOTE: The community users give the newest installation approach after 2023, which is referred to** [here](https://github.com/tjiangHIT/rMETL/issues/8).

---
### Synopsis
Inference of putative MEI loci.

	rMETL.py detection <alignments> <reference> <temp_dir> <output>

Realignment of chimeric read parts.

	rMETL.py realignment <FASTA> <MEREF> <output>

Mobile Element Insertion calling.

	rMETL.py calling <SAM> <reference> <out_type> <output>
	
Strongly recommend making the output directory manually at first.:blush:

---
### Optional Parameters

#### Detection

| Parameters | Descriptions | Defaults |
| :------------ |:---------------|:---------------|
| MIN_SUPPORT   |Mininum number of reads that support a ME.| 5 |
| MIN_LENGTH    | Minimum length of ME to be reported.        |50|
| MIN_DISTANCE  | Minimum distance of two ME clusters. |20|
| THREADS       |Number of threads to use.|1|
| PRESETS       |The sequencing type <pacbio,ont> of the reads.|pacbio|

#### Realignment

| Parameters | Descriptions | Defaults |
| :------------ |:---------------|:---------------|
| THREADS       |Number of threads to use.|1|
| PRESETS       |The sequencing type <pacbio,ont> of the reads.|pacbio|
| SUBREAD_LENGTH       |Length of fragments reads are split into.|128|
| SUBREAD_CORRIDOR       |Length of corridor sub-reads are aligned with.|20|

#### Calling

| Parameters | Descriptions | Defaults |
| :------------ |:---------------|:---------------|
| HOMOZYGOUS       |The minimum score of a genotyping reported as homozygous.|0.8|
| HETEROZYGOUS       |The minimum score of a genotyping reported as a heterozygous.|0.3|
| MIN_MAPQ       |Mininum mapping quality.|20|
| CLIPPING_THRESHOLD  |Mininum threshold of realignment clipping.|0.5|
| SAMPLE       |The name of the sample which is noted.|None|
| MEI       |Enables rMETL to display MEI/MED only.|False|

---
### Citation
If you use rMETL, please cite:
> Tao Jiang *et al*; rMETL: sensitive mobile element insertion detection with long read realignment, *Bioinformatics*, Volume 35, Issue 18, 15 September 2019, Pages 3484–3486, https://doi.org/10.1093/bioinformatics/btz106

---
### Contact
For advising, bug reporting, and requiring help, please post on [Github Issue](https://github.com/tjiangHIT/rMETL/issues) or contact tjiang@hit.edu.cn.
