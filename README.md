# rMETL: 
rMETL - realignment-based Mobile Element insertion detection Tool for Long read

---
### Getting Start
		        __        __   ______   _________   _
		 _ __  |  \      /  | |  ____| |___   ___| | |
		| ^__| |   \    /   | | |___       | |     | |
		| |    | |\ \  / /| | |  ___|      | |     | |
		| |    | | \ \/ / | | | |____      | |     | |____
		|_|    |_|  \__/  |_| |______|     |_|     |______|
     
	
	$ git clone https://github.com/hitbc/rMETL.git (git clone https://github.com/tjiangHIT/rMETL.git)
	$ cd rMETL/
	$ bash INSTALL.sh
	$ ./rMETL.py

---	
### Introduction
Mobile element insertion (MEI) is a major category of structure variations (SVs). The rapid development of long read sequencing provides the opportunity to sensitively discover MEIs. However, the signals of MEIs implied by noisy long reads are highly complex, due to the repetitiveness of mobile elements as well as the serious sequencing errors. Herein, we propose Realignment-based Mobile Element insertion detection Tool for Long read (rMETL). rMETL takes advantage of its novel chimeric read re-alignment approach to well handle complex MEI signals. Benchmarking results on simulated and real datasets demonstrated that rMETL has the ability to more sensitivity discover MEIs as well as prevent false positives. It is suited to produce high quality MEI callsets in many genomics studies.


---
### Simulated datasets

The simulated datasets use for benchmarking are available at: https://drive.google.com/open?id=1ujV2C8e1PNAVhSkh9vKtjWLdG_OHcH-k

---
### Memory usage

The memory usage of rMETL can fit the configurations of most modern servers and workstations.
Its peak memory footprint is about 12.18 Gigabytes (default setting), on a server with Intel Xeon CPU at 2.00 GHz, 1 Terabytes RAM running Linux Ubuntu 14.04. These reads were aligned to human reference genome hs37d5.

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

Current version of rMETL needs to be run on Linux operating system.
The source code is written in python, and can be directly download from: https://github.com/hitbc/rMETL 
A mirror is also in: https://github.com/tjiangHIT/rMETL
The INSTALL.sh is attached. Use the bash command for generating the executable file.

---
### Synopsis
Inference of putative MEI loci.

	rMETL.py detection <alignments> <reference> <temp_dir> <output>

Realignment of chimeric read parts.

	rMETL.py realignment <FASTA> <MEREF> <output>

Mobile Element Insertion calling.

	rMETL.py calling <SAM> <reference> <out_type> <output>
	
Strongly recommend making output directory manually at first.:blush:

---
### Optional Parameters

#### Detection

| Parameters | Descriptions | Defaults |
| :------------ |:---------------|:---------------|
| MIN_SUPPORT   |Mininum number of reads that support a ME.| 5 |
| MIN_LENGTH    | Mininum length of ME to be reported.        |50|
| MIN_DISTANCE  | Mininum distance of two ME clusters. |20|
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
| HOMOZYGOUS       |The mininum score of a genotyping reported as a homozygous.|0.8|
| HETEROZYGOUS       |The mininum score of a genotyping reported as a heterozygous.|0.3|
| MIN_MAPQ       |Mininum mapping quality.|20|
| CLIPPING_THRESHOLD  |Mininum threshold of realignment clipping.|0.5|
| SAMPLE       |The name of the sample which be noted.|None|
| MEI       |Enables rMETL to display MEI/MED only.|False|

---
### Citation
Tao Jiang, Bo Liu, Junyi Li, Yadong Wang; rMETL: sensitive mobile element insertion detection with long read realignment, Bioinformatics, , btz106, https://doi.org/10.1093/bioinformatics/btz106

---
### Contact
For advising, bug reporting and requiring help, please post on Github Issue or contact tjiang@hit.edu.cn.
