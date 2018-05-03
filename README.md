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



---
### Simulated datasets

The simulated datasets use for benchmarking are available at: https://drive.google.com/open?id=1ujV2C8e1PNAVhSkh9vKtjWLdG_OHcH-k

---
### Memory usage

The memory usage of rMETL can fit the configurations of most modern servers and workstations.
Its peak memory footprint is about 16.7G Gigabytes (default setting), on a server with Intel Xeon CPU at 2.00 GHz, 1 Terabytes RAM running Linux Ubuntu 14.04. These reads were aligned to human reference genome hs37d5.

---
### Dependences
	
	1. pysam
	2. Biopython
	3. ngmlr
	4. samtools

	Python version 2.7

---
### Installation

Current version of rMETL needs to be run on Linux operating system.
The source code is written in python, and can be directly download from: https://github.com/hitbc/rMETL 
A mirror is also in: https://github.com/tjiangHIT/rMETL
The INSTALL.sh is attached. Use the make command for generating the executable file.

---
### Synopsis
Detect Mobile Element signals.

	rMETL.py extract <alignments> <reference> <temp_dir> <output>

Realigne signal loci to transposable element concensus librarys.

	rMETL.py map <FASTA> <TEREF> <output>

Mobile Element calling and genotyping.

	rMETL.py call <SAM> <reference> <out_type> <output>

---
### Optional Parameters

#### extract

| Parameters | Descriptions | Defaults |
| :------------ |:---------------|:---------------|
| MIN_SUPPORT   |Mininum number of reads that support a ME.| 5 |
| MIN_LENGTH    | Mininum length of ME to be reported.        |50|
| MIN_DISTANCE  | Mininum distance of two ME clusters. |20|
| THREADS       |Number of threads to use.|1|
| PRESETS       |The sequencing type <pacbio,ont> of the reads.|pacbio|

#### map

| Parameters | Descriptions | Defaults |
| :------------ |:---------------|:---------------|
| THREADS       |Number of threads to use.|1|
| PRESETS       |The sequencing type <pacbio,ont> of the reads.|pacbio|
| SUBREAD_LENGTH       |Length of fragments reads are split into.|128|
| SUBREAD_CORRIDOR       |Length of corridor sub-reads are aligned with.|20|

#### call

| Parameters | Descriptions | Defaults |
| :------------ |:---------------|:---------------|
| HOMOZYGOUS       |The mininum score of a genotyping reported as a homozygous.|0.8|
| HETEROZYGOUS       |The mininum score of a genotyping reported as a heterozygous.|0.3|
| MIN_MAPQ       |Mininum mapping quality.|20|
| SAMPLE       |The name of the sample which be noted.|None|
| MEI       |Enables rMETL to display MEI/MED only.|False|

---
### Reference


---
### Contact
For advising, bug reporting and requiring help, please contact ydwang@hit.edu.cn or tjiang@hit.edu.cn
