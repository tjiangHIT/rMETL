#!/usr/bin/env python
# -*-coding:utf-8-*-

import logging
import os
from Bio import SeqIO
from rMETL_cmdRunner import exe

def load_ref(ref_g):
	logging.info("Loading reference genome...")
	return SeqIO.to_dict(SeqIO.parse(ref_g, "fasta"))

def check_bai(file, tempdir):
	'''
	check the index of a BAM file.
	'''
	if os.path.exists(file+".bai"):
		logging.info("The bam file is legal.")
		return ""
	else:
		logging.info("The bam.bai is missed.")
		logging.info("Running Samtools sort...")
		bam_path = file[:-3] + "sorted.bam"
		cmd = ("samtools sort -@ 4 -O bam -T %s -o %s %s" % (tempdir, \
			bam_path, file))
		r, o, e = exe(cmd)
		if r != 0:
			logging.error("Samtools sort failed!")
			logging.error("RETCODE %d" % (r))
			logging.error("STDOUT %s" % (str(o)))
			logging.error("STDERR %s" % (str(e)))
			logging.error("Exiting")
			exit(r)
		logging.info("Finished Samtools sort.")

		logging.info("Running Samtools index...")
		cmd = ("samtools index %s" % (bam_path))
		r, o, e = exe(cmd)
		if r != 0:
			logging.error("Samtools index failed!")
			logging.error("RETCODE %d" % (r))
			logging.error("STDOUT %s" % (str(o)))
			logging.error("STDERR %s" % (str(e)))
			logging.error("Exiting")
			exit(r)
		logging.info("Finished Samtools index.")
		return bam_path

def call_ngmlr(inFile, ref, presets, nproc, outFile):
	"""
	run ngmlr to generate alignments
	"""
	outFile = outFile + "map.sam"
	logging.info("Running NGMLR...")
	cmd = ("ngmlr -r %s -q %s -o %s -t %d -x %s" % (ref, inFile, outFile, \
		nproc, presets))
	r, o, e = exe(cmd)
	
	if r != 0:
		logging.error("NGMLR mapping failed!")
		logging.error("RETCODE %d" % (r))
		logging.error("STDOUT %s" % (str(o)))
		logging.error("STDERR %s" % (str(e)))
		logging.error("Exiting")
		exit(r)
	logging.info("Finished NGMLR mapping.")
	return outFile

def call_samtools(file, tempdir):
	'''
	run samtools to generate sorted BAM files.
	'''
	logging.info("Running Samtools sort...")
	bam_path = file[:-3] + "bam"
	cmd = ("samtools view -Sb %s | samtools sort -@ 4 -O bam -T %s - > %s" % \
		(file, tempdir, bam_path))
	r, o, e = exe(cmd)
	if r != 0:
		logging.error("Samtools sort failed!")
		logging.error("RETCODE %d" % (r))
		logging.error("STDOUT %s" % (str(o)))
		logging.error("STDERR %s" % (str(e)))
		logging.error("Exiting")
		exit(r)
	logging.info("Finished Samtools sort.")

	logging.info("Running Samtools index...")
	cmd = ("samtools index %s" % (bam_path))
	r, o, e = exe(cmd)
	if r != 0:
		logging.error("Samtools index failed!")
		logging.error("RETCODE %d" % (r))
		logging.error("STDOUT %s" % (str(o)))
		logging.error("STDERR %s" % (str(e)))
		logging.error("Exiting")
		exit(r)
	logging.info("Finished Samtools index.")
	return bam_path
