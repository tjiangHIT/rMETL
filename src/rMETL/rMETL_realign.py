#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  Map.py
 * @Package: argparse, sys, logging
 * @Description: Classify the ME types
  * @author: Jiang Tao (tjiang@hit.edu.cn)
 * @date: Apr 24 2018
 * @version V1.0.2
'''

import argparse
import logging
import sys
import time

from rMETL.rMETL_version import __version__, __author__, __contact__
from rMETL.rMETL_cmdRunner import setupLogging, exe

USAGE="""\
           _  ___  _   _____   _______   _
     _ _  | ^_   _^ | |  ___| |__   __| | |
    | ^_| | | | | | | | |__      | |    | |
    | |   | | | | | | |  __|     | |    | |
    | |   | | | | | | | |___     | |    | |___
    |_|   |_| |_| |_| |_____|    |_|    |_____|

    rMETL - realignment-based Mobile Element insertion detection Tool for Long read

	Realignment of chimeric read parts.

	Aligner: NGMLR version 0.2.6
	TE refs: Alu concensus
		 L1 concensus
		 SVA concensus
	The output is a sam format file called 'cluster.sam'.

	rMETL V%s
	Author: %s
	Contact: %s
"""%(__version__, __author__, __contact__)

# **************************Call-NGMLR********************************
def call_ngmlr(inFile, ref, presets, nproc, outFile, SUBREAD_LENGTH, SUBREAD_CORRIDOR):
	"""
	fq = input file
	automatically search for .sa
	"""
	outFile = outFile + "cluster.sam"
	logging.info("Running NGMLR...")
	cmd = ("ngmlr -r %s -q %s -o %s -t %d -x %s --subread-length %d --subread-corridor %d" \
		% (ref, inFile, outFile, nproc, presets, SUBREAD_LENGTH, SUBREAD_CORRIDOR))
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
# **************************Call-NGMLR********************************
# 
# 
# 
# ************************MAIN_FUNCTION*******************************
def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="rMETL.py realignment", description=USAGE, \
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("input", metavar="FASTA", type=str, help="Input potential_ME.fa on STAGE detection.")
	parser.add_argument("ME_Ref", type=str, help="The transposable element concensus in fasta format.")
	parser.add_argument('output', type=str, help = "Directory to output realignments.")
	parser.add_argument('-t', '--threads', help = "Number of threads to use.[%(default)s]", \
		default = 8, type = int)
	parser.add_argument('-x', '--presets', \
		help = "The sequencing platform <pacbio,ont> of the reads.[%(default)s]", \
		default = "pacbio", type = str)
	parser.add_argument('--subread_length', \
		help = "Length of fragments reads are split into [%(default)s]", \
		default = 128, type = int)
	parser.add_argument('--subread_corridor', \
		help = "Length of corridor sub-reads are aligned with [%(default)s]", \
		default = 20, type = int)
	args = parser.parse_args(argv)
	return args

def run(argv):
	args = parseArgs(argv)
	setupLogging(False)
	starttime = time.time()
	call_ngmlr(args.input, args.ME_Ref, args.presets, args.threads, args.output, \
		args.subread_length, args.subread_corridor)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

if __name__ == '__main__':
    run(sys.argv[:1])
