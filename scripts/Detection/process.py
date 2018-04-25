#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  process.py
 * @Package: argparse, pysam, sys
 * @Description: Control the nTED pipeline
 * @author: tjiang
 * @date: Apr 24 2018
 * @version V1.0     
'''

import argparse
import sys
import pysam
import extract, Map, call_TE
# from process import *

STAGES = {"extract": extract.run, \
          "map": Map.run, \
          "call": call_TE.run}

USAGE = """\
   nTED - Non-reference transposable elements detection with long sequencing reads

   STAGE is one of
   	 extract    detect non-ref TE signals
   	 map        aligne signal seqs to TE librarys
   	 call       cluster TE seqs and genotyping
    
   See Readme.txt for documentation or --help for details\
"""

def parseArgs():

	parser = argparse.ArgumentParser(prog="process.py", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)

	# parser.add_argument("-h", "--help", action="store_true")
	parser.add_argument("stage", metavar="STAGE", choices=STAGES.keys(), type=str, help="Stage to execute")
	parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER, help="Options to pass to the stage")

	args = parser.parse_args()

	STAGES[args.stage](args.options)

if __name__ == '__main__':
	parseArgs()