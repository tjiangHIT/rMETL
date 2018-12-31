#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  process.py
 * @Package: argparse
 * @Description: Control the rMETL pipeline
 * @author: Jiang Tao (tjiang@hit.edu.cn)
 * @date: Apr 24 2018
 * @version V1.0.2   
'''

import argparse
import rMETL_extraction
import rMETL_realign
import rMETL_MEIcalling
from rMETL_version import __version__, __author__, __contact__

STAGES = {'detection': rMETL_extraction.run, \
          'realignment': rMETL_realign.run, \
          'calling': rMETL_MEIcalling.run}

USAGE = '''\
           _  ___  _   _____   _______   _
     _ _  | ^_   _^ | |  ___| |__   __| | |
    | ^_| | | | | | | | |__      | |    | |
    | |   | | | | | | |  __|     | |    | |
    | |   | | | | | | | |___     | |    | |___
    |_|   |_| |_| |_| |_____|    |_|    |_____|

    rMETL - realignment-based Mobile Element insertion detection Tool for Long read

  STAGE is one of
    detection    Inference of putative MEI loci.
    realignment  Realignment of chimeric read parts.
    calling      Mobile Element Insertion calling.
    
  See README.md for documentation or --help for details
  
  rMETL V%s
  Author: %s
  Contact: %s
'''%(__version__, __author__, __contact__)

def parseArgs():
	parser = argparse.ArgumentParser(prog='rMETL.py', description=USAGE, \
    formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('stage', metavar='STAGE', choices=STAGES.keys(), \
    type=str, help='Stage to execute')
	parser.add_argument('options', metavar='OPTIONS', nargs=argparse.REMAINDER, \
    help='Options to pass to the stage')
	args = parser.parse_args()
	STAGES[args.stage](args.options)

if __name__ == '__main__':
	parseArgs()