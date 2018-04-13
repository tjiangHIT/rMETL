import argparse
import sys
import pysam
import tools, Map, call_TE
# from process import *

def main():
	import time
	starttime = time.time()
	args = opt()
	# sam_path = sys.argv[1]
	# signal_path = sys.argv[2]
	# reference_path = sys.argv[3]
	alignment = args.AlignmentFile
	output = args.Output_prefix
	ref = args.Reference

	print("[INFO]: The path of the Alignment File: %s"%(alignment))
	print("[INFO]: The path of the Signal File: %s"%(output))
	print("[INFO]: The path of the Reference File: %s"%(ref))
	# tools.load_sam(args)
	tools.load_sam_multi_processes(args)
	print("[INFO]: Finished in %0.2f seconds."%(time.time() - starttime))
	print("[INFO]: Have a nice day! ^.^ ")

STAGES = {"extract": tools.run, \
          "map": Map.run, \
          "call": call_TE.run}

USAGE = """\
   nTED - Non-reference transposable element detecting with long sequencing reads

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