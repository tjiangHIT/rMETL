import sys
import pysam
import tools
from process import *

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
	tools.load_sam(args)
	print("[INFO]: Finished in %0.2f seconds."%(time.time() - starttime))
	print("[INFO]: Have a nice day! ^.^ ")

if __name__ == '__main__':
	main()