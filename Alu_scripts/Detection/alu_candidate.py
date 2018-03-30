import sys
import pysam
import tools

def main():
	import time
	starttime = time.time()
	sam_path = sys.argv[1]
	signal_path = sys.argv[2]
	reference_path = sys.argv[3]
	print("[INFO]: The path of the Alignment File: %s"%(sam_path))
	print("[INFO]: The path of the Signal File: %s"%(signal_path))
	print("[INFO]: The path of the Reference File: %s"%(reference_path))
	tools.load_sam(sam_path, signal_path, reference_path)
	print("[INFO]: Finished in %0.2f seconds."%(time.time() - starttime))
	print("[INFO]: Have a nice day! ^.^ ")

if __name__ == '__main__':
	main()