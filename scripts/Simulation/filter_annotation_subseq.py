import sys
import time
from Bio.Seq import Seq
from utils import *

def acquire_reads(loci_pos, hg19, chr_list, pre_out):
	# fake_genome = list()
	file = open(pre_out, 'w')
	for key in chr_list:
		if key not in loci_pos:
			continue
		for ele in loci_pos[key]:
			# fake_genome.append(hg19[key].seq[ele[0]:ele[1], ele[0]])
			file.write(">%s*%d*%d*%s*%d\n"%(key, ele[0], ele[1], ele[2], ele[3]))
			file.write("%s\n"%(hg19[key].seq[ele[0]:ele[1]]))
	file.close()


def run(ref_g, Alu_l, pre_out):
	print("[INFO]: Loading reference genome ...")
	hg19 = load_ref(ref_g)
	print("[INFO]: Loading loci file ...")
	chr_list = acquire_chr(ref_g)
	loci_pos = acquire_loci(Alu_l, chr_list)
	acquire_reads(loci_pos, hg19, chr_list, pre_out)

def main():
	starttime = time.time()
	ref_g = sys.argv[1]
	Alu_l = sys.argv[2]
	pre_out = sys.argv[3]
	print("[INFO]: 1. Reference Genome: %s"%(ref_g))
	print("[INFO]: 2. Alu locus annotation: %s"%(Alu_l))
	print("[INFO]: 3. The prefix path of output files: %s"%(pre_out))
	print("[INFO]: This script will waste your a lot of time. sorry!")
	run(ref_g, Alu_l, pre_out)
	print("[INFO]: Finished in %0.2f seconds."%(time.time() - starttime))
	print("[INFO]: Have a nice day! ^.^ ")

if __name__ == '__main__':
	main()