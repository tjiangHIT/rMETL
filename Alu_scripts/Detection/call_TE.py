
import pysam

flag_dic = {0:1, 16:1, 2048:0, 2064:0}

cluster_dic = {}

def parse_name(seq):
	chr = seq.split('_')[0]
	breakpoint = seq.split('_')[1]
	insert_size = seq.split('_')[2]
	return chr, breakpoint, insert_size

def load_sam(p1):
	# samfile = pysam.AlignmentFile(p1)
	AlignmentFile = open(p1, 'r')
	for line in AlignmentFile:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '@':
			continue
		chr, breakpoint, insert_size = seq[0]
		Flag = int(seq[1])
		sub_type = seq[2]
		if flag_dic[Flag] == 1:
			# to do something

	AlignmentFile.close()