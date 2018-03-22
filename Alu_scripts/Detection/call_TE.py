import sys	
import pysam
from collections import Counter

def acquire_count_max(_list_):
	c = Counter(_list_)
	return c.most_common(1)[0][0]

flag_dic = {0:1, 16:1, 2048:0, 2064:0, 4:0}

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
		chr, breakpoint, insert_size = parse_name(seq[0])
		Flag = int(seq[1])
		sub_type = seq[2]
		if flag_dic[Flag] == 1:
			# to do something
			key = "%s_%s_%s"%(chr, breakpoint, insert_size)
			if key not in cluster_dic:
				cluster_dic[key] = list()
			cluster_dic[key].append(sub_type)
	AlignmentFile.close()

	sort_list = list()
	for i in cluster_dic:
		chr, breakpoint, insert_size = parse_name(i)
		final_type = acquire_count_max(cluster_dic[i])
		sort_list.append([chr, breakpoint, insert_size, final_type])
		# print("%s\t%s\t%s\t%s"%(chr, breakpoint, insert_size, final_type))
	sort_list = sorted(sort_list, key = lambda x:(x[0], int(x[1])))
	for i in sort_list:
		# print("%s\t%s\t%s\t%s"%(i[0], breakpoint, insert_size, final_type))
		print "\t".join(i)

if __name__ == '__main__':
	path = sys.argv[1]
	load_sam(path)