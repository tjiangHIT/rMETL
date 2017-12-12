import sys

def analysis_read(path):
	unaligned_num = 0
	aligned_num = 0
	dic = {}
	file = open(path)
	for line in file:
		seq = line.strip()
		if seq[0] == '@':
			continue
		else:
			seq = seq.split('\t')
			cigar = seq[5]
			read_id = seq[0]
			if cigar[0] == '*':
				unaligned_num += 1
				if read_id not in dic:
					dic[read_id] = 1
			else:
				aligned_num += 1
	file.close()
	print unaligned_num, aligned_num
	return dic

def main():
	unaligned_read_id_minimap2 = analysis_read(sys.argv[1])
	unaligned_read_id_nglmr = analysis_read(sys.argv[2])

	for key in unaligned_read_id_nglmr:
		if key not in unaligned_read_id_minimap2:
			print "belongs to nglmr: "+key

	for key in unaligned_read_id_minimap2:
		if key not in unaligned_read_id_nglmr:
			print "belongs to minimap2: "+key


if __name__ == '__main__':
	main()
