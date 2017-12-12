import sys

def acquire_clip(cigar, flag):
	contig_num = len(cigar.split(flag))
	if contig_num == 1:
		return 0
	if contig_num == 2:
		if cigar[-1] == flag:
			return int(cigar.split(flag)[0].split('M')[-1].split('D')[-1].split('I')[-1])
		else:
			return int(cigar.split(flag)[0])
	if contig_num == 3:
		first = int(cigar.split(flag)[0])
		last = int(cigar.split(flag)[1].split('M')[-1].split('D')[-1].split('I')[-1])
		return first + last

def analysis_read(path):
	unaligned_dic = {}
	aligned_dic = {}
	# dic = {}
	file = open(path)
	for line in file:
		seq = line.strip()
		if seq[0] == '@':
			continue
		else:
			seq = seq.split('\t')
			cigar = seq[5]
			read_id = seq[0]
			read_len = int(len(seq[9]))

			if cigar[0] == '*':
				# unaligned_num += 1
				if read_id not in unaligned_dic:
					unaligned_dic[read_id] = 1
			else:
				soft_clip = acquire_clip(cigar, 'S')
				hard_clip = acquire_clip(cigar, 'H')
				clip = soft_clip + hard_clip
				if clip >= read_len * 0.1:
					# print read_id
					# unaligned_num += 1
					if read_id not in unaligned_dic:
						unaligned_dic[read_id] = 1
				else:
					if read_id not in aligned_dic:
						aligned_dic[read_id] = 1
					# aligned_num += 1
				# aligned_num += 1
	file.close()
	# print unaligned_num, aligned_num

	del_read_id = []
	for key in unaligned_dic:
		if key in aligned_dic:
			del_read_id.append(key)
			# unaligned_dic.pop(key)
	for i in xrange(len(del_read_id)):
		unaligned_dic.pop(del_read_id[i])

	# print len(unaligned_dic), len(aligned_dic), len(unaligned_dic)+len(aligned_dic)
	return unaligned_dic

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
