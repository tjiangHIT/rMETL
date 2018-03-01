
import pysam

# list_flag = {1:'I', 4:'S', 5:'H'}
list_flag = {1:'I'}
low_bandary = 20

def revcom_complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)[::-1]

def merge_pos(pos_list):
	start = list()
	end = list()
	for ele in pos_list:
		start.append(ele[0])
		end.append(ele[0] + ele[1])
	return [int(sum(start)/len(pos_list)), int(sum(end)/len(pos_list)) - int(sum(start)/len(pos_list)), len(pos_list)]


def cluster(pos_list):
	_cluster_ = list()
	temp = list()
	temp.append(pos_list[0])
	for pos in pos_list[1:]:
		if temp[-1][0] + temp[-1][1] < pos[0]:
			_cluster_.append(merge_pos(temp))
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	_cluster_.append(merge_pos(temp))
	return _cluster_

def detect_flag(Flag):
	# Signal
	Normal_foward = 1 >> 1
	Abnormal = 1 << 2
	Reverse_complement = 1 << 4
	Supplementary_map = 1 << 11

	signal = {Abnormal: 0, Normal_foward: 1, Reverse_complement: 2, Supplementary_map:3, Reverse_complement | Supplementary_map:4}
	if Flag in signal:
		return signal[Flag]
	else:
		return 0


def parse_read(read):
	'''
	Check:	1.Flag
			2.Supplementary mapping
			3.Seq
	'''
	local_pos = list()
	process_signal = detect_flag(read.flag) 
	if process_signal == 0:
		return local_pos
		# unmapped read

	pos_start = read.reference_start
	shift = 0
	_shift_read_ = 0
	for element in read.cigar:
		if element[0] == 0 or element[0] == 2:
			shift += element[1]
		if element[0] != 2:
			_shift_read_ += element[1]
		if element[0] in list_flag and element[1] > low_bandary:
			shift += 1
			MEI_contig = read.query_sequence[_shift_read_ - element[1]:_shift_read_]
			# if process_signal == 2 or process_signal == 4:
				# MEI_contig = revcom_complement(MEI_contig)
			# if process_signal == 2 or process_signal == 4:
			# 	# strategy 1:
			# 	read_length = len(read.query_sequence)
			# 	# local_SEQ = read.query_sequence[read_length - _shift_read_:read_length - _shift_read_ + element[1]]
			# 	# MEI_contig = revcom_complement(local_SEQ)
			# 	# strategy 2:
			# 	local_SEQ = revcom_complement(read.query_sequence)
			# 	# MEI_contig = local_SEQ[_shift_read_ - element[1]:_shift_read_]
			# 	MEI_contig = local_SEQ[read_length - _shift_read_:read_length - _shift_read_ + element[1]]
			# else:
			# 	MEI_contig = read.query_sequence[_shift_read_ - element[1]:_shift_read_]
			# MEI_contig = read.query_sequence[_shift_read_-element[1]-4:_shift_read_+10]
			# judge flag !!!!!!!!
			local_pos.append([pos_start + shift, element[1]])
			print read.query_name
			print MEI_contig
	# cluster_pos = sorted(cluster_pos, key = lambda x:x[0])
			# return [r_start + shift, element[1]]
	return local_pos

def out_put(chr, cluster):
	for i in cluster:
		print("%s\t%d\t%d\t%d"%(chr, i[0], i[1], i[2]))

def load_sam(path):
	'''
	Load_BAM_File
	library:	pysam.AlignmentFile
	'''
	samfile = pysam.AlignmentFile(path)
	# print(samfile.get_index_statistics())
	contig_num = len(samfile.get_index_statistics())
	# Acquire_Chr_name
	for _num_ in xrange(contig_num):
		Chr_name = samfile.get_reference_name(_num_)
		cluster_pos = list()
		for read in samfile.fetch(Chr_name):
			feed_back = parse_read(read)

			if len(feed_back) > 0:
				for i in feed_back:
					cluster_pos.append(i)
		cluster_pos = sorted(cluster_pos, key = lambda x:x[0])
		Cluster = cluster(cluster_pos)
		out_put(Chr_name, Cluster)
		break
	samfile.close()