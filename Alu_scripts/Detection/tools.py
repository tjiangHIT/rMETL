
import pysam

# list_flag = {1:'I', 4:'S', 5:'H'}
list_flag = {1:'I'}
low_bandary = 20

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

def parse_cigar(read):
	local_pos = list()
	r_start = read.reference_start
	shift = 0
	for element in read.cigar:
		if element[0] == 0 or element[0] == 2:
			shift += element[1]

		if element[0] in list_flag and element[1] > low_bandary:
			shift += 1
			local_pos.append([r_start + shift, element[1]])
			# print(list_flag[element[0]], r_start + shift, element[1])
			# print(chr, read.query_name)
	# cluster_pos = sorted(cluster_pos, key = lambda x:x[0])
			# return [r_start + shift, element[1]]
	return local_pos

def out_put(chr, cluster):
	for i in cluster:
		print chr, i[0], i[1], i[2]

def load_sam(path):
	samfile = pysam.AlignmentFile(path)
	# print(samfile.get_index_statistics())
	contig_num = len(samfile.get_index_statistics())
	for _num_ in xrange(contig_num):
		Chr_name = samfile.get_reference_name(_num_)
		cluster_pos = list()
		for read in samfile.fetch(Chr_name):
			feed_back = parse_cigar(read)
			if len(feed_back) > 0:
				for i in feed_back:
					cluster_pos.append(i)
		cluster_pos = sorted(cluster_pos, key = lambda x:x[0])
		Cluster = cluster(cluster_pos)
		out_put(Chr_name, Cluster)
		break
	samfile.close()