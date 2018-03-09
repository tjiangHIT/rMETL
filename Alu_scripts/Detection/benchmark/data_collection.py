import sys

def collect_Tea(p1, p2):
	dic_tea = dict()
	tea_alu = open(p1, 'r')
	for line in tea_alu:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		if len(chr) == 3:
			continue
		chr = chr[3:]
		start_pos = int(seq[1])
		end_pos = int(seq[2])
		# length = 
		subtype = seq[3]
		if chr not in dic_tea:
			dic_tea[chr] = list()
		dic_tea[chr].append([start_pos, end_pos, 'alu', 0])
	tea_alu.close()

	tea_L1 = open(p2, 'r')
	for line in tea_L1:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		if len(chr) == 3:
			continue
		chr = chr[3:]
		start_pos = int(seq[1])
		end_pos = int(seq[2])
		# length = 
		subtype = seq[3]
		if chr not in dic_tea:
			dic_tea[chr] = list()
		dic_tea[chr].append([start_pos, end_pos, 'L1', 0])	
	tea_L1.close()
	return dic_tea

def collect_Tangram(p):
	tangram = open(p, 'r')
	dic_tangram = dict()
	for line in tangram:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		start_pos = int(seq[1])
		subtype = seq[4][-3:-1]
		# print seq[7].split(';')[2]
		# break
		# info_len = int(seq[7].split(';')[2].split('=')[1])
		info_len = int(seq[7].split('=')[-1])

		if chr not in dic_tangram:
			dic_tangram[chr] = list()
		dic_tangram[chr].append([start_pos, info_len, subtype, 0])
	tangram.close()
	return dic_tangram

def collect_RetroSeq(p):
	retroseq = open(p, 'r')
	dic_retroseq = dict()
	for line in retroseq:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		pos = int(seq[1])
		subtype = seq[7].split('=')[1].split(',')[0]
		if chr not in dic_retroseq:
			dic_retroseq[chr] = list()
		dic_retroseq[chr].append([pos, pos+1, subtype, 0])
	retroseq.close()
	return dic_retroseq

def collection(p1, p2, p3):
	_benchmark_ = dict()
	tea = open(p1, 'r')
	dic_tea = dict()
	for line in tea:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		if len(chr) == 3:
			continue
		start_pos = int(seq[1])
		end_pos = int(seq[2])
		# length = 
		subtype = seq[3]
		if chr not in dic_tea:
			dic_tea[chr] = list()
		dic_tea[chr].append([start_pos, end_pos])
	tea.close()

	tangram = open(p2, 'r')
	dic_tangram = dict()
	for line in tangram:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		start_pos = int(seq[1])
		subtype = seq[4][-3:-1]
		# print seq[7].split(';')[2]
		# break
		# info_len = int(seq[7].split(';')[2].split('=')[1])
		info_len = int(seq[7].split('=')[-1])

		if subtype == 'AL':
			if chr not in dic_tangram:
				dic_tangram[chr] = list()
			dic_tangram[chr].append([start_pos, info_len])
	tangram.close()

	retroseq = open(p3, 'r')
	dic_retroseq = dict()
	for line in retroseq:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		pos = int(seq[1])
		subtype = seq[7].split('=')[1].split(',')[0]
		if subtype == 'Alu':
			if chr not in dic_retroseq:
				dic_retroseq[chr] = list()
			dic_retroseq[chr].append([pos])
	retroseq.close()
	_benchmark_['Tea'] = dic_tea
	_benchmark_['Tangram'] = dic_tangram
	_benchmark_['RetroSeq'] = dic_retroseq
	return _benchmark_

def test(dic):
	# print len(dic)
	for key in dic:
		# print key
		for k2 in dic[key]:
			# print key, k2	
			for i in xrange(len(dic[key][k2])):
				if len(dic[key][k2][i]) == 1:
					print("%s\t%s\t%d"%(key, k2, dic[key][k2][i][0]))
				else:
					print("%s\t%s\t%d\t%d"%(key, k2, dic[key][k2][i][0], dic[key][k2][i][1]))


def collect():
	Tea_path = sys.argv[1]
	Tangram_path = sys.argv[2]
	RetroSeq_path = sys.argv[3]

	answer = collection(Tea_path, Tangram_path, RetroSeq_path)
	test(answer)


# if __name__ == '__main__':
# 	collect()