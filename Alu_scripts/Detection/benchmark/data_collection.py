import sys

def parse_Tea_subtype(sub_list):
	seq = sub_list.split(',')
	out_sub_list = list()
	for i in xrange(len(seq)):
		local_subtype = seq[i]
		if len(local_subtype.split('_')) == 2:
			if local_subtype.split('_')[1][0] == 'A':
				out_sub_list.append(local_subtype.split('_')[1])
		else:
			out_sub_list.append(local_subtype)
	return out_sub_list

def collect_Tea_plus(p1, p2):
	dic_tea = dict()
	alu_Tea = open(p1, 'r')
	for line in alu_Tea:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		if len(chr) == 3:
			continue
		chr = chr[3:]
		breakpoint = int(seq[1])
		# subtype = parse_Tea_subtype(seq[3])
		subtype = "Alu"
		if chr not in dic_tea:
			dic_tea[chr] = list()
		dic_tea[chr].append([breakpoint, subtype, 0])
	alu_Tea.close()

	L1_Tea = open(p2, 'r')
	for line in L1_Tea:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		if len(chr) == 3:
			continue
		chr = chr[3:]
		breakpoint = int(seq[1])
		# subtype = parse_Tea_subtype(seq[3])
		subtype = "L1"
		if chr not in dic_tea:
			dic_tea[chr] = list()
		dic_tea[chr].append([breakpoint, subtype, 0])
	L1_Tea.close()
	return dic_tea

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

def collect_Tangram_plus(p):
	tangram = open(p, 'r')
	dic_tangram = dict()
	for line in tangram:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		breakpoint = int(seq[1])
		subtype = seq[4][-3:-1]
		if subtype[0] == 'A':
			subtype = "Alu"
		elif subtype[0] == 'L':
			subtype = "L1"
		elif subtype[0] == "S":
			subtype = "SVA"
		else:
			subtype = "HERV"
		# strand = seq[7].split()
		# genotype = seq[9]

		if chr not in dic_tangram:
			dic_tangram[chr] = list()
		dic_tangram[chr].append([breakpoint, subtype, 0])
	tangram.close()
	return dic_tangram

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

def collect_RetroSeq_plus(p):
	retroseq = open(p, 'r')
	dic_retroseq = dict()
	for line in retroseq:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		breakpoint = int(seq[1])
		subtype = seq[7].split('=')[1].split(',')[0]
		if subtype[0] == 'A':
			subtype = "Alu"
		if subtype[0] == "L":
			subtype = "L1"
		if chr not in dic_retroseq:
			dic_retroseq[chr] = list()
		dic_retroseq[chr].append([breakpoint, subtype, 0])
	retroseq.close()
	return dic_retroseq

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

def parse_genotype(genotype):
	left = genotype.split(':')[0].split('/')[0]
	if left == '.':
		left = 0
	else:
		left = int(left)
	right = genotype.split(':')[0].split('/')[1]
	if right == '.':
		right = 0
	else:
		right = int(right)

	if left + right > 0:
		return 1
	else:
		return 0

def parse_info(info):
	svtype = info.split(';')[3].split('=')[1]
	svlen = int(info.split(';')[4].split('=')[1])
	subtype = info.split(';')[5].split('=')[1].split(',')[0]
	polarity = info.split(';')[5].split('=')[1].split(',')[3]
	return [svlen, subtype, polarity]

def collect_1000G(p):
	NA12878_ans = list()
	file = open(p, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		pos = int(seq[1])
		info = parse_info(seq[7])

		genotype = parse_genotype(seq[432])
		if genotype == 0:
			continue
		else:
			NA12878_ans.append([chr, pos, info])
	file.close()
	return NA12878_ans

def collect_1kg(p):
	dic_1kg = dict()
	local_list = collect_1000G(p)
	for i in local_list:
		if i[0] not in dic_1kg:
			dic_1kg[i[0]] = list()
		dic_1kg[i[0]].append([i[1], 0])
	return dic_1kg

def collect_1kg_plus(p):
	dic_1kg = dict()
	file = open(p, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		breakpoint = int(seq[1])
		# subtype = seq[7].split(';')[3].split('=')[1]
		subtype = seq[4].split(':')[2][0]+seq[4].split(':')[0][1]
		# if subtype[0] == "A":
		# 	subtype = "Alu"
		# if subtype[0] == "L":
		# 	subtype = "L1"
		# if subtype[0] == "S":
		# 	subtype = "SVA"
		genotype = parse_genotype(seq[432])
		if genotype == 1:
			if chr not in dic_1kg:
				dic_1kg[chr] = list()
			dic_1kg[chr].append([breakpoint, subtype, 0])
	file.close()
	return dic_1kg

def collect_Mobster_plus(p):
	dic_Mobster = dict()
	file = open(p, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0][3:]
		breakpoint = int(seq[4])
		subtype = seq[3]
		if subtype[0] == "A":
			subtype = "Alu"
		if subtype[0] == "L":
			subtype = "L1"
		if subtype[0] == "S":
			subtype = "SVA"
		if chr not in dic_Mobster:
			dic_Mobster[chr] = list()
		dic_Mobster[chr].append([breakpoint, subtype, 0])
	file.close()
	return dic_Mobster

def collect_MELT_plus(p):
	dic_MELT = dict()
	file = open(p, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue

		chr = seq[0]
		breakpoint = int(seq[1])
		subtype = seq[4].split(':')[2][0]+seq[4].split(':')[0][1]
		if chr not in dic_MELT:
			dic_MELT[chr] = list()
		dic_MELT[chr].append([breakpoint,subtype, 0])
	file.close()
	return dic_MELT
