# import sys
from collections import Counter

def acquire_count_max(_list_):
	c = Counter(_list_)
	return c.most_common(1)[0]
	# this is a tuple

def revcom_complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)[::-1]

def construct_concensus_info(Ins_list, Clip_list):
	total_count = len(Ins_list) + len(Clip_list)
	if total_count < 5:
		return 0
	breakpoint = list()
	insert_size = list()

	boundary = list()

	for i in Ins_list:
		breakpoint.append(i[0])
		insert_size.append(i[1])
		# boundary.append(i[0])
	for i in Clip_list:
		# boundary.append(i[0])
		if i[2] == 1:
			breakpoint.append(i[0])

	# ==============method_1=====================
	Prob_pos_1 = Counter(breakpoint).most_common(1)[0][0]
	# ==============method_2=====================
	Prob_pos_2 = sum(breakpoint)/len(breakpoint)
	Average_size = int(sum(insert_size)/len(insert_size))

	if Average_size < 50:
		return 0

	# print Average_size
	local_info = list()
	# local_name = "_%d_%d_"%(Prob_pos_2, Average_size)
	local_name = [Prob_pos_2, Average_size]
	local_id = 0
	for i in Ins_list:
		# info = local_name + str(local_id) + '\n' + i[2] + '\n'
		info = local_name + [str(local_id), i[2]]
		# info = local_name + ["%d:%d:%d"%(local_id, min(boundary), max(boundary)), i[2]]
		local_id += 1
		local_info.append(info)
	for i in Clip_list:
		# info = local_name + str(local_id) + '\n' + i[1] + '\n'
		info = local_name + [str(local_id), i[1]]
		# info = local_name + ["%d:%d:%d"%(local_id, min(boundary), max(boundary)), i[1]]
		local_id += 1
		local_info.append(info)
			# print(">%d\n%s"%(i[0], i[2])) 

	# for i in local_info:
	# 	print i

	return local_info


def construct_concensus_seq(Ins_list, Clip_list):
	'''
	Ins_list: 	start position on reference genome
				Insertion size
				Insertion sequence
	Clip_list:	clip position on reference genome
				clip sequence
				clip type(0 for left and 1 for right)
	'''
	breakpoint = list()
	insert_size = list()
	for i in Ins_list:
		breakpoint.append(i[0])
		insert_size.append(i[1])
	for i in Clip_list:
		if i[2] == 1:
			breakpoint.append(i[0])

	# Low_bound = min(breakpoint)
	# Up_bound = max(breakpoint)
	# ==============method_1=====================
	Prob_pos_1 = Counter(breakpoint).most_common(1)[0][0]
	# ==============method_2=====================
	Prob_pos_2 = sum(breakpoint)/len(breakpoint)
	# exist a bug with start position
	# print Prob_pos_2
	# for i in xrange(len(Ins_list)):
	# 	bias = Ins_list[i][0] - Prob_pos_2
	# 	insert_size.append(Ins_list[i][1]+bias)
	# 	print Ins_list[i][0], Ins_list[i][1]+bias
	Max_size = max(insert_size)
	Min_size = min(insert_size)
	Average_size = int(sum(insert_size)/len(insert_size))
	# print Average_size

	Seq = dict()
	for i in Ins_list:
		for j in xrange(i[1]):
			pos = i[0] + j
			ch = i[2][j]
			if pos not in Seq:
				Seq[pos] = list()
			Seq[pos].append(ch)

	for i in Clip_list:
		if Average_size <= len(i[1]):
			boundary = Average_size
		else:
			boundary = len(i[1])

		if i[2] == 0:
			local_clip_seq = i[1][len(i[1])-boundary:]
			for j in xrange(boundary):
				pos = i[0] + j
				ch = local_clip_seq[j]
				if pos not in Seq:
					Seq[pos] = list()
				Seq[pos].append(ch)
		else:
			for j in xrange(boundary):
				pos = i[0] + j
				ch = i[1][j]
				if pos not in Seq:
					Seq[pos] = list()
				Seq[pos].append(ch)

	# concensus has a bug with length also
	Seq_trans = list()
	for key in Seq:
		if len(Seq[key]) < 5:
			continue
		Seq_trans.append([key, acquire_count_max(Seq[key])[0]])
	Seq_trans = sorted(Seq_trans, key = lambda x:x[0])
	final_consensus = str()
	for i in Seq_trans:
		if i[0] < Prob_pos_1:
			continue
		final_consensus += i[1]
		if len(final_consensus) > Average_size:
			break

	# if len(Ins_list) + len(Clip_list) >= 5:
	# 	# print Prob_pos_1, Prob_pos_2, Average_size, len(final_consensus)
	# 	# print Ins_list
	# 	for i in Ins_list:
	# 		print(">%d\n%s"%(i[0], i[2]))
	# 		# print i[0]
	# 		# print i[2]
	# 	# print Clip_list
	# 	for i in Clip_list:
	# 		print(">%d\n%s"%(i[0], i[1]))
			# print i[0]
			# print i[1]
	# print Prob_pos
	# print final_consensus
	# print revcom_complement(final_consensus)

	return final_consensus, Prob_pos_1