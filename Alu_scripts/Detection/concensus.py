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
		breakpoint.append(i[0])

	# Low_bound = min(breakpoint)
	# Up_bound = max(breakpoint)
	Prob_pos = Counter(breakpoint).most_common(1)[0][0]
	# print Prob_pos
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

	Seq_trans = list()
	for key in Seq:
		if len(Seq[key]) < 5:
			continue
		Seq_trans.append([key, acquire_count_max(Seq[key])[0]])
	Seq_trans = sorted(Seq_trans, key = lambda x:x[0])
	final_consensus = str()
	for i in Seq_trans:
		if i[0] < Prob_pos:
			continue
		final_consensus += i[1]
		if len(final_consensus) > Average_size:
			break

	# print Prob_pos
	# print final_consensus
	# print revcom_complement(final_consensus)

	return final_consensus, Prob_pos