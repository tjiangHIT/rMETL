import sys
from data_collection import *

# standard = 20
dataset_name = ["Tea.alu.bed", "Tea.L1.bed", "Tangram.vcf", "RetroSeq.vcf", "1kg.vcf", "Mobster.txt", "MELT.vcf"]
Alu = dict()
L1 = dict()
SVA = dict()
Ans = list()
Ans_tag = ["Tea", "Tangram", "RetroSeq", "1kG", "Mobster", "MELT"]

def process_path(path):
	load_path = list()
	for name in dataset_name:
		newpath = path + name
		load_path.append(newpath)
	return load_path

def acquire_count_max(_list_):
	c = Counter(_list_)
	return c.most_common(1)[0][0]

def merge_pos(pos_list):
	Name_list = list()
	for ele in pos_list:
		Name_list.append(ele[1])
	start = pos_list[0][0]
	end = pos_list[-1][0]
	result = [start, end, list(set(Name_list)), 0]
	return result

def cluster(pos_list):
	_cluster_ = list()
	temp = list()
	temp.append(pos_list[0])
	for pos in pos_list[1:]:
		# if temp[-1][0] + temp[-1][1] < pos[0]:
		if temp[-1][0] + 50 < pos[0]:
			result = merge_pos(temp)
			# if result != 0:
			_cluster_.append(result)
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	result = merge_pos(temp)
	# if result != 0:
	_cluster_.append(result)

	return _cluster_

def adjust():
	# print "Alu"
	for i in Alu:
		Alu[i] = sorted(Alu[i], key = lambda x:x[0])
		Alu[i] = cluster(Alu[i])
		# for j in Alu[i]:
			# print("chr%s\t%d\t%s"%(i, j[0], j[1]))
			# print i, j[0], j[1], j[2]

	# print "L1"
	for i in L1:
		L1[i] = sorted(L1[i], key = lambda x:x[0])
		L1[i] = cluster(L1[i])
		# for j in L1[i]:
			# print i, j[0], j[1], j[2]

	# print "SVA"
	for i in SVA:
		SVA[i] = sorted(SVA[i], key = lambda x:x[0])
		SVA[i] = cluster(SVA[i])
		# for j in SVA[i]:
		# 	print i, j[0], j[1], j[2]

def load_data(path):
	# Tea = collect_Tea_plus(path[0], path[1])
	# Tangram = collect_Tangram_plus(path[2])
	# RetroSeq = collect_RetroSeq_plus(path[3])
	# KG = collect_1kg_plus(path[4])
	# Mobster = collect_Mobster_plus(path[5])

	Ans.append(collect_Tea_plus(path[0], path[1]))
	Ans.append(collect_Tangram_plus(path[2]))
	Ans.append(collect_RetroSeq_plus(path[3]))
	Ans.append(collect_1kg_plus(path[4]))
	Ans.append(collect_Mobster_plus(path[5]))
	Ans.append(collect_MELT_plus(path[6]))

	for i in xrange(len(Ans)):
		Name = Ans_tag[i]
		dic = Ans[i]
		for chr in dic:
			for element in dic[chr]:
				if element[1][0] == "A":
					if chr not in Alu:
						Alu[chr] = list()
					Alu[chr].append([element[0], Name])
				if element[1][0] == "L":
					if chr not in L1:
						L1[chr] = list()
					L1[chr].append([element[0], Name])
				if element[1][0] == 'S':
					if chr not in SVA:
						SVA[chr] = list()
					SVA[chr].append([element[0], Name])
	adjust()

def compare(chr, pos, subtype):
	if subtype[0] == 'A':
		data_ans = Alu
		standard = 20
	if subtype[0] == 'L':
		data_ans = L1
		standard = 50
	if subtype[0] == 'S':
		data_ans = SVA
		standard = 50

	if chr in data_ans:
		for i in xrange(len(data_ans[chr])):
		# for ele in data_ans[chr]:
			if data_ans[chr][i][0] - standard <= pos and pos <= data_ans[chr][i][1] + standard:
				data_ans[chr][i][3] = 1

def statics():
	Talu_0 = [0]*len(Ans)
	Talu_1 = [0]*len(Ans)
	for key in Alu:
		for ele in Alu[key]:
			if ele[3] == 0:
				Talu_0[len(ele[2])-1] +=1
				# if len(ele[2]) == 5:
				# 	print key, ele
			else:
				# print ele
				Talu_1[len(ele[2])-1] += 1
	TL1_0 = [0]*len(Ans)
	TL1_1 = [0]*len(Ans)
	for key in L1:
		for ele in L1[key]:
			if ele[3] == 0:
				TL1_0[len(ele[2])-1] +=1
				# if len(ele[2]) == 2:
				# 	print key, ele
			else:
				TL1_1[len(ele[2])-1] += 1 
	TSVA_0 = [0]*len(Ans)
	TSVA_1 = [0]*len(Ans)
	for key in SVA:
		for ele in SVA[key]:
			if ele[3] == 0:
				TSVA_0[len(ele[2])-1] +=1
			else:
				TSVA_1[len(ele[2])-1] += 1

	print("\tAlu\tL1\tSVA")
	for i in xrange(len(Ans)):
		print("%d\t%d\t%d\t%d"%(i + 1, Talu_1[i], TL1_1[i], TSVA_1[i]))
	for i in xrange(len(Ans)):
		print("%d\t%d\t%d\t%d"%(i + 1, Talu_0[i], TL1_0[i], TSVA_0[i]))

	for i in xrange(len(Ans)):
		Tr_A = 0
		Tw_A = 0
		Tr_L = 0
		Tw_L = 0
		Tr_S = 0
		Tw_S = 0
		for chr in Ans[i]:
			for j in Ans[i][chr]:
				if j[2] == 1:
					Tr_A += 1
				if j[2] == -1:
					Tw_A += 1
				if j[2] == 2:
					Tr_L += 1
				if j[2] == -2:
					Tw_L += 1
				if j[2] == 3:
					Tr_S += 1
				if j[2] == -3:
					Tw_S += 1
		if Tr_A != 0:
			print("[INFO]: Alu\tFor %s, True is %d, False is %d.\t%0.5f"%(Ans_tag[i], Tr_A, Tw_A, Tr_A*1.0/(Tr_A+Tw_A)))
		else:
			print("[INFO]: Alu\tFor %s, True is %d, False is %d.\t%0.5f"%(Ans_tag[i], Tr_A, Tw_A, 0.0))
		if Tr_L != 0:
			print("[INFO]: L1\tFor %s, True is %d, False is %d.\t%0.5f"%(Ans_tag[i], Tr_L, Tw_L, Tr_L*1.0/(Tr_L+Tw_L)))
		else:
			print("[INFO]: L1\tFor %s, True is %d, False is %d.\t%0.5f"%(Ans_tag[i], Tr_L, Tw_L, 0.0))
		if Tr_S != 0:
			print("[INFO]: SVA\tFor %s, True is %d, False is %d.\t%0.5f"%(Ans_tag[i], Tr_S, Tw_S, Tr_S*1.0/(Tr_S+Tw_S)))
		else:
			print("[INFO]: SVA\tFor %s, True is %d, False is %d.\t%0.5f"%(Ans_tag[i], Tr_S, Tw_S, 0.0))

def compare_each_base(chr, breakpoint, subtype):
	if subtype[0] == 'A':
		standard = 20
	else:
		standard = 50
	for i in xrange(len(Ans)):
		if chr in Ans[i]:
			for j in xrange(len(Ans[i][chr])):
				if Ans[i][chr][j][1][0] == 'A' and Ans[i][chr][j][2] == 0:
					Ans[i][chr][j][2] = -1
				if Ans[i][chr][j][1][0] == 'S' and Ans[i][chr][j][2] == 0:
					Ans[i][chr][j][2] = -3
				if Ans[i][chr][j][1][0] == 'L' and Ans[i][chr][j][2] == 0:
					Ans[i][chr][j][2] = -2

				if Ans[i][chr][j][0] - standard <= breakpoint and Ans[i][chr][j][0] + standard >= breakpoint and Ans[i][chr][j][1][0] == subtype[0]:
					if subtype[0] == 'A':
						Ans[i][chr][j][2] = 1
					if subtype[0] == 'L':
						Ans[i][chr][j][2] = 2
					if subtype[0] == 'S':
						Ans[i][chr][j][2] = 3

def evaluation(p):
	file = open(p, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		breakpoint = int(seq[1])
		subtype = seq[3]
		subtype = subtype.split(':')[2]
		compare(chr, breakpoint, subtype)
		compare_each_base(chr, breakpoint, subtype)
	file.close()
	statics()


def main():
	dataset_prefix = sys.argv[1]
	load_path = process_path(dataset_prefix)
	load_data(load_path)
	call_path = sys.argv[2]
	evaluation(call_path)


if __name__ == '__main__':
	main()