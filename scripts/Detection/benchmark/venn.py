import sys
from data_collection import *

# standard = 20
# dataset_name = ["Tea.alu.bed", "Tea.L1.bed", "Tangram.vcf", "RetroSeq.vcf", "1kg.vcf", "Mobster.txt", "MELT.vcf"]
dataset_name = ["Tea.alu.bed", "Tea.L1.bed", "Tangram.vcf", "1kg.vcf", "Mobster.txt"]
Alu = dict()
Alu_D = dict()
L1 = dict()
L1_D = dict()
SVA = dict()
SVA_D = dict()
Ans = list()
# Ans_tag = ["Tea", "Tangram", "RetroSeq", "1kG", "Mobster", "MELT"]
Ans_tag = ["Tea", "Tangram", "1kG", "Mobster"]

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
	for i in Alu_D:
		Alu_D[i] = sorted(Alu_D[i], key = lambda x:x[0])
		Alu_D[i] = cluster(Alu_D[i])

	# print "L1"
	for i in L1:
		L1[i] = sorted(L1[i], key = lambda x:x[0])
		L1[i] = cluster(L1[i])
		# for j in L1[i]:
			# print i, j[0], j[1], j[2]
	for i in L1_D:
		L1_D[i] = sorted(L1_D[i], key = lambda x:x[0])
		L1_D[i] = cluster(L1_D[i])
	# print "SVA"
	for i in SVA:
		SVA[i] = sorted(SVA[i], key = lambda x:x[0])
		SVA[i] = cluster(SVA[i])
		# for j in SVA[i]:
		# 	print i, j[0], j[1], j[2]
	for i in SVA_D:
		SVA_D[i] = sorted(SVA_D[i], key = lambda x:x[0])
		SVA_D[i] = cluster(SVA_D[i])

def load_data(path):
	# Tea = collect_Tea_plus(path[0], path[1])
	# Tangram = collect_Tangram_plus(path[2])
	# RetroSeq = collect_RetroSeq_plus(path[3])
	# KG = collect_1kg_plus(path[4])
	# Mobster = collect_Mobster_plus(path[5])

	Ans.append(collect_Tea_plus(path[0], path[1]))
	Ans.append(collect_Tangram_plus(path[2]))
	# Ans.append(collect_RetroSeq_plus(path[3]))
	Ans.append(collect_1kg_plus(path[3]))
	Ans.append(collect_Mobster_plus(path[4]))
	# Ans.append(collect_MELT_plus(path[6]))

	for i in xrange(len(Ans)):
		Name = Ans_tag[i]
		dic = Ans[i]
		for chr in dic:
			for element in dic[chr]:
				if element[1] == "AI":
					if chr not in Alu:
						Alu[chr] = list()
					Alu[chr].append([element[0], Name])
				elif element[1] == "AD":
					if chr not in Alu_D:
						Alu_D[chr] = list()
					Alu_D[chr].append([element[0], Name])
				elif element[1] == "LI":
					if chr not in L1:
						L1[chr] = list()
					L1[chr].append([element[0], Name])
				elif element[1] == "LD":
					if chr not in L1_D:
						L1_D[chr] = list()
					L1_D[chr].append([element[0], Name])
				elif element[1] == 'SI':
					if chr not in SVA:
						SVA[chr] = list()
					SVA[chr].append([element[0], Name])
				elif element[1] == "SD":
					if chr not in SVA_D:
						SVA_D[chr] = list()
					SVA_D[chr].append([element[0], Name])
	adjust()

# def compare(chr, pos, subtype):
# 	if subtype == 'AI':
# 		data_ans = Alu
# 		standard = 20
# 	if subtype == 'AD':
# 		data_ans = Alu_D
# 		standard = 20
# 	if subtype == 'LI':
# 		data_ans = L1
# 		standard = 50
# 	if subtype == 'LD':
# 		data_ans = L1_D
# 		standard = 50
# 	if subtype == 'SI':
# 		data_ans = SVA
# 		standard = 50
# 	if subtype == 'SD':
# 		data_ans = SVA_D
# 		standard = 50

# 	if chr in data_ans:
# 		for i in xrange(len(data_ans[chr])):
# 		# for ele in data_ans[chr]:
# 			if data_ans[chr][i][0] - standard <= pos and pos <= data_ans[chr][i][1] + standard:
# 				data_ans[chr][i][3] = 1
# 				return 1
# 			else:
# 				return 0

def compare(chr, pos, subtype):
	if subtype == 'AI':
		data_ans = Alu
		standard = 50
	if subtype == 'AD':
		data_ans = Alu_D
		standard = 50
	if subtype == 'LI':
		data_ans = L1
		standard = 50
	if subtype == 'LD':
		data_ans = L1_D
		standard = 50
	if subtype == 'SI':
		data_ans = SVA
		standard = 50
	if subtype == 'SD':
		data_ans = SVA_D
		standard = 50

	if chr in data_ans:
		for i in xrange(len(data_ans[chr])):
		# for ele in data_ans[chr]:
			if data_ans[chr][i][0] - standard <= pos and pos <= data_ans[chr][i][1] + standard:
				data_ans[chr][i][3] = 1

# def compare_sniffles(chr, pos, subtype):
# 	if subtype == 'XI':
# 		data_ans = [Alu, L1, SVA]
# 	if subtype == 'XD':
# 		data_ans = [Alu_D, L1_D, SVA_D]

# 	for k in xrange(len(data_ans)):
# 		if k == 0:
# 			standard = 20
# 		else:
# 			standard = 50
# 		if chr in data_ans[k]:
# 			for i in xrange(len(data_ans[k][chr])):
# 		# for ele in data_ans[chr]:
# 				if data_ans[k][chr][i][0] - standard <= pos and pos <= data_ans[k][chr][i][1] + standard:
# 					data_ans[k][chr][i][3] = 1
# 					return 1
# 				else:
# 					return 0

def compare_sniffles(chr, pos, subtype):
	if subtype == 'XI':
		data_ans = [Alu, L1, SVA]
	if subtype == 'XD':
		data_ans = [Alu_D, L1_D, SVA_D]

	for k in xrange(len(data_ans)):
		if k == 0:
			standard = 50
		else:
			standard = 50
		if chr in data_ans[k]:
			for i in xrange(len(data_ans[k][chr])):
		# for ele in data_ans[chr]:
				if data_ans[k][chr][i][0] - standard <= pos and pos <= data_ans[k][chr][i][1] + standard:
					data_ans[k][chr][i][3] = 1

def statics():
	Talu_0 = [0]*len(Ans)
	Talu_1 = [0]*len(Ans)
	for key in Alu:
		for ele in Alu[key]:
			if ele[3] == 0:
				Talu_0[len(ele[2])-1] +=1
				# if len(ele[2]) == 6:
				# 	print key, ele
			else:
				# print ele
				Talu_1[len(ele[2])-1] += 1
	for key in Alu_D:
		for ele in Alu_D[key]:
			if ele[3] == 0:
				Talu_0[len(ele[2])-1] +=1
			else:
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

	for key in L1_D:
		for ele in L1_D[key]:
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
	for key in SVA_D:
		for ele in SVA_D[key]:
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

def statics_new():
	test_list = [Alu, Alu_D, L1, L1_D, SVA, SVA_D]
	sta = 0
	A = "Tea"
	B = "Tangram"
	C = "1kG"
	D = "Mobster"
	for i in test_list:
		for chr in i:
			for j in i[chr]:
				# if C in j[2] and A in j[2] and D in j[2] and j[3] == 1:
				# if j[3] == 1 or j[3] == 0:
					# print chr, j
				# if A in j[2] and B in j[2] and  D in j[2] and j[3] == 1:
				if B in j[2] and C in j[2] and D in j[2] and j[3] == 1:
					sta += 1
					# print chr, j	
	# print A,B,D, "tjiang"
	# print A, "rMETL"
	print sta



	# Talu_0 = [0]*len(Ans)
	# Talu_1 = [0]*len(Ans)
	
	# for key in Alu:
	# 	for ele in Alu[key]:
	# 		if ele[3] == 0:
	# 			Talu_0[len(ele[2])-1] +=1
	# 			# if A in ele[2] and B in ele[2] and len(ele[2]) == 2:
	# 			# if A in ele[2] and B in ele[2] and C in ele[2] and len(ele[2]) == 3:
	# 			if A in ele[2] and B in ele[2] and C in ele[2] and D in ele[2] and len(ele[2]) == 4:
	# 				sta += 1
	# 		else:
	# 			Talu_1[len(ele[2])-1] += 1
	# 			# if "Mobster" in ele[2] and len(ele[2]) == 1:
	# 			# if "Tea" in ele[2] and "Mobster" in ele[2] and len(ele[2]) == 2:
	# 			# if "Tangram" in ele[2] and "Mobster" in ele[2] and "Tea" in ele[2] and len(ele[2]) == 3:
	# 			# if "Tea" in ele[2] and "Tangram" in ele[2] and "1kG" in ele[2] and "Mobster" in ele[2] and len(ele[2]) == 4:
	# 				# sta += 1
	# for key in Alu_D:
	# 	for ele in Alu_D[key]:
	# 		if ele[3] == 0:
	# 			Talu_0[len(ele[2])-1] +=1
	# 			# if A in ele[2] and B in ele[2] and len(ele[2]) == 2:
	# 			# if A in ele[2] and B in ele[2] and C in ele[2] and len(ele[2]) == 3:
	# 			if A in ele[2] and B in ele[2] and C in ele[2] and D in ele[2] and len(ele[2]) == 4:
	# 				sta += 1
	# 		else:
	# 			Talu_1[len(ele[2])-1] += 1
	# 			# if "Mobster" in ele[2] and len(ele[2]) == 1:
	# 			# if "Tea" in ele[2] and "Mobster" in ele[2] and len(ele[2]) == 2:
	# 			# if "Tangram" in ele[2] and "Mobster" in ele[2] and "Tea" in ele[2] and len(ele[2]) == 3:
	# 			# if "Tea" in ele[2] and "Tangram" in ele[2] and "1kG" in ele[2] and "Mobster" in ele[2] and len(ele[2]) == 4:
	# 				# sta += 1
	# TL1_0 = [0]*len(Ans)
	# TL1_1 = [0]*len(Ans)
	# for key in L1:
	# 	for ele in L1[key]:
	# 		if ele[3] == 0:
	# 			TL1_0[len(ele[2])-1] +=1
	# 			# if A in ele[2] and B in ele[2] and len(ele[2]) == 2:
	# 			# if A in ele[2] and B in ele[2] and C in ele[2] and len(ele[2]) == 3:
	# 			if A in ele[2] and B in ele[2] and C in ele[2] and D in ele[2] and len(ele[2]) == 4:
	# 				sta += 1
	# 		else:
	# 			TL1_1[len(ele[2])-1] += 1 
	# 			# if "Mobster" in ele[2] and len(ele[2]) == 1:
	# 			# if "Tea" in ele[2] and "Mobster" in ele[2] and len(ele[2]) == 2:
	# 			# if "Tangram" in ele[2] and "Mobster" in ele[2] and "Tea" in ele[2] and len(ele[2]) == 3:
	# 			# if "Tea" in ele[2] and "Tangram" in ele[2] and "1kG" in ele[2] and "Mobster" in ele[2] and len(ele[2]) == 4:
	# 				# sta += 1
	# for key in L1_D:
	# 	for ele in L1_D[key]:
	# 		if ele[3] == 0:
	# 			TL1_0[len(ele[2])-1] +=1
	# 			# if A in ele[2] and B in ele[2] and len(ele[2]) == 2:
	# 			# if A in ele[2] and B in ele[2] and C in ele[2] and len(ele[2]) == 3:
	# 			if A in ele[2] and B in ele[2] and C in ele[2] and D in ele[2] and len(ele[2]) == 4:
	# 				sta += 1
	# 		else:
	# 			TL1_1[len(ele[2])-1] += 1  
	# 			# if "Mobster" in ele[2] and len(ele[2]) == 1:
	# 			# if "Tea" in ele[2] and "Mobster" in ele[2] and len(ele[2]) == 2:
	# 			# if "Tangram" in ele[2] and "Mobster" in ele[2] and "Tea" in ele[2] and len(ele[2]) == 3:
	# 			# if "Tea" in ele[2] and "Tangram" in ele[2] and "1kG" in ele[2] and "Mobster" in ele[2] and len(ele[2]) == 4:
	# 				# sta += 1
	# TSVA_0 = [0]*len(Ans)
	# TSVA_1 = [0]*len(Ans)
	# for key in SVA:
	# 	for ele in SVA[key]:
	# 		if ele[3] == 0:
	# 			TSVA_0[len(ele[2])-1] +=1
	# 			# if A in ele[2] and B in ele[2] and len(ele[2]) == 2:
	# 			# if A in ele[2] and B in ele[2] and C in ele[2] and len(ele[2]) == 3:
	# 			if A in ele[2] and B in ele[2] and C in ele[2] and D in ele[2] and len(ele[2]) == 4:
	# 				sta += 1
	# 		else:
	# 			TSVA_1[len(ele[2])-1] += 1
	# 			# if "Mobster" in ele[2] and len(ele[2]) == 1:
	# 			# if "Tea" in ele[2] and "Mobster" in ele[2] and len(ele[2]) == 2:
	# 			# if "Tangram" in ele[2] and "Mobster" in ele[2] and "Tea" in ele[2] and len(ele[2]) == 3:
	# 			# if "Tea" in ele[2] and "Tangram" in ele[2] and "1kG" in ele[2] and "Mobster" in ele[2] and len(ele[2]) == 4:
	# 				# sta += 1
	# for key in SVA_D:
	# 	for ele in SVA_D[key]:
	# 		if ele[3] == 0:
	# 			TSVA_0[len(ele[2])-1] +=1
	# 			# if A in ele[2] and B in ele[2] and len(ele[2]) == 2:
	# 			# if A in ele[2] and B in ele[2] and C in ele[2] and len(ele[2]) == 3:
	# 			if A in ele[2] and B in ele[2] and C in ele[2] and D in ele[2] and len(ele[2]) == 4:
	# 				sta += 1
	# 		else:
	# 			TSVA_1[len(ele[2])-1] += 1
	# 			# if "Mobster" in ele[2] and len(ele[2]) == 1:
	# 			# if "Tea" in ele[2] and "Mobster" in ele[2] and len(ele[2]) == 2:
	# 			# if "Tangram" in ele[2] and "Mobster" in ele[2] and "Tea" in ele[2] and len(ele[2]) == 3:
	# 			# if "Tea" in ele[2] and "Tangram" in ele[2] and "1kG" in ele[2] and "Mobster" in ele[2] and len(ele[2]) == 4:
	# 				# sta += 1
	# # print "Mobster"+"_"+"tjiang"
	# # print "Tea"+"_"+"Mobster"+"_"+"tjiang"
	# # print "Tangram"+"_"+"Mobster"+"_"+"Tea"+"_"+"tjiang"
	# # print "Tea"+"_"+"Tangram"+"_"+"1kG"+"_"+"Mobster"+"_"+"tjiang"
	# print A, B, C, D
	# print sta

def compare_each_base(chr, breakpoint, subtype):
	if subtype[0] == 'A':
		standard = 50
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

				if Ans[i][chr][j][0] - standard <= breakpoint and Ans[i][chr][j][0] + standard >= breakpoint and Ans[i][chr][j][1] == subtype:
					if subtype[0] == 'A':
						Ans[i][chr][j][2] = 1
					if subtype[0] == 'L':
						Ans[i][chr][j][2] = 2
					if subtype[0] == 'S':
						Ans[i][chr][j][2] = 3

def compare_each_base_sniffles(chr, breakpoint, subtype):

	for i in xrange(len(Ans)):
		if chr in Ans[i]:
			for j in xrange(len(Ans[i][chr])):
				if Ans[i][chr][j][1][0] == 'A' and Ans[i][chr][j][2] == 0:
					Ans[i][chr][j][2] = -1
				if Ans[i][chr][j][1][0] == 'S' and Ans[i][chr][j][2] == 0:
					Ans[i][chr][j][2] = -3
				if Ans[i][chr][j][1][0] == 'L' and Ans[i][chr][j][2] == 0:
					Ans[i][chr][j][2] = -2

				if Ans[i][chr][j][1][0] == 'A':
					standard = 50
				else:
					standard = 50

				if Ans[i][chr][j][0] - standard <= breakpoint and Ans[i][chr][j][0] + standard >= breakpoint and Ans[i][chr][j][1][1] == subtype[1]:
					if Ans[i][chr][j][1][0] == 'A':
						Ans[i][chr][j][2] = 1
					if Ans[i][chr][j][1][0] == 'L':
						Ans[i][chr][j][2] = 2
					if Ans[i][chr][j][1][0] == 'S':
						Ans[i][chr][j][2] = 3

def evaluation(p):
	file = open(p, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		breakpoint = int(seq[1])
		# subtype = seq[3]
		# subtype = subtype.split(':')[2]
		subtype = seq[3][8]+seq[3][1]
		ans = compare(chr, breakpoint, subtype)
		# compare_each_base(chr, breakpoint, subtype)
		# if ans == 0:
		# 	print("%s\t%d\t%s"%(chr, breakpoint, subtype))
	file.close()
	statics_new()

def evaluation_sniffles(p):
	file = open(p, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		breakpoint = int(seq[13])
		# subtype = seq[3]
		# subtype = subtype.split(':')[2]
		# if seq[10][1:4] == "DEL":
		# 	subtype = 'XD'
		# if seq[10][1:4] == "INS" or seq[4][1:4] == "DUP":
		# 	subtype = 'XI'
		# if seq[4][1:4] == "DEL":
		# 	subtype = 'XD'
		# if seq[4][1:4] == "INS" or seq[4][1:4] == "DUP":
		# 	subtype = 'XI'
		if seq[10][:3] == "DEL":
			subtype = 'XD'
		if seq[10][:3] == "INS" or seq[10][:3] == "DUP":
			subtype = 'XI'
		compare_sniffles(chr, breakpoint, subtype)
		# if ans == 0:
		# 	print("%s\t%d\t%s"%(chr, breakpoint, seq[4])) 
		# compare_each_base_sniffles(chr, breakpoint, subtype)
	file.close()
	statics_new()

def evaluation_tag(p):
	file = open(p, 'r')
	for line in file:
		seq = line.strip('\n').split(' ')
		if seq[0][0] == '#':
			continue
		chr = seq[1]
		breakpoint = int(seq[2])
		# subtype = seq[3]
		# subtype = subtype.split(':')[2]
		if seq[0][1:4] == "DEL":
			subtype = 'XD'
		if seq[0][1:4] == "INS":
			subtype = 'XI'
		compare_sniffles(chr, breakpoint, subtype)
		compare_each_base_sniffles(chr, breakpoint, subtype)
	file.close()
	statics()


def main():
	dataset_prefix = sys.argv[1]
	load_path = process_path(dataset_prefix)
	load_data(load_path)
	call_path = sys.argv[2]
	# evaluation(call_path)
	# evaluation_tag(call_path)
	evaluation_sniffles(call_path)



	# sta = 0
	# for i in Alu:
	# 	for j in Alu[i]:
	# 		# if '1kG' in j[2] and 'Mobster' in j[2] and len(j[2]) == 2:
	# 		# if '1kG' in j[2] and 'Tea' in j[2] and "Mobster" in j[2] and len(j[2]) == 3:
	# 		if 'Tea' in j[2] and 'Tangram' in j[2] and '1kG' in j[2] and 'Mobster' in j[2] and len(j[2]) == 4:
	# 			sta += 1
	# for i in Alu_D:
	# 	for j in Alu_D[i]:
	# 		# if '1kG' in j[2] and 'Mobster' in j[2] and len(j[2]) == 2:
	# 		# if '1kG' in j[2] and 'Tea' in j[2] and "Mobster" in j[2] and len(j[2]) == 3:
	# 		if 'Tea' in j[2] and 'Tangram' in j[2] and '1kG' in j[2] and 'Mobster' in j[2] and len(j[2]) == 4:
	# 			sta += 1
	# for i in L1:
	# 	for j in L1[i]:
	# 		# if '1kG' in j[2] and 'Mobster' in j[2] and len(j[2]) == 2:
	# 		# if '1kG' in j[2] and 'Tea' in j[2] and "Mobster" in j[2] and len(j[2]) == 3:
	# 		if 'Tea' in j[2] and 'Tangram' in j[2] and '1kG' in j[2] and 'Mobster' in j[2] and len(j[2]) == 4:
	# 			sta += 1
	# for i in L1_D:
	# 	for j in L1_D[i]:
	# 		# if '1kG' in j[2] and 'Mobster' in j[2] and len(j[2]) == 2:
	# 		# if '1kG' in j[2] and 'Tea' in j[2] and "Mobster" in j[2] and len(j[2]) == 3:
	# 		if 'Tea' in j[2] and 'Tangram' in j[2] and '1kG' in j[2] and 'Mobster' in j[2] and len(j[2]) == 4:
	# 			sta += 1
	# for i in SVA:
	# 	for j in SVA[i]:
	# 		# if '1kG' in j[2] and 'Mobster' in j[2] and len(j[2]) == 2:
	# 		# if '1kG' in j[2] and 'Tea' in j[2] and "Mobster" in j[2] and len(j[2]) == 3:
	# 		if 'Tea' in j[2] and 'Tangram' in j[2] and '1kG' in j[2] and 'Mobster' in j[2] and len(j[2]) == 4:
	# 			sta += 1
	# for i in SVA_D:
	# 	for j in SVA_D[i]:
	# 		# if '1kG' in j[2] and 'Mobster' in j[2] and len(j[2]) == 2:
	# 		# if '1kG' in j[2] and 'Tea' in j[2] and "Mobster" in j[2] and len(j[2]) == 3:
	# 		if 'Tea' in j[2] and 'Tangram' in j[2] and '1kG' in j[2] and 'Mobster' in j[2] and len(j[2]) == 4:
	# 			sta += 1
	# # print "1kG"+"_"+"Mobster"
	# # print "1kG"+"_"+"Tea"+"_"+"Mobster"
	# print "Tea"+"_"+"Tangram"+"_"+"1kG"+"_"+"Mobster"
	# print sta

def cmp(p1, p2):
	d1 = dict()
	file = open(p1, 'r')
	for line in file:
		seq = line.strip().split('\t')
		chr = seq[0]
		breakpoint = int(seq[1])
		type = seq[2][1:4]
		if type == "DEL":
			type = 'D'
		elif type == "INS":
			type = 'I'
		elif type == "DUP":
			type = 'I'
		else:
			continue
		if chr not in d1:
			d1[chr] = list()
		d1[chr].append([breakpoint, type])
	file.close()

	file = open(p2, 'r')
	for line in file:
		flag = 0
		seq = line.strip().split('\t')
		chr = seq[0]
		breakpoint = int(seq[1])
		type = seq[2][1]
		if chr in d1:
			for i in d1[chr]:
				if type == d1[chr][1]:
					if breakpoint - 50 <= d1[chr][0] and d1[chr][0] <= breakpoint + 50:
						flag = 1
		if flag == 0:
			print line,
	file.close()

if __name__ == '__main__':
	main()
	# cmp(sys.argv[1], sys.argv[2])