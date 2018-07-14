'''
usage:
	python2 (dbRIP)ALU-L1-SVA > filter.txt
	filter.txt: 1. target sequence duplication
				2. element sequence
'''
from Bio import SeqIO
from pyfasta import Fasta
import random
# import sys
def revcom_complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)[::-1]

def judge_str(seq):
	if seq.find('n') == -1 and seq.find('N') == -1:
		return 1
	else:
		return 0

def acquire_TSD_realELE(seq):
	split_TSD_signal = "<font color=\"green\"><u>"
	split_realELE_signal = "<font color=\"red\">"
	seq = seq.split(split_TSD_signal)
	TSD = "".join(seq[2].split('</u>')[0].split("<br>"))
	realELE = "".join(seq[1].split(split_realELE_signal)[1].split("</font>")[0].split("<br>"))
	return TSD, realELE

def read_dbRIP(Alu_dbRIP_path):
	dbRIP_Alu = list()
	dbRIP_Alu_rc = list()
	file = open(Alu_dbRIP_path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		TSD, realELE = acquire_TSD_realELE(seq[13])
		# seq_2 = seq[13].split(flag_1)
		# for i in seq_2:
		# 	print(i)
		# if TSD.find('n') == -1 and realELE.find('N') == -1:

		# for SVA and L1
		if judge_str(TSD) == 1 and judge_str(realELE) == 1:
			# if len(TSD) <= 25:
			element_seq = TSD+realELE+TSD
			# dbRIP_Alu.append(TSD+"\t"+realELE)
			try:
				dbRIP_Alu_rc.append(revcom_complement(element_seq))
				dbRIP_Alu.append(element_seq)
			except:
				continue
				'''
		for ALU:
		if judge_str(TSD) == 1 and judge_str(realELE) == 1:
			if len(TSD) <= 25:
				print(TSD+"\t"+realELE)
		'''
		# break
	file.close()
	# print("FFFFFFFFFFF%d\t%d"%(len(dbRIP_Alu),len(dbRIP_Alu_rc)))
	return dbRIP_Alu, dbRIP_Alu_rc

def check_length(start, end):
	length = end - start
	if length >= 250 and length <= 400:
		return 1
	else:
		return 0

def check_family(family_name):
	# if family_name[0:2] != "L1" and family_name[0:3] != "Alu" and family_name[0:3] != "SVA" and family_name[0:4] != "HERV":
	# 	return 0
	# else:
	# 	return 1
	if family_name[0:2] == "L1":
		return 2
	elif family_name[0:3] == "Alu":
		return 1
	elif family_name[0:3] == "SVA":
		return 3
	# elif family_name[0:4] == "HERV":
	# 	return 4
	else:
		return 0

def acquire_chr(ref_genome):
	file = Fasta(ref_genome)
	return sorted(file.keys())

def acquire_locus(genome_locus_path, chr_list):
	locus_dic = dict()
	file = open(genome_locus_path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		locus_chr = seq[0]
		start_pos = int(seq[1])
		end_pos = int(seq[2])
		class_ele = seq[3]
		strand = seq[-1]

		flag_1 = check_length(start_pos, end_pos)
		flag_2 = check_family(class_ele)

		# if flag_1 == 0 or flag_2 == 0:
		if flag_2 == 0:
			continue
		# if flag_2 == 1 and flag_1 == 0:
		# 	continue

		if locus_chr in chr_list:
			# statics_class[flag_2 - 1] += 1
			if locus_chr not in locus_dic:
				locus_dic[locus_chr] = list()
			# random int
			# random_locus = random.randint(start_pos, end_pos)
			# locus_dic[locus_chr].append([random_locus, strand, flag_2])
			locus_dic[locus_chr].append([start_pos, end_pos, strand, flag_2])

	# Sum = 0
	for i in locus_dic:
		locus_dic[i] = sorted(locus_dic[i], key=lambda x:x[0])
		# print(i, len(locus_dic[i]))
		# Sum += len(locus_dic[i])
	# print(Sum, statics_class[0], statics_class[1], statics_class[2], statics_class[3])
	file.close()
	return locus_dic

def acquire_locus_random_del(genome_locus_path, chr_list):
	locus_dic = dict()
	locus_dic_ = dict()

	file = open(genome_locus_path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		locus_chr = seq[0]
		start_pos = int(seq[1])
		end_pos = int(seq[2])
		class_ele = seq[3]
		strand = seq[-1]

		flag = check_family(class_ele)
		# if flag == 1:
		if flag == 0:
			if locus_chr not in chr_list:
				continue

			if locus_chr not in locus_dic:
				locus_dic[locus_chr] = list()
				locus_dic_[locus_chr] = list()
			locus_dic[locus_chr].append([start_pos, end_pos, class_ele])
	# random deletion
	for key in locus_dic:
		# kept = int(0.2*len(locus_dic[key]))
		# locus_dic[key] = random.sample(locus_dic[key],kept)
		print("[INFO] Before filtering is %d in %s."%(len(locus_dic[key]), key))
		# locus_dic[key] = random.sample(locus_dic[key],20670)
		locus_dic[key] = sorted(locus_dic[key], key=lambda x:x[0])
		te = 0
		for i in locus_dic[key]:
			if te < i[0] and i[1]-i[0] > 50:
				te = i[1]
				locus_dic_[key].append(i)
			# if len(locus_dic_[key]) == 20000:
			# 	break
		# print len(locus_dic_[key])
		print("[INFO] After filtering is %d in %s."%(len(locus_dic[key]), key))
		if len(locus_dic[key]) < 20000:
			sample_size = len(locus_dic[key])
		else:
			sample_size = 20000
		locus_dic_[key] = random.sample(locus_dic_[key],sample_size)
		locus_dic_[key] = sorted(locus_dic_[key], key=lambda x:x[0])
	file.close()
	return locus_dic_

def acquire_locus_random_del_NonME_chr1(ref_genome, chr_list):
	locus_dic = dict()
	locus_dic_ = dict()

	L_Alu = [random.gauss(280,10) for i in xrange(250000)]
	L_L1 = [random.gauss(6000,50) for i in xrange(50000)]
	L_SVA = [random.gauss(1500,30) for i in xrange(10000)]
	L_total = L_Alu + L_L1 + L_SVA
	random.shuffle(L_total)

	file = Fasta(ref_genome)

	loci_list = random.sample(xrange(len(file[chr_list[0]])), 300000)
	locus_dic[chr_list[0]] = list()
	locus_dic_[chr_list[0]] = list()
	for i in xrange(300000):
		locus_dic[chr_list[0]].append([loci_list[i], loci_list[i]+int(L_total[i]), ""])
	locus_dic[chr_list[0]] = sorted(locus_dic[chr_list[0]], key=lambda x:x[0])
	print("[INFO] Before filtering is %d in %s."%(len(locus_dic[chr_list[0]]), chr_list[0]))
		
	te = 0
	for i in locus_dic[chr_list[0]]:
		if te < i[0]:
			te = i[1]
			locus_dic_[chr_list[0]].append(i)
	print("[INFO] After filtering is %d in %s."%(len(locus_dic_[chr_list[0]]), chr_list[0]))
	if len(locus_dic[chr_list[0]]) < 20000:
		sample_size = len(locus_dic[chr_list[0]])
	else:
		sample_size = 20000
	locus_dic_[chr_list[0]] = random.sample(locus_dic_[chr_list[0]],sample_size)
	locus_dic_[chr_list[0]] = sorted(locus_dic_[chr_list[0]], key=lambda x:x[0])

	return locus_dic_

def acquire_loci(genome_locus_path, chr_list):
	locus_dic = dict()
	_id_ = 0
	file = open(genome_locus_path, 'r')
	for line in file:
		_id_ += 1
		seq = line.strip('\n').split('\t')
		locus_chr = seq[0]
		start_pos = int(seq[1])
		end_pos = int(seq[2])
		class_ele = seq[3]
		strand = seq[-1]

		flag = check_family(class_ele)
		# if flag == 1:
		if flag != 0:
			if locus_chr not in chr_list:
				continue

			if locus_chr not in locus_dic:
				locus_dic[locus_chr] = list()
			locus_dic[locus_chr].append([start_pos, end_pos, class_ele, _id_])
	# random deletion
	for key in locus_dic:
		locus_dic[key] = sorted(locus_dic[key], key=lambda x:x[0])
	file.close()
	return locus_dic

def sum_dict(dic):
	total = 0 
	for i in dic:
		total += len(dic[i])
	return total

def acquire_locus_random(genome_locus_path, chr_list):
	locus_dic = dict()
	file = open(genome_locus_path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		locus_chr = seq[0]
		start_pos = int(seq[1])
		end_pos = int(seq[2])
		class_ele = seq[3]
		strand = seq[-1]

		flag_1 = check_length(start_pos, end_pos)
		flag_2 = check_family(class_ele)

		# if flag_1 == 0 or flag_2 == 0:
		if flag_2 == 0:
			continue
		# if flag_2 == 1 and flag_1 == 0:
		# 	continue

		if locus_chr in chr_list:
			# statics_class[flag_2 - 1] += 1
			if locus_chr not in locus_dic:
				locus_dic[locus_chr] = list()
			# random int
			random_locus = random.randint(start_pos, end_pos)
			locus_dic[locus_chr].append([random_locus, strand, flag_2])
			# locus_dic[locus_chr].append([start_pos, end_pos, strand, flag_2])

	# Sum = 0
	for i in locus_dic:
		locus_dic[i] = sorted(locus_dic[i], key=lambda x:x[0])
		# print(i, len(locus_dic[i]))
		# Sum += len(locus_dic[i])
	# print(Sum, statics_class[0], statics_class[1], statics_class[2], statics_class[3])
	file.close()
	return locus_dic

def load_ref(ref_g):
	return SeqIO.to_dict(SeqIO.parse(ref_g, "fasta"))

def file_out(path, data, Lf):
	file = open(path, 'w')
	for line in data:
		if Lf == 0:
			file.write(line)
		else:
			file.write(('\t'.join(line)) + "\n")
	file.close()

def slice_2_fa(s, standard):
	out_s = str()
	line_id = len(s) / standard
	for i in xrange(line_id):
		out_s = out_s + s[i*standard:(i+1)*standard] + "\n"
	if line_id * standard < len(s):
		out_s = out_s + s[line_id*standard:] + "\n"
	return out_s

# if __name__ == '__main__':
# 	main()