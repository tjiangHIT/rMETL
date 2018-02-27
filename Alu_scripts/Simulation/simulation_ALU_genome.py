from pyfasta import Fasta
import sys
import random
import time
from Bio.Seq import Seq

statics_class = [0, 0, 0, 0]

def acquire_chr(ref_genome):
	file = Fasta(ref_genome)
	return sorted(file.keys()), file

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
			statics_class[flag_2 - 1] += 1
			if locus_chr not in locus_dic:
				locus_dic[locus_chr] = list()
			# random int
			random_locus = random.randint(start_pos, end_pos)
			# locus_dic[locus_chr].append([random_locus, strand, flag_2])
			locus_dic[locus_chr].append([start_pos, end_pos, strand, flag_2])

	# Sum = 0
	for i in locus_dic:
		locus_dic[i] = sorted(locus_dic[i], key=lambda x:x[0])
		# print(i, len(locus_dic[i]))
		# Sum += len(locus_dic[i])
	# print(Sum, statics_class[0], statics_class[1], statics_class[2], statics_class[3])
	file.close()
	print("[INFO]: Loaded Alu annotations.")
	return locus_dic

def revcom_complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)[::-1]

def slice_2_fa(s, standard):
	out_s = str()
	line_id = len(s) / standard
	for i in xrange(line_id):
		out_s = out_s + s[i*standard:(i+1)*standard] + "\n"
	if line_id * standard < len(s):
		out_s = out_s + s[line_id*standard:] + "\n"
	return out_s

def pipeline(genome_locus_path, ALU_seq_path, ref_genome, L1_seq_path, SVA_seq_path):
	Answer_List = list()
	# Answer_pointer = 0
	fake_genome = list()
	chr_list, ref_sequence = acquire_chr(ref_genome)
	Locus_pos = acquire_locus(genome_locus_path, chr_list)
	ALU_list, rc_ALU_list = read_ele_seq(ALU_seq_path)
	print("[INFO]: %d Alu elements will be simulated."%(len(ALU_list)))
	# SVA_list = read_ele_seq(SVA_seq_path)
	# L1_list = read_ele_seq(L1_seq_path)

	# for i in chr_list: print(i)
	# for j in Locus_pos: print(j)

	print("[INFO]: Begin to copy-paste Alu elements.")
	_count_ = 0
	_process_ = 0 
	for key in chr_list:
		Answer_pointer = 0
		# local_Answer = ">" + key + "\n"
		local_Answer = str()
		if key not in Locus_pos:
			local_Answer = ref_sequence[key][:]
			fake_genome.append(">" + key + "\n" + local_Answer + "\n")
			# fake_genome.append(">" + key + "\n" + slice_2_fa(local_Answer, 60))
			continue

		for needs in Locus_pos[key]:
			if needs[3] == 1:
				_count_ += 1
				if _count_ == 500:
					_process_ += 1
					print("[INFO]: Round %d, finished %d Alus insertion."%(_process_, _process_*_count_))
					_count_ = 0
				# store alu_ID
				# test_pos_1 = time.time()
				seed = random.randint(0, len(ALU_list) - 1) 
				# test_pos_2 = time.time()
				if needs[2] == "+":
					this_alu = ALU_list[seed]
				else:
					this_alu = rc_ALU_list[seed]
					# try:
					# 	this_alu = Seq(ALU_list[seed]).reverse_complement()
					# 	# this_alu = revcom_complement(ALU_list[seed])
					# except:
					# 	continue
					# # except Exception as e:
					# # 	raise e
					# # this_alu = revcom_complement(ALU_list[seed])

				# test_pos_3 = time.time()
				if Answer_pointer < needs[0]:
					local_Answer = local_Answer + ref_sequence[key][Answer_pointer:needs[0]] + this_alu
					Answer_pointer = needs[1]
					Answer_List.append([key, str(needs[0]), str(needs[1]), this_alu])
				# test_pos_4 = time.time()
				# print(key, test_pos_2 - test_pos_1)
				# print(key, test_pos_3 - test_pos_2)
				# print(key, test_pos_4 - test_pos_3)
				# break
		# 	print(key, "a", len(Answer_List))
		# print(key, "b", len(Answer_List))
		print("[INFO]: Construct %s fake genome."%(key))
		fake_genome.append(">" + key + "\n" + local_Answer + "\n")
		# fake_genome.append(">" + key + "\n" + slice_2_fa(local_Answer, 60))
	# print(chr_list)
	print("[INFO]: Insert %d Alu elements."%(len(Answer_List)))
	return fake_genome, Answer_List

def read_ele_seq(path):
	seq_list = list()
	rc_seq_list = list()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		try:
			rc_seq_list.append(revcom_complement(seq[0]+seq[1]+seq[0]))
			seq_list.append(seq[0]+seq[1]+seq[0])
		except:
			continue
		# seq_list.append(seq[0]+seq[1]+seq[0])
		# rc_seq_list.append(Seq(seq[0]+seq[1]+seq[0]).reverse_complement())
	file.close()
	return seq_list, rc_seq_list

def file_out(path, data, Lf):
	file = open(path, 'w')
	for line in data:
		if Lf == 0:
			file.write(line)
		else:
			file.write(('\t'.join(line)) + "\n")
	file.close()

def main():
	genome_locus_path = sys.argv[1]
	ALU_seq_path = sys.argv[3]
	L1_seq_path = sys.argv[4]
	SVA_seq_path = sys.argv[5]
	ref_genome = sys.argv[2]
	output_dir = sys.argv[6]
	print("[INFO]: This script is uesd for simulate Alus in reference genome.")
	print("[INFO]: This script will run for several hours. sorry!")
	print("[INFO]: Loading files ...")
	fake_genome, truth = pipeline(genome_locus_path, ALU_seq_path, ref_genome, L1_seq_path, SVA_seq_path)

	# file_fake = open(, 'w')
	# for chrom in fake_genome:
	# 	file_fake.write(chrom)
	# file_fake.close()
	print("[INFO]: Write fake genome and grand truth on disk.")
	file_out(output_dir+"fake_genome", fake_genome, 0)
	file_out(output_dir+"grand_truth", truth, 1)
	print("[INFO]: Finished.")
	print("[INFO]: Have a nice day! ^.^ ")

if __name__ == '__main__':
	main()