# from pyfasta import Fasta
import sys
import random
import time
from Bio.Seq import Seq
from utils import *

statics_class = [0, 0, 0, 0]

def pipeline(genome_locus_path, ALU_seq_path, ref_genome, L1_seq_path, SVA_seq_path):
	Answer_List = list()
	fake_genome = list()
	chr_list, ref_sequence = acquire_chr(ref_genome)
	Locus_pos = acquire_locus(genome_locus_path, chr_list)
	ALU_list, rc_ALU_list = read_ele_seq(ALU_seq_path)
	print("[INFO]: %d Alu elements will be simulated."%(len(ALU_list)))
	print("[INFO]: Begin to copy-paste Alu elements.")
	_count_ = 0
	_process_ = 0 
	for key in chr_list:
		Answer_pointer = 0

		local_Answer = str()
		if key not in Locus_pos:
			local_Answer = ref_sequence[key][:]
			fake_genome.append(">" + key + "\n" + local_Answer + "\n")
			continue

		for needs in Locus_pos[key]:
			if needs[3] == 1:
				_count_ += 1
				if _count_ == 500:
					_process_ += 1
					print("[INFO]: Round %d, finished %d Alus insertion."%(_process_, _process_*_count_))
					_count_ = 0
				seed = random.randint(0, len(ALU_list) - 1) 
				if needs[2] == "+":
					this_alu = ALU_list[seed]
				else:
					this_alu = rc_ALU_list[seed]
				if Answer_pointer < needs[0]:
					local_Answer = local_Answer + ref_sequence[key][Answer_pointer:needs[0]] + this_alu
					Answer_pointer = needs[1]
					Answer_List.append([key, str(needs[0]), str(needs[1]), this_alu])
		print("[INFO]: Construct %s fake genome."%(key))
		fake_genome.append(">" + key + "\n" + local_Answer + "\n")
	print("[INFO]: Insert %d Alu elements."%(len(Answer_List)))
	return fake_genome, Answer_List

def simulation(Alu_list, rc_Alu_list, Locus_pos, hg19, chr_list, pre_out):
	# fake_genome = SeqIO.SeqRecord("","")
	fake_genome = list()
	# test_num = 0
	_count_ = 0
	_process_ = 0 
	Answer = list()
	for key in chr_list:
		_pointer_ = 0
		_offect_ = 0 
		# _local_fake_genome_ = SeqIO.SeqRecord("","")
		_local_fake_genome_ = list()
		if key not in Locus_pos:
			print("[INFO]: Construct %s fake genome."%(key))
			fake_genome.append(hg19[key])
			# SeqIO.write(hg19[key], pre_out + key + ".fasta", "fasta")
			continue

		for ele in Locus_pos[key]:
			# identify alu
			if ele[3] == 1:
				_count_ += 1
				if _count_ == 10000:
					_process_ += 1
					print("[INFO]: Round %d, finished %d Alus insertion."%(_process_, _process_*_count_))
					_count_ = 0

				seed = random.randint(0, len(Alu_list) - 1) 
				if ele[2] == "+":
					alu_seq = Alu_list[seed]
				else:
					alu_seq = rc_Alu_list[seed]

				if _pointer_ < ele[0]:
					_local_fake_genome_.append(str(hg19[key].seq[_pointer_:ele[0]]))
					_local_fake_genome_.append(alu_seq)
					_pointer_ = ele[1]
					fp_s = ele[0] + _offect_
					# fp_e = fp_s + ele[1] - ele[0]
					_offect_ = _offect_ + len(alu_seq) - ele[1] + ele[0]
					fp_e = fp_s + len(alu_seq)
					Answer.append([key, str(ele[0]), str(ele[1]), str(fp_s), str(fp_e), alu_seq])
					# print(alu_seq)
					# print(Answer)
					# test_num += 1
					# if test_num == 5:
					# 	break
		if _pointer_ < len(hg19[key].seq):
			_local_fake_genome_.append(str(hg19[key].seq[_pointer_:]))
		print("[INFO]: Construct %s fake genome."%(key))
		_convert_local_fake_genome_ = SeqIO.SeqRecord(seq = str(), id = key, name = key, description = key)
		_convert_local_fake_genome_.seq = Seq("".join(_local_fake_genome_))
		# SeqIO.write(_convert_local_fake_genome_, pre_out + key + ".fasta", "fasta")
		fake_genome.append(_convert_local_fake_genome_)
	print("[INFO]: Write fake genome on disk.")
	SeqIO.write(fake_genome, pre_out + "simulation.fa", "fasta")
		# break
	print("[INFO]: Insert %d Alu elements."%(len(Answer)))
	return Answer

def simulation_random(Alu_list, rc_Alu_list, Locus_pos, hg19, chr_list, pre_out):
	# fake_genome = SeqIO.SeqRecord("","")
	fake_genome = list()
	# test_num = 0
	_count_ = 0
	_process_ = 0 
	Answer = list()
	for key in chr_list:
		_pointer_ = 0
		_offect_ = 0 
		# _local_fake_genome_ = SeqIO.SeqRecord("","")
		_local_fake_genome_ = list()
		if key not in Locus_pos:
			print("[INFO]: Construct %s fake genome."%(key))
			fake_genome.append(hg19[key])
			# SeqIO.write(hg19[key], pre_out + key + ".fasta", "fasta")
			continue

		for ele in Locus_pos[key]:
			# identify alu
			if ele[2] == 1:
				_count_ += 1
				if _count_ == 10000:
					_process_ += 1
					print("[INFO]: Round %d, finished %d Alus insertion."%(_process_, _process_*_count_))
					_count_ = 0

				seed = random.randint(0, len(Alu_list) - 1) 
				if ele[1] == "+":
					alu_seq = Alu_list[seed]
				else:
					alu_seq = rc_Alu_list[seed]

				if _pointer_ < ele[0]:
					_local_fake_genome_.append(str(hg19[key].seq[_pointer_:ele[0]]))
					_local_fake_genome_.append(alu_seq)
					_pointer_ = ele[0]
					fp_s = ele[0] + _offect_
					# fp_e = fp_s + ele[1] - ele[0]
					_offect_ = _offect_ + len(alu_seq)
					fp_e = fp_s + len(alu_seq)
					Answer.append([key, str(ele[0]), str(fp_s), str(fp_e), alu_seq])
					# print(alu_seq)
					# print(Answer)
					# test_num += 1
					# if test_num == 5:
					# 	break
		if _pointer_ < len(hg19[key].seq):
			_local_fake_genome_.append(str(hg19[key].seq[_pointer_:]))
		print("[INFO]: Construct %s fake genome."%(key))
		_convert_local_fake_genome_ = SeqIO.SeqRecord(seq = str(), id = key, name = key, description = key)
		_convert_local_fake_genome_.seq = Seq("".join(_local_fake_genome_))
		# SeqIO.write(_convert_local_fake_genome_, pre_out + key + ".fasta", "fasta")
		fake_genome.append(_convert_local_fake_genome_)
	print("[INFO]: Write fake genome on disk.")
	SeqIO.write(fake_genome, pre_out + "simulation.fa", "fasta")
		# break
	print("[INFO]: Insert %d Alu elements."%(len(Answer)))
	return Answer

def simulation_random_del(Locus_pos, hg19, chr_list, pre_out):
	# # fake_genome = SeqIO.SeqRecord("","")
	fake_genome = list()
	# # test_num = 0
	_count_ = 0
	_process_ = 0 
	Answer = list()
	for key in chr_list:
		_pointer_ = 0
		_offect_ = 0 
	# 	# _local_fake_genome_ = SeqIO.SeqRecord("","")
		_local_fake_genome_ = list()
		if key not in Locus_pos:
			print("[INFO]: Construct %s fake genome."%(key))
			fake_genome.append(hg19[key])
			# SeqIO.write(hg19[key], pre_out + key + ".fasta", "fasta")
			continue

		for ele in Locus_pos[key]:
			# identify alu
			_count_ += 1
			if _count_ == 10000:
				_process_ += 1
				print("[INFO]: Round %d, finished %d Alus insertion."%(_process_, _process_*_count_))
				_count_ = 0

	# 			seed = random.randint(0, len(Alu_list) - 1) 
	# 			if ele[1] == "+":
	# 				alu_seq = Alu_list[seed]
	# 			else:
	# 				alu_seq = rc_Alu_list[seed]

			if _pointer_ < ele[0]:
				_local_fake_genome_.append(str(hg19[key].seq[_pointer_:ele[0]]))
				# _local_fake_genome_.append(alu_seq)
				_pointer_ = ele[1]
				fp_s = ele[0] - _offect_
				# fp_e = fp_s + ele[1] - ele[0]
				_offect_ = _offect_ + ele[1] - ele[0]
				# fp_e = fp_s + len(alu_seq)
				Answer.append([key, str(fp_s), str(ele[1] - ele[0]), ele[2]])
					# print(alu_seq)
					# print(Answer)
					# test_num += 1
					# if test_num == 5:
					# 	break
		if _pointer_ < len(hg19[key].seq):
			_local_fake_genome_.append(str(hg19[key].seq[_pointer_:]))
		print("[INFO]: Construct %s fake genome."%(key))
		_convert_local_fake_genome_ = SeqIO.SeqRecord(seq = str(), id = key, name = key, description = key)
		_convert_local_fake_genome_.seq = Seq("".join(_local_fake_genome_))
		# SeqIO.write(_convert_local_fake_genome_, pre_out + key + ".fasta", "fasta")
		fake_genome.append(_convert_local_fake_genome_)
	print("[INFO]: Write fake genome on disk.")
	SeqIO.write(fake_genome, pre_out + "simulation.fa", "fasta")
		# break
	print("[INFO]: Insert %d Alu elements."%(len(Answer)))
	return Answer

def run(ref_g, Alu_l, alu_s, pre_out, choice):
	print("[INFO]: Loading reference genome ...")
	hg19 = load_ref(ref_g)
	if choice == "Trans":
		print("[INFO]: Loading Alu sequences file ...")
		Alu_list, rc_Alu_list = read_dbRIP(alu_s)
		print("[INFO]: Loading Alu locus file ...")
		chr_list = acquire_chr(ref_g)
		Locus_pos = acquire_locus(Alu_l, chr_list)	
		print("[INFO]: %d Alu elements will be simulated."%(len(Alu_list)))
		print("[INFO]: Begin to copy-paste Alu elements.")
		truth = simulation(Alu_list, rc_Alu_list, Locus_pos, hg19, chr_list, pre_out)

	if choice == "Insert":
		print("[INFO]: Loading Alu sequences file ...")
		Alu_list, rc_Alu_list = read_dbRIP(alu_s)
		print("[INFO]: Loading Alu locus file ...")
		chr_list = acquire_chr(ref_g)
		Locus_pos_random = acquire_locus_random(Alu_l, chr_list)	
		print("[INFO]: %d Alu elements will be simulated."%(len(Alu_list)))
		print("[INFO]: Begin to copy-paste Alu elements.")
		truth = simulation_random(Alu_list, rc_Alu_list, Locus_pos_random, hg19, chr_list, pre_out)

	if choice == "random_delete":
		# Delete Alu sequence randomly (10%)
		print("[INFO]: Loading Alu locus file ...")
		chr_list = acquire_chr(ref_g)
		Locus_pos_random_del = acquire_locus_random_del(Alu_l, chr_list)
		total = sum_dict(Locus_pos_random_del)
		# Locus_pos_random = acquire_locus_random(Alu_l, chr_list)	
		print("[INFO]: %d Alu elements will be deleted."%(total))
		print("[INFO]: Begin to delete Alu elements.")
		truth = simulation_random_del(Locus_pos_random_del, hg19, chr_list, pre_out)
		# truth = simulation_random(Alu_list, rc_Alu_list, Locus_pos_random, hg19, chr_list, pre_out)

	print("[INFO]: Write grand truth on disk.")
	# file_out(pre_out+"fake_genome", fake_genome, 0)
	file_out(pre_out+"grand_truth.txt", truth, 1)

def main():
	starttime = time.time()
	ref_g = sys.argv[1]
	Alu_l = sys.argv[2]
	alu_s = sys.argv[3]
	pre_out = sys.argv[4]
	print("[INFO]: 1. Reference Genome: %s"%(ref_g))
	print("[INFO]: 2. Alu locus annotation: %s"%(Alu_l))
	print("[INFO]: 3. Alu sequence annotatio: %s"%(alu_s))
	print("[INFO]: 4. The prefix path of output files: %s"%(pre_out))
	print("[INFO]: This script is uesd for simulate Alus on the reference genome.")
	print("[INFO]: This script will waste your a lot of time. sorry!")
	# run(ref_g, Alu_l, alu_s, pre_out, "Trans")
	run(ref_g, Alu_l, alu_s, pre_out, "random_delete")
	print("[INFO]: Finished in %0.2f seconds."%(time.time() - starttime))
	print("[INFO]: Have a nice day! ^.^ ")

if __name__ == '__main__':
	main()