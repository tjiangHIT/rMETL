
import pysam
import Queue as Q
from concensus import *
import cigar
from Bio import SeqIO
from Bio.Seq import Seq
from genotype import *
from CommandRunner import *
from multiprocessing import Pool
import os
import argparse
import logging
import tempfile
import sys, time

# list_flag = {1:'I', 4:'S', 5:'H'}
INS_flag = {1:'I'}
DEL_flag = {2:'D'}
clip_flag = {4:'S', 5:'H'}
low_bandary = 20
P_homozygous = 0.8
P_heterozygous = 0.3

CLIP_note = dict()
total_signal = list()
global_ref = list()
# clip_store = Q.PriorityQueue()

def revcom_complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'T', 'c': 'G', 'g': 'C', 't': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)[::-1]

def detect_flag(Flag):
	# Signal
	Normal_foward = 1 >> 1
	Abnormal = 1 << 2
	Reverse_complement = 1 << 4
	Supplementary_map = 1 << 11

	signal = {Abnormal: 0, Normal_foward: 1, Reverse_complement: 2, Supplementary_map:3, Reverse_complement | Supplementary_map:4}
	if Flag in signal:
		return signal[Flag]
	else:
		return 0

# def acquire_ins_part():
# 	# about insert sequence

def store_clip_pos(locus, chr, seq, flag):
	# about collecting breakpoint from clipping 
	hash_1 = int(locus /10000)
	mod = locus % 10000
	hash_2 = int(mod / 50)
	element = [locus, seq, flag]

	if hash_1 not in CLIP_note[chr]:
		CLIP_note[chr][hash_1] = dict()
		CLIP_note[chr][hash_1][hash_2] = list()
		CLIP_note[chr][hash_1][hash_2].append(element)
	else:
		if hash_2 not in CLIP_note[chr][hash_1]:
			CLIP_note[chr][hash_1][hash_2] = list()
			CLIP_note[chr][hash_1][hash_2].append(element)
		else:
			CLIP_note[chr][hash_1][hash_2].append(element)

def acquire_clip_pos(deal_cigar):
	seq = list(cigar.Cigar(deal_cigar).items())
	if seq[0][1] == 'S':
		first_pos = seq[0][0]
	else:
		first_pos = 0
	if seq[-1][1] == 'S':
		last_pos = seq[-1][0]
	else:
		last_pos = 0
	# seq = cigar.split('S')
	# if len(seq) == 3:
	# 	first_pos = int(seq[0])
	# 	last_pos = int(seq[1].split('M')[-1])
	# 	return [first_pos, last_pos]
	# if len(seq) == 1:
	# 	return []
	# if len(seq) == 2:
	# 	if seq[1] == '':
	# 		return []
	# 	first_pos = int(seq[0])
	# 	last_pos = 0
	bias = 0
	for i in seq:
		if i[1] == 'M' or i[1] == 'D':
			bias += i[0]

	return [first_pos, last_pos, bias]

def organize_split_signal(chr, primary_info, Supplementary_info, total_L):
	overlap = list()
	for i in Supplementary_info:
		seq = i.split(',')
		local_chr = seq[0]
		local_start = int(seq[1])
		local_cigar = seq[3]

		dic_starnd = {1:'+', 2: '-'}

		if dic_starnd[primary_info[4]] != seq[2]:
			continue
		if chr != local_chr:
			# return overlap
			continue
		local_set = acquire_clip_pos(local_cigar)
		# if len(local_set) == 0:
		# 	continue
		if primary_info[0] < local_start:
			if primary_info[3]+local_set[0]-total_L > 20:
				overlap.append([total_L - primary_info[3], local_set[0], primary_info[1]])
		else:
			if local_set[1]+primary_info[2]-total_L > 20:
				overlap.append([total_L - local_set[1], primary_info[2], local_start+local_set[2]-1])
			# exist some bugs

		# if local_start <= primary_end + 50 and local_start >= primary_end - 50:
		# 	local_set = acquire_clip_pos(local_cigar)
		# 	if len(local_set) != 0:
		# 		overlap.append([primary_clip, sum(local_set) - primary_clip])
	return overlap



def parse_read(read, Chr_name):
	'''
	Check:	1.Flag
			2.Supplementary mapping
			3.Seq
	'''
	DEL_ME_pos = list()
	INS_ME_pos = list()
	process_signal = detect_flag(read.flag)
	if process_signal == 0:
		return INS_ME_pos, DEL_ME_pos

	# Add DEL:ME type call signal
	# DEL_ME_pos = list()
	pos_start = read.reference_start
	shift = 0
	# temp_test = dict()
	for element in read.cigar:
		# if element[0] not in temp_test:
		# 	temp_test[element[0]] = 0
		# temp_test[element[0]] += element[1]
		if element[0] == 0:
			shift += element[1]
		# if element[0] == 1:
		# 	shift += 1
		if element[0] in DEL_flag and element[1] <= low_bandary:
			shift += element[1]
		if element[0] in DEL_flag and element[1] > low_bandary:

			# +++++++++++++++++++++++BUG++++++++++++++++++++
			# if 'chr'+Chr_name in Ref:
			# 	MEI_contig = str(Ref['chr'+Chr_name].seq[pos_start+shift:shift+element[1]+pos_start])
			# Normal condation
			# MEI_contig = str(Ref[Chr_name].seq[pos_start+shift:shift+element[1]+pos_start])
			# ++++++++++++++++++++++++++++++++++++++++++++++

			# MEI_contig = read.query_alignment_sequence[shift:shift+element[1]]
			DEL_ME_pos.append([pos_start+shift, element[1]])
			# print read.query_name, read.reference_length, len(read.query_alignment_sequence)
			# print temp_test, shift
			# print read.query_name, pos_start, pos_start+shift, element[1], MEI_contig
			shift += element[1]

	# Add INS:ME type call signal
	INS_ME_pos = list()
	# process_signal = detect_flag(read.flag) 
	# if process_signal == 0:
	# 	return local_pos
		# unmapped read

	pos_start = read.reference_start
	shift = 0
	_shift_read_ = 0
	pos_end = read.reference_end
	primary_clip_0 = 0
	primary_clip_1 = 0
	for element in read.cigar:
		if element[0] == 0 or element[0] == 2:
			shift += element[1]
		if element[0] != 2:
			_shift_read_ += element[1]
		if element[0] in INS_flag and element[1] > low_bandary:
			shift += 1
			MEI_contig = read.query_sequence[_shift_read_ - element[1]:_shift_read_]
			# if process_signal == 2 or process_signal == 4:
				# MEI_contig = revcom_complement(MEI_contig)
			# if process_signal == 2 or process_signal == 4:
			# 	# strategy 1:
			# 	read_length = len(read.query_sequence)
			# 	# local_SEQ = read.query_sequence[read_length - _shift_read_:read_length - _shift_read_ + element[1]]
			# 	# MEI_contig = revcom_complement(local_SEQ)
			# 	# strategy 2:
			# 	local_SEQ = revcom_complement(read.query_sequence)
			# 	# MEI_contig = local_SEQ[_shift_read_ - element[1]:_shift_read_]
			# 	MEI_contig = local_SEQ[read_length - _shift_read_:read_length - _shift_read_ + element[1]]
			# else:
			# 	MEI_contig = read.query_sequence[_shift_read_ - element[1]:_shift_read_]
			# MEI_contig = read.query_sequence[_shift_read_-element[1]-4:_shift_read_+10]
			# judge flag !!!!!!!!
			INS_ME_pos.append([pos_start + shift, element[1], MEI_contig])

			# print read.query_name, "I", pos_start + shift
			# print MEI_contig

		if element[0] in clip_flag:

			if shift == 0:
				primary_clip_0 = element[1]
			else:
				primary_clip_1 = element[1]

			if element[1] > low_bandary:
				if shift == 0:
					clip_pos = pos_start - 1
					clip_contig = read.query_sequence[:element[1]]
					store_clip_pos(clip_pos, Chr_name, clip_contig, 0)

					# primary_clip_0 = element[1]
					# left clip size

				else:
					clip_pos = pos_start + shift - 1
					# primary_clip = read.query_length - element[1]
					clip_contig = read.query_sequence[read.query_length - element[1]:]
					store_clip_pos(clip_pos, Chr_name, clip_contig, 1)

					# primary_clip_1 = read.query_length - element[1]
					# right clip size


				# store_clip_pos(clip_pos, Chr_name, clip_contig)
				# print read.query_name, "S"
				# print clip_contig
				# hash_1 = int(clip_pos / 10000)
				# mod_1 = int(clip_pos % 10000)
				# hash_2 = int(mod_1 / 50)
				# if hash_1 not in CLIP_note[Chr_name]:
				# 	CLIP_note[Chr_name][hash_1] = dict()
				# 	CLIP_note[Chr_name][hash_1][hash_2] = list()
				# 	CLIP_note[Chr_name][hash_1][hash_2].append(clip_pos)
				# else:
				# 	if hash_2 not in CLIP_note[Chr_name][hash_1]:
				# 		CLIP_note[Chr_name][hash_1][hash_2]
				# # if clip_pos not in CLIP_note[Chr_name]:
				# # 	CLIP_note[Chr_name][clip_pos] = 1
				# # else:
				# # 	CLIP_note[Chr_name][clip_pos] += 1
				# CLIP_note[Chr_name].put(clip_pos)
				# print Chr_name, clip_pos
	# cluster_pos = sorted(cluster_pos, key = lambda x:x[0])
			# return [r_start + shift, element[1]]

	if process_signal == 1 or process_signal == 2:
		Tags = read.get_tags()
		chr = Chr_name
		# primary_clip = pos_start
		primary_info = [pos_start, pos_end, primary_clip_0, primary_clip_1, process_signal]

		for i in Tags:
			if i[0] == 'SA':
				Supplementary_info = i[1].split(';')[:-1]
				# print process_signal
				# print chr, primary_info, read.query_length
				# print read.cigar
				# print i[1].split(';')[-1]
				overlap = organize_split_signal(chr, primary_info, Supplementary_info, read.query_length)
				for k in overlap:
					# print k
					MEI_contig = read.query_sequence[k[0]:k[1]]
					INS_ME_pos.append([k[2], k[1] - k[0], MEI_contig])

	return INS_ME_pos, DEL_ME_pos

def merge_siganl(chr, cluster):
	# for i in cluster:
	# 	if i[2] >= 5:
	# 		total_signal.append("%s\t%d\t%d\t%d\t%s\n"%(chr, i[0], i[1], i[2], i[3]))
	# 		# print("%s\t%d\t%d\t%d"%(chr, i[0], i[1], i[2]))
	for i in cluster:
		for j in i:
			total_signal.append(j)

def acquire_clip_locus(down, up, chr):
	list_clip = list()
	if int(up/10000) == int(down/10000):
		key_1 = int(down/10000)
		if key_1 not in CLIP_note[chr]:
			return list_clip
		for i in xrange(int((up%10000)/50)-int((down%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = int((down%10000)/50)+i
			if key_2 not in CLIP_note[chr][key_1]:
				continue
			for ele in CLIP_note[chr][key_1][key_2]:
				if ele[0] >= down and ele[0] <= up:
					list_clip.append(ele)
	else:
		key_1 = int(down/10000)
		if key_1 in CLIP_note[chr]:
			for i in xrange(200-int((down%10000)/50)):
				# exist a bug ***********************************
				key_2 = int((down%10000)/50)+i
				if key_2 not in CLIP_note[chr][key_1]:
					continue
				for ele in CLIP_note[chr][key_1][key_2]:
					if ele[0] >= down and ele[0] <= up:
						list_clip.append(ele)
		key_1 += 1
		if key_1 not in CLIP_note[chr]:
			return list_clip
		for i in xrange(int((up%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = i
			if key_2 not in CLIP_note[chr][key_1]:
				continue
			for ele in CLIP_note[chr][key_1][key_2]:
				if ele[0] >= down and ele[0] <= up:
					list_clip.append(ele)
	return list_clip

def merge_pos(pos_list, chr):
	start = list()
	end = list()
	for ele in pos_list:
		start.append(ele[0])
		end.append(ele[0] + ele[1])

	search_down = min(start) - 10
	search_up = max(start) + 10
	temp_clip = acquire_clip_locus(search_down, search_up, chr)

	# concensus, ref_pos = construct_concensus_seq(pos_list, temp_clip)
	result = construct_concensus_info(pos_list, temp_clip)
	if result != 0:
		for i in xrange(len(result)):
			# result[i] = ">INS_" + chr + result[i]
			result[i] = ["INS", chr] + result[i] + [len(result)]
		return result
	else:
		return 0

	# print ref_pos
	# print concensus

	'''
	TE detect
	'''

	# total_breakpoint = int((sum(start) + sum(temp_clip)) / (len(pos_list) + len(temp_clip)))
	# total_length = int(sum(end)/len(pos_list)) - total_breakpoint
	# total_read_count = len(pos_list) + len(temp_clip)
	# return	[ref_pos, len(concensus), total_read_count, concensus]
	# return [int(sum(start)/len(pos_list)), int(sum(end)/len(pos_list)) - int(sum(start)/len(pos_list)), len(pos_list)]

def cluster(pos_list, chr):
	_cluster_ = list()
	temp = list()
	temp.append(pos_list[0])
	for pos in pos_list[1:]:
		# if temp[-1][0] + temp[-1][1] < pos[0]:
		if temp[-1][0] + 20 < pos[0]:
			result = merge_pos(temp, chr)
			if result != 0:
				_cluster_.append(result)
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	result = merge_pos(temp, chr)
	if result != 0:
		_cluster_.append(result)
	# _cluster_.append(merge_pos(temp, chr))
	return _cluster_

def merge_pos_del(pos_list, chr, Ref):
	start = list()
	end = list()
	for ele in pos_list:
		start.append(ele[0])
		end.append(ele[0] + ele[1])

	breakpoint = sum(start)/len(start)
	size = sum(end)/len(end) - breakpoint

	result = list()
	if len(pos_list) < 5:
		return result
	else:
		if chr in Ref and size >= 50:
		# if 'chr'+chr in Ref and size >= 50:
			# result.append(">DEL_%s_%d_%d_%d\n%s\n"%(chr, breakpoint, size, len(pos_list), str(Ref['chr'+chr].seq[breakpoint:breakpoint+size])))
			result.append(['DEL', chr, breakpoint, size, len(pos_list), str(Ref[chr].seq[breakpoint:breakpoint+size])])
			# result.append(['DEL', chr, breakpoint, size, len(pos_list), str(Ref['chr'+chr].seq[breakpoint:breakpoint+size])])
	# for i in xrange(len(pos_list)):
	# 	result.append(">DEL_%s_%d_%d_%d\n%s\n"%(chr, breakpoint, size, i, pos_list[i][2]))
	return result


def cluster_del(pos_list, chr, Ref):
	_cluster_ = list()
	temp = list()
	temp.append(pos_list[0])
	for pos in pos_list[1:]:
		if temp[-1][0] + 20 < pos[0]:
			result = merge_pos_del(temp, chr, Ref)
			if len(result) != 0:
				_cluster_.append(result)
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	result = merge_pos_del(temp, chr, Ref)
	if len(result) != 0:
		_cluster_.append(result)
	return _cluster_

def load_ref(ref_g):
	logging.info("Loading reference genome...")
	return SeqIO.to_dict(SeqIO.parse(ref_g, "fasta"))

def combine_result(INS, DEL):
	result = list()
	# INS_chr_pos_len_#_seq_rc_dp
	for i in INS:
		for j in i:
			if len(j) != 8:
				continue
			key = "%s_%s_%d_%d_%s_%d_%d"%(j[0], j[1], j[2], j[3], j[4], j[6], j[7])
			fake_seq = SeqIO.SeqRecord(seq = str(), id = key, name = key, description = key)
			fake_seq.seq = Seq(j[5])
			result.append(fake_seq)
	# DEL_chr_pos_len_seq_rc_dp
	for i in DEL:
		for j in i:
			if len(j) != 7:
				continue
			key = "%s_%s_%d_%d_%d_%d"%(j[0], j[1], j[2], j[3], j[4], j[6])
			fake_seq = SeqIO.SeqRecord(seq = str(), id = key, name = key, description = key)
			fake_seq.seq = Seq(j[5])
			result.append(fake_seq)
	# print INS+DEL
	# Temp = sorted(INS + DEL, key = lambda x:x[2])
	# result = list()
	# for i in Temp:
	# 	key = "%s_%s_%d_%d_%s"%(i[0], i[1], i[2], i[3], i[4])
	# 	fake_seq = SeqIO.SeqRecord(seq = str(), id = key, name = key, description = key)
	# 	fake_seq.seq = Seq(i[5])
	# 	result.append(fake_seq)
	return result
	# return INS+DEL

def count_coverage(chr, s, e, f):
	total = 0
	for i in f.fetch(chr, s, e):
		total += 1
	return total

def add_genotype(info_list, file):
	for i in xrange(len(info_list)):
		if info_list[i][0][0] == 'INS':
			chr = info_list[i][0][1]
			'''
			method_1
			'''
			start = info_list[i][0][2]-20
			# end = info_list[i][0][2] + info_list[i][0][3]
			end = info_list[i][0][2] + 20
			# evidence = len(info_list[i])
			# '''
			# method_2
			# '''
			# start = int(info_list[i][0][4].split(':')[1])
			# end = int(info_list[i][0][4].split(':')[2])
			locus_cov = count_coverage(chr, start, end, file)
			# # GT, GL = caculate_genotype_likelyhood(evidence, locus_cov)
			# result = simple_call_genotype(evidence, locus_cov, P_heterozygous, P_homozygous)
			# if result != 0:
			# 	# GT, GL = simple_call_genotype(evidence, locus_cov)
			# 	GT, GL = result[0], result[1]
			# 	for j in xrange(len(info_list[i])):
			# 		info_list[i][j].append(GT)
			# 		info_list[i][j].append(GL)
			for j in xrange(len(info_list[i])):
				info_list[i][j].append(locus_cov)

		else:
			# DEL_chr_pos_len_seq_rc_dp
			for j in xrange(len(info_list[i])):
			# if info_list[i][j][0] == 'DEL':
				chr = info_list[i][j][1]
				start = info_list[i][j][2]
				end = info_list[i][j][2]+info_list[i][j][3]
				# evidence = info_list[i][j][4]
				locus_cov = count_coverage(chr, start, end, file)
				# # GT, GL = caculate_genotype_likelyhood(evidence, locus_cov)
				# # GT, GL = simple_call_genotype(evidence, locus_cov)
				# result = simple_call_genotype(evidence, locus_cov, P_heterozygous, P_homozygous)
				# if result != 0:
				# 	GT, GL = result[0], result[1]
				# 	info_list[i][j].append(GT)
				# 	info_list[i][j].append(GL)
				info_list[i][j].append(locus_cov)
	return info_list

def load_sam(args):
	'''
	Load_BAM_File
	library:	pysam.AlignmentFile

	load_Ref_Genome
	library:	Bio
	'''
	p1 = args.AlignmentFile
	p2 = args.Output_prefix
	p3 = args.Reference
	# global low_bandary
	# low_bandary = args.min_distance
	# global P_homozygous
	# P_homozygous = args.homozygous
	# global P_heterozygous
	# P_heterozygous = args.heterozygous


	Ref = load_ref(p3)

	samfile = pysam.AlignmentFile(p1)
	# print(samfile.get_index_statistics())
	contig_num = len(samfile.get_index_statistics())
	logging.info("The total number of chromsomes: %d"%(contig_num))
	# Acquire_Chr_name
	for _num_ in xrange(contig_num):
		Chr_name = samfile.get_reference_name(_num_)
		logging.info("Resolving the chromsome %s."%(Chr_name))
		# Chr_length = samfile.lengths[_num_]
		if Chr_name not in CLIP_note:
			# CLIP_note[Chr_name] = [0] * Chr_length
			# CLIP_note[Chr_name] = Q.PriorityQueue()
			CLIP_note[Chr_name] = dict()

		cluster_pos_INS = list()
		cluster_pos_DEL = list()
		for read in samfile.fetch(Chr_name):
			feed_back, feed_back_del = parse_read(read, Chr_name)

			if len(feed_back) > 0:
				for i in feed_back:
					cluster_pos_INS.append(i)

			if len(feed_back_del) > 0:
				for i in feed_back_del:
					cluster_pos_DEL.append(i)
		# while not CLIP_note[Chr_name].empty():
		# 	print Chr_name, CLIP_note[Chr_name].get()
		# print CLIP_note[Chr_name][6]
		cluster_pos_INS = sorted(cluster_pos_INS, key = lambda x:x[0])
		cluster_pos_DEL = sorted(cluster_pos_DEL, key = lambda x:x[0])
		if len(cluster_pos_INS) == 0:
			Cluster_INS = list()
		else:
			Cluster_INS = cluster(cluster_pos_INS, Chr_name)

		if len(cluster_pos_DEL) == 0:
			Cluster_DEL = list()
		else:
			Cluster_DEL = cluster_del(cluster_pos_DEL, Chr_name, Ref)

		logging.info("%d Alu signal locuses in the chromsome %s."%(len(Cluster_INS)+len(Cluster_DEL), Chr_name))

		# print Cluster_DEL
		# merge_siganl(Chr_name, Cluster)
		# break

		# merge step
		Final_result = combine_result(add_genotype(Cluster_INS, samfile), add_genotype(Cluster_DEL,samfile))

		# print Final_result

		# out_signal = open(p2, 'a+')
		# for i in Final_result:
		# 	for j in i:
		# 		out_signal.write(j)
		# out_signal.close()

		path = p2+Chr_name+'.fa'
		SeqIO.write(Final_result, path, "fasta")

	# out_signal = open(p2, 'w')
	# for i in total_signal:
	# 	out_signal.write(i)
	# out_signal.close()
	samfile.close()

def call_ngmlr(inFile, ref, presets, nproc, outFile):
	"""
	fq = input file
	automatically search for .sa
	"""
	outFile = outFile + "map.sam"
	logging.info("Running NGMLR...")
	cmd = ("ngmlr -r %s -q %s -o %s -t %d -x %s" % (ref, inFile, outFile, nproc, presets))
	r, o, e = exe(cmd)
	
	if r != 0:
		logging.error("NGMLR mapping failed!")
		logging.error("RETCODE %d" % (r))
		logging.error("STDOUT %s" % (str(o)))
		logging.error("STDERR %s" % (str(e)))
		logging.error("Exiting")
		exit(r)
	
	logging.info(str([r, o, e]))
	return outFile

def call_samtools(file, tempdir):
	logging.info("Running Samtools sort...")
	# samtools view -Sb ${WORK_PATH}_0001.sam | samtools sort -O bam -T /data/tjiang/ - > ${WORK_PATH}_0001.bam
	# samtools index ${WORK_PATH}_0001.bam
	bam_path = file[:-3] + "bam"
	cmd = ("samtools view -Sb %s | samtools sort -O bam -T %s - > %s" % (file, tempdir, bam_path))
	r, o, e = exe(cmd)
	if r != 0:
		logging.error("Samtools sort failed!")
		logging.error("RETCODE %d" % (r))
		logging.error("STDOUT %s" % (str(o)))
		logging.error("STDERR %s" % (str(e)))
		logging.error("Exiting")
		exit(r)
	logging.info(str([r, o, e]))

	logging.info("Running Samtools index...")
	cmd = ("samtools index %s" % (bam_path))
	r, o, e = exe(cmd)
	if r != 0:
		logging.error("Samtools index failed!")
		logging.error("RETCODE %d" % (r))
		logging.error("STDOUT %s" % (str(o)))
		logging.error("STDERR %s" % (str(e)))
		logging.error("Exiting")
		exit(r)
	logging.info(str([r, o, e]))


def single_pipe(out_path, chr, bam_path):
	# print out_path, chr_id, Ref, bam_path
	samfile = pysam.AlignmentFile(bam_path)
	# chr = samfile.get_reference_name(chr_id)
	logging.info("Resolving the chromsome %s."%(chr))
	if chr not in CLIP_note:
		CLIP_note[chr] = dict()
	cluster_pos_INS = list()
	cluster_pos_DEL = list()
	for read in samfile.fetch(chr):
		feed_back, feed_back_del = parse_read(read, chr)
		if len(feed_back) > 0:
			for i in feed_back:
				cluster_pos_INS.append(i)

		if len(feed_back_del) > 0:
			for i in feed_back_del:
				cluster_pos_DEL.append(i)
	cluster_pos_INS = sorted(cluster_pos_INS, key = lambda x:x[0])
	cluster_pos_DEL = sorted(cluster_pos_DEL, key = lambda x:x[0])
	if len(cluster_pos_INS) == 0:
		Cluster_INS = list()
	else:
		Cluster_INS = cluster(cluster_pos_INS, chr)

	if len(cluster_pos_DEL) == 0:
		Cluster_DEL = list()
	else:
		Ref = global_ref[0]
		Cluster_DEL = cluster_del(cluster_pos_DEL, chr, Ref)

	# print("[INFO]: %d Alu signal locuses in the chromsome %s."%(len(Cluster_INS)+len(Cluster_DEL), chr))
	logging.info("%d Alu signal locuses in the chromsome %s."%(len(Cluster_INS)+len(Cluster_DEL), chr))
	Final_result = combine_result(add_genotype(Cluster_INS, samfile), add_genotype(Cluster_DEL, samfile))
	samfile.close()
	# path = out_path+chr+'.fa'
	# SeqIO.write(Final_result, path, "fasta")
	return Final_result

def multi_run_wrapper(args):
   return single_pipe(*args)

def load_sam_multi_processes(args):
	'''
	Load_BAM_File
	library:	pysam.AlignmentFile

	load_Ref_Genome
	library:	Bio
	'''
	p1 = args.input
	p2 = args.temp_dir
	p3 = args.Reference
	threads = args.threads

	# Major Steps:
	# loading alignment file: bam format
	samfile = pysam.AlignmentFile(p1)
	# loading reference genome
	Ref = load_ref(p3)
	global_ref.append(Ref)

	# acquire the total numbers of the ref contigs
	contig_num = len(samfile.get_index_statistics())
	# print("[INFO]: The total number of chromsomes: %d"%(contig_num))
	logging.info("The total number of chromsomes: %d"%(contig_num))

	# Thread scheduling
	process_list = list()
	for i in samfile.get_index_statistics():
		process_list.append([i[0], i[3]])
		# #chr #read
	process_list = sorted(process_list, key = lambda x:-x[1])

	# start to establish multiprocesses
	analysis_pools = Pool(processes=threads)
	# Acquire_Chr_name
	Final_result = list()
	for i in process_list:
		# single_pipe(out_path, chr_id, Ref, samfile):
		# para = [(p2, _num_, Ref, samfile)]
		# print i
		para = [(p2, i[0], p1)]
		Final_result.append(analysis_pools.map_async(multi_run_wrapper, para))
	analysis_pools.close()
	analysis_pools.join()

	TE_list = list()
	for res in Final_result:
		# print res.get()[0]
		TE_list += res.get()[0]

	logging.info("Writing into disk")
	path = p2+'potential_TE.fa'
	# print TE_list
	SeqIO.write(TE_list, path, "fasta")

	samfile.close()

VERSION="1.0"

USAGE="""\

	some descriptions

"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="process.py extract", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	# parser.add_argument("AlignmentFile", type=str, help="the bam format file generated by ngmlr, within a '.bai' index file")
	parser.add_argument("input", metavar="[SAM,BAM,FASTA,FASTQ,FOFN]", type=str, help="Input reads to be mapped")
	parser.add_argument("Reference", metavar="REFERENCE", type=str, help="the reference genome(fasta format)")
	parser.add_argument('temp_dir', type=str, help = "temporary directory to use for distributed jobs")
	# parser.add_argument("--temp", type=str, default=tempfile.gettempdir(), help="Where to save temporary files")

	parser.add_argument('-s', '--min_support', help = "Mininum number of reads that support a TE.[%(default)s]", default = 5, type = int)
	parser.add_argument('-l', '--min_length', help = "Mininum length of TE to be reported.[%(default)s]", default = 50, type = int)
	parser.add_argument('-d', '--min_distance', help = "Mininum distance of two TE clusters.[%(default)s]", default = 20, type = int)
	# parser.add_argument('-hom', '--homozygous', help = "The mininum score of a genotyping reported as a homozygous.[%(default)s]", default = 0.8, type = float)
	# parser.add_argument('-het','--heterozygous', help = "The mininum score of a genotyping reported as a heterozygous.[%(default)s]", default = 0.3, type = float)
	# parser.add_argument('-q', '--min_mapq', help = "Mininum mapping quality.[20]", default = 20, type = int)
	parser.add_argument('-t', '--threads', help = "Number of threads to use.[%(default)s]", default = 1, type = int)
	parser.add_argument('-x', '--presets', help = "The sequencing type <pacbio,ont> of the reads.[%(default)s]", default = "pacbio", type = str)
	# parser.add_argument("--temp", type=str, default=tempfile.gettempdir(), help="Where to save temporary files")
	# parser.add_argument("--chunks", type=int, default=0, help="Create N scripts containing commands to each input of the fofn (%(default)s)")
	# parser.add_argument("--debug", action="store_true")
	args = parser.parse_args(argv)

	# setupLogging(args.debug)
	# checkBlasrParams(args.params)
	
	# if args.output is None:
	# 	ext =  args.input[args.input.rindex('.'):]
	# 	main = args.input[:args.input.rindex('.')]
	# 	if ext in [".sam", ".bam"]:
	# 		args.output = main + ".tails" + ext
	# 	else:
	# 		args.output = main + ".tails.sam"
	return args


def decipherInput(input):
    """
    returns True if initial map needs to happen
    and list of inputFileNames 
    in input.fofn and chunks, you'll have lists of lists
    """
    extension = input.split('.')[-1].lower()

    if extension == "bam":
    	# if chunks != 0:
    	# 	logging.error("chunks not applicable to %s files" % extension)
    	return 0

    if extension == "sam":
    	# if chunks != 0:
    	# 	logging.error("chunks not applicable to %s files" % extension)
    	return 1

    if extension in ["fastq", "fasta", "fa", "fq"]:
        # if chunks != 0:
        #     logging.error("chunks not applicable to %s files" % extension)
        #     exit(1)
        return 2
    
    # if extension == "fofn":
    #     inputs = [x.strip() for x in open(input).readlines()]
    #     if chunks != 0:
    #         return True, partition(inputs, chunks)
    #     else:
    #         return True, [input]

def run(argv):
    args = parseArgs(argv)
    setupLogging(False)
    # print args
    starttime = time.time()
    flag = decipherInput(args.input)
    
    if flag == 0:
    	# detection
    	load_sam_multi_processes(args)

    elif flag == 1:
    	call_samtools(args.input, args.temp_dir)
    	# samtools transfer
    else:
    	# inFile, ref, seq_type, nproc=1, outFile="map.sam", presets="pacbio"
    	file = call_ngmlr(args.input, args.Reference, args.presets, args.threads, args.temp_dir)
    	call_samtools(file, args.temp_dir)
    	# ngmlr mapping
    logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

if __name__ == '__main__':
    run(sys.argv[:1])