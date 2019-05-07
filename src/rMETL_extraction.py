#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  extract_simple.py
 * @Package: argparse, pysam, sys, Bio, os, logging
 * @Description: Parse the ME signatures from alignments
 * @author: Jiang Tao (tjiang@hit.edu.cn)
 * @date: Apr 24 2018
 * @version V1.0.2     
'''

import pysam
import cigar
import os
import argparse
import logging
import sys
import time
import gc

from multiprocessing import Pool
from rMETL_version import __version__, __author__, __contact__
from rMETL_concensus import construct_concensus_info
from rMETL_genotype import add_genotype
from rMETL_utils import load_ref, check_bai, call_ngmlr, call_samtools
from rMETL_cmdRunner import setupLogging, exe

USAGE="""\
           _  ___  _   _____   _______   _
     _ _  | ^_   _^ | |  ___| |__   __| | |
    | ^_| | | | | | | | |__      | |    | |
    | |   | | | | | | |  __|     | |    | |
    | |   | | | | | | | |___     | |    | |___
    |_|   |_| |_| |_| |_____|    |_|    |_____|

    rMETL - realignment-based Mobile Element insertion detection Tool for Long read


	Support reads aligned with Ngmlr and sorted with Samtools

	If input is a fastq or fasta format file, rMETL generates
	alignments with Ngmlr at first;

	If input is a sam format file, rMETL converts and sorts it
	to be a bam format file;

	If your input is a bam format file with index, rMETL extracts
	the ME signatures and collects the sub-sequence of them.

	The output is a fasta format file called 'potential.fa' 
	contains potentials non-reference ME clusters.

	rMETL V%s
	Author: %s
	Contact: %s
"""%(__version__, __author__, __contact__)

INS_flag = {1:'I'}
DEL_flag = {2:'D'}
clip_flag = {4:'S', 5:'H'}
global_ref = list()

# **********************check-input-format****************************
def decipherInput(input):
	"""
	resolution input format
	"""
	extension = input.split('.')[-1].lower()
	choice = {"bam": 0, \
			  "sam": 1, \
			  "fasta": 2, \
			  "fastq": 2, \
			  "fa": 2, \
			  "fq": 2}
	return choice[extension]
# **********************check-input-format****************************
# 
# 
# 
# ************************mini-operations*****************************
def revcom_complement(s): 
	'''
	generation reverse complementary sequence.
	all of the lowercase will be changed to capital letter.
	'''
	basecomplement = {'A': 'T', \
					  'C': 'G', \
					  'G': 'C', \
					  'T': 'A', \
					  'a': 'T', \
					  'c': 'G', \
					  'g': 'C', \
					  't': 'A'} 
	letters = list(s) 
	letters = [basecomplement[base] for base in letters] 
	return ''.join(letters)[::-1]

def detect_flag(Flag):
	'''
	identification of flag in BAM
		0 means unmapped read
		1 & 2 means primary mapping read with normal strand or reverse strand 
		3 & 4 means supplementary mapping read with normal strand or reverse strand
	'''
	Normal_foward = 1 >> 1
	Abnormal = 1 << 2
	Reverse_complement = 1 << 4
	Supplementary_map = 1 << 11

	signal = {Abnormal: 0, \
			  Normal_foward: 1, \
			  Reverse_complement: 2, \
			  Supplementary_map: 3, \
			  Reverse_complement | Supplementary_map: 4}

	back_sig = signal[Flag] if Flag in signal else 0
	return back_sig

def acquire_clip_pos(deal_cigar):
	'''
	resolution of cigar in supplementary mapping
	'''
	seq = list(cigar.Cigar(deal_cigar).items())
	first_pos = seq[0][0] if seq[0][1] == 'S' else 0
	last_pos = seq[-1][0] if seq[-1][1] == 'S' else 0
	bias = 0
	for i in seq:
		if i[1] in ['M', 'D']:
			bias += i[0]
	return [first_pos, last_pos, bias]
# ************************mini-operations*****************************
# 
# 
# 
# ************************soft-clippings******************************
def store_clip_pos(locus, chr, seq, flag, CLIP_note):
	'''
	A data structure store info of soft-clippings.
	It has two hashtables:
		1. key1 is an integer calculated by coordinate / 10000;
		2. key2 is an integer calculated by (coordinate % 10000) / 50.
	'''
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

def acquire_clip_locus(down, up, chr, CLIP_note):
	'''
	search soft-clippings within limited region
	'''
	list_clip = list()
	if int(up/10000) == int(down/10000):
		key_1 = int(down/10000)
		if key_1 not in CLIP_note[chr]:
			return list_clip
		for i in xrange(int((up%10000)/50)-int((down%10000)/50)+1):
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
			key_2 = i
			if key_2 not in CLIP_note[chr][key_1]:
				continue
			for ele in CLIP_note[chr][key_1][key_2]:
				if ele[0] >= down and ele[0] <= up:
					list_clip.append(ele)
	return list_clip
# ************************soft-clippings******************************
# 
# 
# 
# ***********************resolution-reads*****************************
def organize_split_signal(chr, primary_info, Supplementary_info, \
	total_L, low_bandary):
	'''
	resolution split alignments
	'''
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
			continue
		local_set = acquire_clip_pos(local_cigar)
		if primary_info[0] < local_start:
			if primary_info[3]+local_set[0]-total_L > low_bandary:
				overlap.append([total_L - primary_info[3], \
					local_set[0], primary_info[1]])
		else:
			if local_set[1]+primary_info[2]-total_L > low_bandary:
				overlap.append([total_L - local_set[1], \
					primary_info[2], local_start+local_set[2]-1])
	return overlap

def parse_read(read, Chr_name, low_bandary, CLIP_note):
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
	pos_start = read.reference_start
	shift = 0
	for element in read.cigar:
		if element[0] == 0:
			shift += element[1]
		if element[0] in DEL_flag and element[1] <= low_bandary:
			shift += element[1]
		if element[0] in DEL_flag and element[1] > low_bandary:
			DEL_ME_pos.append([pos_start+shift, element[1]])
			shift += element[1]

	# Add INS:ME type call signal
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
			MEI_contig = read.query_sequence[_shift_read_ - \
			element[1]:_shift_read_]
			INS_ME_pos.append([pos_start + shift, element[1], \
				MEI_contig])
		if element[0] in clip_flag:
			if shift == 0:
				primary_clip_0 = element[1]
			else:
				primary_clip_1 = element[1]
			if element[1] > low_bandary:
				if shift == 0:
					clip_pos = pos_start - 1
					clip_contig = read.query_sequence[:element[1]]
					store_clip_pos(clip_pos, Chr_name, clip_contig, \
						0, CLIP_note)
				else:
					clip_pos = pos_start + shift - 1
					clip_contig = read.query_sequence[read.query_length \
					- element[1]:]
					store_clip_pos(clip_pos, Chr_name, clip_contig, 1, \
						CLIP_note)

	if process_signal in [1, 2]:
		Tags = read.get_tags()
		chr = Chr_name
		primary_info = [pos_start, pos_end, primary_clip_0, primary_clip_1, \
		process_signal]
		for i in Tags:
			if i[0] == 'SA':
				Supplementary_info = i[1].split(';')[:-1]
				overlap = organize_split_signal(chr, primary_info, \
					Supplementary_info, read.query_length, low_bandary)
				for k in overlap:
					MEI_contig = read.query_sequence[k[0]:k[1]]
					INS_ME_pos.append([k[2], k[1] - k[0], MEI_contig])
	return INS_ME_pos, DEL_ME_pos
# ***********************resolution-reads*****************************
# 
# 
# 
# ***********************Cluster-Function*****************************
def merge_pos(pos_list, chr, evidence_read, SV_size, CLIP_note):
	'''
	INS: inner-cluster function
	'''
	start = list()
	end = list()
	for ele in pos_list:
		start.append(ele[0])
		end.append(ele[0] + ele[1])
	search_down = min(start) - 10
	search_up = max(start) + 10
	temp_clip = acquire_clip_locus(search_down, search_up, chr, CLIP_note)
	result = construct_concensus_info(pos_list, temp_clip, evidence_read, \
		SV_size)
	if result != 0:
		for i in xrange(len(result)):
			result[i] = ["INS", chr] + result[i] + [len(result)]
		return result
	else:
		return 0

def cluster(pos_list, chr, evidence_read, SV_size, low_bandary, CLIP_note):
	'''
	INS: outer-cluster function
	'''
	_cluster_ = list()
	temp = list()
	temp.append(pos_list[0])
	for pos in pos_list[1:]:
		if temp[-1][0] + low_bandary < pos[0]:
			result = merge_pos(temp, chr, evidence_read, SV_size, CLIP_note)
			if result != 0:
				_cluster_.append(result)
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	result = merge_pos(temp, chr, evidence_read, SV_size, CLIP_note)
	if result != 0:
		_cluster_.append(result)
	return _cluster_

def merge_pos_del(pos_list, chr, Ref, evidence_read, SV_size):
	'''
	DEL: inner-cluster function
	'''
	start = list()
	end = list()
	for ele in pos_list:
		start.append(ele[0])
		end.append(ele[0] + ele[1])
	breakpoint = sum(start)/len(start)
	size = sum(end)/len(end) - breakpoint
	result = list()
	if len(pos_list) < evidence_read:
		return result
	else:
		if chr in Ref and size >= SV_size:
			result.append(['DEL', chr, breakpoint, size, len(pos_list), \
				str(Ref[chr].seq[breakpoint:breakpoint+size])])
	return result

def cluster_del(pos_list, chr, Ref, evidence_read, SV_size, low_bandary):
	'''
	DEL: outer-cluster function
	'''
	_cluster_ = list()
	temp = list()
	temp.append(pos_list[0])
	for pos in pos_list[1:]:
		if temp[-1][0] + low_bandary < pos[0]:
			result = merge_pos_del(temp, chr, Ref, evidence_read, SV_size)
			if len(result) != 0:
				_cluster_.append(result)
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	result = merge_pos_del(temp, chr, Ref, evidence_read, SV_size)
	if len(result) != 0:
		_cluster_.append(result)
	return _cluster_
# ***********************Cluster-Function*****************************
# 
# 
# 
# ***********************Output-Function******************************
def combine_result(INS, DEL, path, chr):
	'''
	Merge results into one list and output it.
	'''
	output = "%ssignatures/%s_sig.fa"%(path, chr)
	file = open(output, 'w')
	for i in INS:
		for j in i:
			if len(j) != 8:
				continue
			key = "%s*%s*%d*%d*%s*%d*%d"%(j[0], j[1], j[2], j[3], j[4], j[6], \
				j[7])
			file.write(">"+key+'\n')
			file.write(j[5]+'\n')
	del INS
	gc.collect()
	for i in DEL:
		for j in i:
			if len(j) != 7:
				continue
			key = "%s*%s*%d*%d*%d*%d"%(j[0], j[1], j[2], j[3], j[4], j[6])
			file.write(">%s\n"%(key))
			file.write('%s\n'%(j[5]))
	del DEL
	gc.collect()
	file.close()
# ***********************Output-Function******************************
# 
# 
# 
# ********************Signatures-extraction***************************
def single_pipe(out_path, chr, bam_path, low_bandary, evidence_read, SV_size):
	'''
	resolution signatures
	'''
	samfile = pysam.AlignmentFile(bam_path)
	CLIP_note = dict()
	logging.info("Resolving chromsome %s."%(chr))
	if chr not in CLIP_note:
		CLIP_note[chr] = dict()
	cluster_pos_INS = list()
	cluster_pos_DEL = list()
	for read in samfile.fetch(chr):
		feed_back, feed_back_del = parse_read(read, chr, low_bandary, CLIP_note)
		cluster_pos_INS += feed_back
		cluster_pos_DEL += feed_back_del
	cluster_pos_INS = sorted(cluster_pos_INS, key = lambda x:x[0])
	cluster_pos_DEL = sorted(cluster_pos_DEL, key = lambda x:x[0])
	if len(cluster_pos_INS) == 0:
		Cluster_INS = list()
	else:
		Cluster_INS = cluster(cluster_pos_INS, chr, evidence_read, SV_size, \
			low_bandary, CLIP_note)
		del cluster_pos_INS
		del CLIP_note[chr]
		gc.collect()
	if len(cluster_pos_DEL) == 0:
		Cluster_DEL = list()
	else:
		Ref = global_ref[0]
		Cluster_DEL = cluster_del(cluster_pos_DEL, chr, Ref, evidence_read, \
			SV_size, low_bandary)
		del cluster_pos_DEL
		gc.collect()
	logging.info("%d MEI/MED signal loci in the chromsome %s."%(len(Cluster_INS)+\
		len(Cluster_DEL), chr))
	combine_result(add_genotype(Cluster_INS, samfile, low_bandary), \
		add_genotype(Cluster_DEL, samfile, low_bandary), out_path, chr)
	samfile.close()

def multi_run_wrapper(args):
   return single_pipe(*args)

def load_sam_multi_processes(args):
	'''
	task scheduling
	'''
	temporary_dir = args.temp_dir if args.temp_dir.endswith('/') else \
	"%s/"%(args.temp_dir)
	os.mkdir("%ssignatures"%temporary_dir)
	# Major Steps:
	# loading alignment file: bam format
	samfile = pysam.AlignmentFile(args.input)
	# loading reference genome
	Ref = load_ref(args.Reference)
	global_ref.append(Ref)
	# acquire the total numbers of the ref contigs
	contig_num = len(samfile.get_index_statistics())
	logging.info("The total number of chromsomes: %d"%(contig_num))
	# Thread scheduling
	process_list = list()
	for i in samfile.get_index_statistics():
		process_list.append([i[0], i[3]])
		# #chr #read
	process_list = sorted(process_list, key = lambda x:-x[1])
	# start to establish multiprocesses
	analysis_pools = Pool(processes = args.threads)
	# Acquire_Chr_name
	for i in process_list:
		para = [(temporary_dir, i[0], args.input, args.min_distance, \
			args.min_support, args.min_length)]
		analysis_pools.map_async(multi_run_wrapper, para)
	analysis_pools.close()
	analysis_pools.join()
	samfile.close()

	output_p = args.output if args.output.endswith('/') else "%s/"%(args.output)
	if not os.path.exists(output_p):
		os.mkdir(output_p)
	merge_cmd = ("cat %ssignatures/* > %spotential_ME.fa"%(temporary_dir, output_p))
	r, o, e = exe(merge_cmd)
	if r != 0:
		logging.error("Merging ME signatures failed!")
		logging.error("RETCODE %d" % (r))
		logging.error("STDOUT %s" % (str(o)))
		logging.error("STDERR %s" % (str(e)))
		logging.error("Exiting")
		exit(r)
	logging.info("Cleaning temporary files.")
	cmd_remove_tempfile = ("rm -r %ssignatures"%(temporary_dir))
	r, o, e = exe(cmd_remove_tempfile)
	if r != 0:
		logging.error("Cleaning temporary files failed!")
		logging.error("RETCODE %d" % (r))
		logging.error("STDOUT %s" % (str(o)))
		logging.error("STDERR %s" % (str(e)))
		logging.error("Exiting")
		exit(r)
# ********************Signatures-extraction***************************
# 
# 
# 
# ************************MAIN_FUNCTION*******************************
def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="rMETL.py detection", \
		description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("input", metavar="[SAM,BAM,FASTA,FASTQ]", type=str, \
		help="Input reads with/without alignment.")
	parser.add_argument("Reference", metavar="REFERENCE", type=str, \
		help="The reference genome in fasta format.")
	parser.add_argument('temp_dir', type=str, \
		help = "Temporary directory to use for distributed jobs.")
	parser.add_argument('output_dir', type=str, \
		help = "Directory to output potential ME loci.")
	parser.add_argument('-s', '--min_support',\
	 help = "Mininum number of reads that support a ME.[%(default)s]", \
	 default = 5, type = int)
	parser.add_argument('-l', '--min_length', \
		help = "Mininum length of ME to be reported.[%(default)s]", \
		default = 50, type = int)
	parser.add_argument('-d', '--min_distance', \
		help = "Mininum distance of two ME signatures to be intergrated.[%(default)s]", \
		default = 20, type = int)
	parser.add_argument('-t', '--threads', \
		help = "Number of threads to use.[%(default)s]", default = 8, \
		type = int)
	parser.add_argument('-x', '--presets', \
		help = "The sequencing platform <pacbio,ont> of the reads.[%(default)s]", \
		default = "pacbio", type = str)
	args = parser.parse_args(argv)
	return args

def run(argv):
	args = parseArgs(argv)
	setupLogging(False)
	starttime = time.time()
	flag = decipherInput(args.input)
	
	if flag == 0:
		# detection
		result = check_bai(args.input, args.temp_dir)
		if len(result) == 0:
			load_sam_multi_processes(args)
		else:
			args.input = result
			load_sam_multi_processes(args)
	elif flag == 1:
		bam_path = call_samtools(args.input, args.temp_dir)
		args.input = bam_path
		load_sam_multi_processes(args)
		# samtools transfer
	else:
		# inFile, ref, seq_type, nproc=1, outFile="map.sam", presets="pacbio"
		file = call_ngmlr(args.input, args.Reference, args.presets, \
			rgs.threads, args.temp_dir)
		bam_path = call_samtools(file, args.temp_dir)
		args.input = bam_path
		load_sam_multi_processes(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

if __name__ == '__main__':
	run(sys.argv[:1])
