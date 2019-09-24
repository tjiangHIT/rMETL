#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  call_TE.py
 * @Package: argparse, sys, logging, pysam, Bio
 * @Description: Establish the ME callset
 * @author: Jiang Tao (tjiang@hit.edu.cn)
 * @date: Apr 24 2018
 * @version V1.0.2  
'''

import argparse
import logging
import sys
import time
import cigar

from collections import Counter
from rMETL.rMETL_version import __version__, __author__, __contact__
from rMETL.rMETL_genotype import simple_call_genotype
from rMETL.rMETL_cmdRunner import setupLogging
from rMETL.rMETL_utils import load_ref

USAGE="""\
           _  ___  _   _____   _______   _
     _ _  | ^_   _^ | |  ___| |__   __| | |
    | ^_| | | | | | | | |__      | |    | |
    | |   | | | | | | |  __|     | |    | |
    | |   | | | | | | | |___     | |    | |___
    |_|   |_| |_| |_| |_____|    |_|    |_____|

    rMETL - realignment-based Mobile Element insertion detection Tool for Long read

	Generate final MEI/MED callset in bed or vcf file.
	
	The output file called 'calling.bed' or 'calling.vcf'
	stores in output directory.
	
	rMETL V%s
	Author: %s
	Contact: %s
"""%(__version__, __author__, __contact__)

def acquire_count_max(_list_):
	c = Counter(_list_)
	return c.most_common(1)[0][0]

flag_dic = {0:1, \
			16:2, \
			256:0, \
			272:0, \
			2048:0, \
			2064:0, \
			4:0}

STRAND = {'1':'+', \
		  '2':'-', \
		  '*':'+-'}

cluster_dic = {}

strand_dic = {1:'+', \
			  2:'-'}

class R_INFO(object):
	"""store the infomation of the signal sequence"""
	def __init__(self, Type, Chr, Pos, Len, GT):
		self.Type = Type
		self.Chr = Chr
		self.Pos = Pos
		self.Len = Len
		self.GT = GT

def parse_name(seq):
	chr = seq.split('*')[0]
	breakpoint = seq.split('*')[1]
	insert_size = seq.split('*')[2]
	GT = seq.split('*')[3]
	return chr, breakpoint, insert_size, GT

def parse_name_tp(line):
	'''
	resolution signatures for bed format
	'''
	seq = line.split('*')
	Type = seq[0]
	chr = seq[1]
	pos = seq[2]
	len = seq[3]
	if Type == 'DEL':
		rc = seq[4]
		cov = seq[5]
	else:
		rc = seq[5]
		cov = seq[6]
	GT = rc+':'+cov
	local_info = R_INFO(Type, chr, pos, len, GT)
	return local_info

def clip_analysis(deal_cigar, clipping_threshold):
	'''
	resolution cogar
	'''
	seq = list(cigar.Cigar(deal_cigar).items())
	if seq[0][1] == 'S':
		first_pos = seq[0][0]
	else:
		first_pos = 0
	if seq[-1][1] == 'S':
		last_pos = seq[-1][0]
	else:
		last_pos = 0
	total_len = first_pos + last_pos
	signal_len = 0
	for i in seq:
		signal_len += i[0]
	if signal_len == 0:
		return 0
	if total_len*1.0 / signal_len >= clipping_threshold:
		return 0
	else:
		return 1

def print_vcf_head(ref, sample):
	'''
	generation of VCF head
	'''
	import time
	Date = time.strftime("%Y%m%d")
	head = list()
	head.append("##fileformat=VCFv4.2\n")
	head.append("##fileDate=%s\n"%(Date))
	head.append("##source=rMETL\n")
	for i in ref:
		head.append("##contig=<ID=%s,length=%d>\n"%(i, len(ref[i])))
	head.append("##ALT=<ID=<DEL>,Description=\"Deletion relative to the reference\">\n")
	head.append("##ALT=<ID=<INS>,Description=\"Insertion of sequence relative to the reference\">\n")
	head.append("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
	head.append("##INFO=<ID=SVLEN,Number=.,Type=String,Description=\"Difference in length between REF and ALT alleles\">\n")
	head.append("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n")
	head.append("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n")
	head.append("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count'\">\n")
	head.append("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele frequency'\">\n")
	head.append("##INFO=<ID=AN,Number=.,Type=String,Description=\"Allele name'\">\n")
	head.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
	head.append("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"#High-quality variant reads\">\n")
	head.append("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"#Reference reads\">\n")
	head.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n"%(sample))
	return head

def parse_seq_head(line):
	'''
	resolution signatures for vcf format
	'''
	seq = line.split('*')
	Type = seq[0]
	chr = seq[1]
	pos = seq[2]
	len = seq[3]
	if Type == 'DEL':
		rc = seq[4]
		cov = seq[5]
	else:
		rc = seq[5]
		cov = seq[6]
	GT = rc+':'+cov
	local_info = R_INFO(Type, chr, pos, len, GT)
	return local_info

# ************************BED_FUNCTION*******************************
def call_bed(args):
	path = args.input
	out_path = args.output + "calling.bed"
	AlignmentFile = open(path, 'r')
	logging.info("Loading ME realignmets...")
	for line in AlignmentFile:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '@':
			continue
		local_info = parse_name_tp(seq[0])
		Flag = int(seq[1])
		sub_type = seq[2]
		MAPQ = int(seq[4])
		cigar = seq[5]
		cigar_flag = clip_analysis(cigar, args.clipping_threshold)
		if flag_dic[Flag] != 0 and MAPQ >= args.min_mapq and cigar_flag == 1:
			key = "%s*%s*%s*%s"%(local_info.Chr, local_info.Pos, local_info.Len, \
				local_info.GT)
			if key not in cluster_dic:
				cluster_dic[key] = list()
			cluster_dic[key].append("<%s:ME:%s>"%(local_info.Type, sub_type))
	AlignmentFile.close()

	if args.MEI == 'False':
		AlignmentFile = open(path, 'r')
		for line in AlignmentFile:
			seq = line.strip('\n').split('\t')
			if seq[0][0] == '@':
				continue
			local_info = parse_name_tp(seq[0])
			Flag = int(seq[1])
			sub_type = seq[2]
			if sub_type == '*':
				key = "%s*%s*%s*%s"%(local_info.Chr, local_info.Pos, local_info.Len, \
					local_info.GT)
				if key not in cluster_dic:
					cluster_dic[key] = list()
					cluster_dic[key].append("<%s>"%(local_info.Type))
		AlignmentFile.close()

	sort_list = list()
	for i in cluster_dic:
		chr, breakpoint, insert_size, GT = parse_name(i)
		final_type = acquire_count_max(cluster_dic[i])
		sort_list.append([chr, breakpoint, insert_size, final_type])
	sort_list = sorted(sort_list, key = lambda x:(x[0], int(x[1])))
	file = open(out_path, 'w')
	logging.info("Writing results into disk...")
	file.write("# Chromsome\tBreakpoint\tSV length\tMEI Type")
	for i in sort_list:
		file.write("\t".join(i)+"\n")
	file.close()
# ************************BED_FUNCTION*******************************
# 
# 
# 
# ************************VCF_FUNCTION*******************************
def call_vcf(args):
	path = args.input
	out_path = args.output + "calling.vcf"
	ref = load_ref(args.Reference)
	AlignmentFile = open(path, 'r')
	logging.info("Loading ME realignmets...")
	for line in AlignmentFile:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '@':
			continue
		local_info = parse_seq_head(seq[0])
		Flag = int(seq[1])
		sub_type = seq[2]
		MAPQ = int(seq[4])
		cigar = seq[5]
		cigar_flag = clip_analysis(cigar, args.clipping_threshold)
		if flag_dic[Flag] != 0 and MAPQ >= args.min_mapq and cigar_flag == 1:
			key = "%s*%s*%s*%s"%(local_info.Chr, local_info.Pos, local_info.Len, \
				local_info.GT)
			if key not in cluster_dic:
				cluster_dic[key] = list()
			cluster_dic[key].append("<%s:ME:%s>\t%d"%(local_info.Type, sub_type, \
				flag_dic[Flag]))
	AlignmentFile.close()

	if args.MEI == 'False':
		AlignmentFile = open(path, 'r')
		for line in AlignmentFile:
			seq = line.strip('\n').split('\t')
			if seq[0][0] == '@':
				continue
			local_info = parse_name_tp(seq[0])
			Flag = int(seq[1])
			sub_type = seq[2]
			if sub_type == '*':
				key = "%s*%s*%s*%s"%(local_info.Chr, local_info.Pos, local_info.Len, \
					local_info.GT)
				if key not in cluster_dic:
					cluster_dic[key] = list()
					cluster_dic[key].append("<%s>\t%s"%(local_info.Type, '*'))
		AlignmentFile.close()

	sort_list = list()
	for i in cluster_dic:
		chr, breakpoint, insert_size, GT = parse_name(i)
		final_type = acquire_count_max(cluster_dic[i]).split('\t')[0]
		strand = STRAND[acquire_count_max(cluster_dic[i]).split('\t')[1]]
		sort_list.append([chr, breakpoint, insert_size, final_type, GT, strand])
	sort_list = sorted(sort_list, key = lambda x:(x[0], int(x[1])))
	head_info = print_vcf_head(ref, args.sample)

	file = open(out_path, 'w')
	logging.info("Writing results into disk...")

	for line in head_info:
		file.write(line)

	ID = 0
	for i in sort_list:
		concordant = int(i[4].split(':')[0])
		discordant = int(i[4].split(':')[1]) - int(i[4].split(':')[0])
		if discordant < 0:
			discordant = 0
		GT, GL, reliability = simple_call_genotype(concordant, concordant+discordant, \
			args.heterozygous, args.homozygous)

		if reliability == 1:
			INFO = "PRECISE;SVTYPE=%s;SVLEN=%d;END=%d;SAMPLE=%s;STRAND=%s"%(i[3][1:4], \
				int(i[2]), int(i[1])+int(i[2])-1, args.sample, i[5])
		else:
			INFO = "IMPRECISE;SVTYPE=%s;SVLEN=%d;END=%d;SAMPLE=%s;STRAND=%s"%(i[3][1:4], \
				int(i[2]), int(i[1])+int(i[2])-1, args.sample, i[5])
		try:
			REF = ref[i[0]][int(i[1])-1]
		except:
			REF = "N"
		file.write("%s\t%s\t%d\t%s\t%s\t.\tPASS\t%s\tGT:DV:DR\t%s:%s\n"%(i[0], i[1], \
			ID, REF, i[3], INFO, GT, GL))
		ID += 1
	file.close()
# *************************VCF_FUNCTION*******************************
# 
# 
# 
# ************************MAIN_FUNCTION*******************************
def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="rMETL.py calling", description=USAGE, \
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("input", metavar="SAM", type=str, help="Input cluster.sam on STAGE realignment.")
	parser.add_argument("Reference", metavar="REFERENCE", type=str, \
		help="The reference genome in fasta format.")
	parser.add_argument("format", metavar="[BED,VCF]", type=str, \
		help="The format of the output file. [%(default)s]", default = "bed")
	parser.add_argument('output', type=str, help = "Directory to output final callset.")
	parser.add_argument('-hom', '--homozygous', \
		help = "The mininum score of a genotyping reported as a homozygous.[%(default)s]", \
		default = 0.8, type = float)
	parser.add_argument('-het','--heterozygous', \
		help = "The mininum score of a genotyping reported as a heterozygous.[%(default)s]", \
		default = 0.3, type = float)
	parser.add_argument('-q', '--min_mapq', help = "Mininum mapping quality.[%(default)s]", \
		default = 20, type = int)
	parser.add_argument('-c', '--clipping_threshold', \
		help = "Mininum threshold of realignment clipping.[%(default)s]", \
		default = 0.5, type = float)
	parser.add_argument('--sample', help = "Sample description", \
		default = "None", type = str)
	parser.add_argument('--MEI', help = "Enables rMETL to display MEI/MED only.[%(default)s]", \
		default = "True", type = str)
	args = parser.parse_args(argv)
	return args

def run(argv):
    args = parseArgs(argv)
    setupLogging(False)
    starttime = time.time()
    if args.format == "bed":
    	call_bed(args)
    elif args.format == "vcf":
    	call_vcf(args)
    else:
    	logging.error("Invalid format.")
    	exit(1)
    logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

if __name__ == '__main__':
	run(sys.argv[:1])
