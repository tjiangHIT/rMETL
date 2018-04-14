#!/usr/bin/env python

import pysam
from collections import Counter
from genotype import *
import os
import argparse
import logging
import tempfile
import sys, time
from CommandRunner import *
from Bio import SeqIO

def acquire_count_max(_list_):
	c = Counter(_list_)
	return c.most_common(1)[0][0]

flag_dic = {0:1, 16:2, 2048:0, 2064:0, 4:0}

STRAND = {'1':'+', '2':'-'}

cluster_dic = {}

strand_dic = {1:'+', 2:'-'}

class R_INFO(object):
	"""store the infomation of the signal sequence"""
	def __init__(self, Type, Chr, Pos, Len, GT):
		# super(ClassName, self).__init__()
		self.Type = Type
		self.Chr = Chr
		self.Pos = Pos
		self.Len = Len
		self.GT = GT

def parse_name(seq):
	chr = seq.split('_')[0]
	breakpoint = seq.split('_')[1]
	insert_size = seq.split('_')[2]
	GT = seq.split('_')[3]
	return chr, breakpoint, insert_size, GT

def load_sam(p1):
	# samfile = pysam.AlignmentFile(p1)
	AlignmentFile = open(p1, 'r')
	for line in AlignmentFile:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '@':
			continue
		chr, breakpoint, insert_size = parse_name(seq[0])
		Flag = int(seq[1])
		sub_type = seq[2]
		if flag_dic[Flag] != 0 and int(insert_size) >= 50:
			# to do something
			key = "%s_%s_%s"%(chr, breakpoint, insert_size)
			if key not in cluster_dic:
				cluster_dic[key] = list()
			# cluster_dic[key].append(sub_type+"&"+str(flag_dic[Flag]))
			cluster_dic[key].append(sub_type)
	AlignmentFile.close()

	sort_list = list()
	for i in cluster_dic:
		chr, breakpoint, insert_size = parse_name(i)
		final = acquire_count_max(cluster_dic[i])
		final_type = final.split('&')[0]
		# final_strand = final.split('&')[1]
		# final_strand = acquire_count_max(cluster_dic[i][1])
		# sort_list.append([chr, breakpoint, insert_size, final_type, strand_dic[int(final_strand)]])
		sort_list.append([chr, breakpoint, insert_size, final_type])
		# print("%s\t%s\t%s\t%s"%(chr, breakpoint, insert_size, final_type))
	sort_list = sorted(sort_list, key = lambda x:(x[0], int(x[1])))
	for i in sort_list:
		# print("%s\t%s\t%s\t%s"%(i[0], breakpoint, insert_size, final_type))
		print "\t".join(i)

def parse_name_tp(line):
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

def call_bed(args):
	# samfile = pysam.AlignmentFile(p1)
	path = args.input
	out_path = args.temp_dir + "calling.bed"
	AlignmentFile = open(path, 'r')
	logging.info("Loading TE alignmets...")
	for line in AlignmentFile:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '@':
			continue

		local_info = parse_name_tp(seq[0])
		Flag = int(seq[1])
		sub_type = seq[2]
		MAPQ = int(seq[4])
		if flag_dic[Flag] != 0 and MAPQ >= args.min_mapq:
			# to do something
			# key = "%s_%s_%s"%(chr, breakpoint, insert_size)
			key = "%s_%s_%s_%s"%(local_info.Chr, local_info.Pos, local_info.Len, local_info.GT)
			if key not in cluster_dic:
				cluster_dic[key] = list()
			# cluster_dic[key].append("%s:ME:%s"%(Stype,sub_type))
			cluster_dic[key].append("<%s:ME:%s>"%(local_info.Type, sub_type))

			# if key not in cluster_dic:
			# 	cluster_dic[key] = [MAPQ, sub_type]
			# else:
			# 	if MAPQ > cluster_dic[key][0]:
			# 		cluster_dic[key] = [MAPQ, sub_type]

	AlignmentFile.close()

	sort_list = list()
	for i in cluster_dic:
		chr, breakpoint, insert_size, GT = parse_name(i)
		# print cluster_dic[i]
		final_type = acquire_count_max(cluster_dic[i])
		# final_type = cluster_dic[i][1]
		# final_MAPQ = cluster_dic[i][0]
		# final_strand = final.split('&')[1]
		# final_strand = acquire_count_max(cluster_dic[i][1])
		# sort_list.append([chr, breakpoint, insert_size, final_type, strand_dic[int(final_strand)]])

		# if final_MAPQ >= 20:
		# sort_list.append([chr, breakpoint, insert_size, final_type])
		# concordant = int(GT.split(':')[0])
		# coverage = int(GT.split(':')[1])
		# flag = simple_filter_genotype(concordant, coverage, 0.2)
		# if flag == 0:
		# 	continue
		sort_list.append([chr, breakpoint, insert_size, final_type])

		# sort_list.append([chr, breakpoint, insert_size, final_type, str(final_MAPQ)])
		# print("%s\t%s\t%s\t%s"%(chr, breakpoint, insert_size, final_type))
	sort_list = sorted(sort_list, key = lambda x:(x[0], int(x[1])))
	file = open(out_path, 'w')
	logging.info("Writing results into disk...")
	for i in sort_list:
		# print("%s\t%s\t%s\t%s"%(i[0], breakpoint, insert_size, final_type))
		# print "\t".join(i)
		file.write("\t".join(i)+"\n")
	file.close()


def print_vcf_head():
	import time
	Date = time.strftime("%Y%m%d")
	head = list()
	head.append("##fileformat=VCFv4.2\n")
	head.append("##fileDate=%s\n"%(Date))
	head.append("##source=tjiang_scripts\n")
	head.append("##reference=Grch37\n")
	head.append("##ALT=<ID=<DEL>,Description=\"Deletion relative to the reference\">\n")
	head.append("##ALT=<ID=<INS>,Description=\"Insertion of sequence relative to the reference\">\n")
	head.append("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
	head.append("##INFO=<ID=SVLEN,Number=.,Type=String,Description=\"Difference in length between REF and ALT alleles\">\n")
	head.append("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count'\">\n")
	head.append("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele frequency'\">\n")
	head.append("##INFO=<ID=AN,Number=.,Type=String,Description=\"Allele name'\">\n")
	head.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	return head

def load_ref(ref_g):
	logging.info("Loading reference genome...")
	return SeqIO.to_dict(SeqIO.parse(ref_g, "fasta"))

def parse_seq_head(line):
	# seq = line.split('_')
	# Type = seq[0]
	# chr = seq[1]
	# breakpoint = seq[2]
	# size = seq[3]
	# GT = seq[5]+':'+seq[6]
	# local_info = R_INFO(Type, chr, breakpoint, size, GT)
	# return local_info

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

def call_vcf(args):
	path = args.input
	out_path = args.temp_dir + "calling.vcf"

	ref = load_ref(args.Reference)

	AlignmentFile = open(path, 'r')
	logging.info("Loading TE alignmets...")
	for line in AlignmentFile:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '@':
			continue

		local_info = parse_seq_head(seq[0])
		Flag = int(seq[1])
		sub_type = seq[2]
		MAPQ = int(seq[4])
		if flag_dic[Flag] != 0 and MAPQ >= args.min_mapq:
			# to do something
			key = "%s_%s_%s_%s"%(local_info.Chr, local_info.Pos, local_info.Len, local_info.GT)
			if key not in cluster_dic:
				cluster_dic[key] = list()
			cluster_dic[key].append("<%s:ME:%s>\t%d"%(local_info.Type, sub_type, flag_dic[Flag]))
	AlignmentFile.close()

	sort_list = list()
	for i in cluster_dic:
		chr, breakpoint, insert_size, GT = parse_name(i)
		final_type = acquire_count_max(cluster_dic[i]).split('\t')[0]
		# print acquire_count_max(cluster_dic[i])
		strand = STRAND[acquire_count_max(cluster_dic[i]).split('\t')[1]]
		# final_type = cluster_dic[i][1]
		# final_MAPQ = cluster_dic[i][0]
		# final_strand = final.split('&')[1]
		# final_strand = acquire_count_max(cluster_dic[i][1])
		# sort_list.append([chr, breakpoint, insert_size, final_type, strand_dic[int(final_strand)]])

		# if final_MAPQ >= 20:
		sort_list.append([chr, breakpoint, insert_size, final_type, GT, strand])

		# sort_list.append([chr, breakpoint, insert_size, final_type, str(final_MAPQ)])
		# print("%s\t%s\t%s\t%s"%(chr, breakpoint, insert_size, final_type))
	sort_list = sorted(sort_list, key = lambda x:(x[0], int(x[1])))
	head_info = print_vcf_head()

	file = open(out_path, 'w')
	logging.info("Writing results into disk...")

	for line in head_info:
		file.write(line)

	ID = 0
	for i in sort_list:
		# print("%s\t%s\t%s\t%s"%(i[0], breakpoint, insert_size, final_type))
		# print "\t".join(i)
		INFO = "SVTYPE=%s;SVLEN=%d;END=%d;SAMPLE=%s;STRAND=%s"%(i[3][1:4], int(i[2]), int(i[1])+int(i[2])-1, args.sample, i[5])
		concordant = int(i[4].split(':')[0])
		discordant = int(i[4].split(':')[1]) - int(i[4].split(':')[0])
		if discordant < 0:
			discordant = 0
		# print concordant, discordant
		# GT, GL = genotype_call_with_read_pair(concordant, discordant)
		# print GT, GL
		GT, GL = simple_call_genotype(concordant, concordant+discordant, args.heterozygous, args.homozygous)
		# print("%s\t%s\t%d\tN\t%s\t.\t.\t%s\tGT:DV:DR\t%s:%s"%(i[0], i[1], ID, i[3], INFO, GT, GL))
		try:
			REF = ref[i[0]][int(i[1])-1]
		except:
			REF = "N"
		file.write("%s\t%s\t%d\t%s\t%s\t.\t.\t%s\tGT:DV:DR\t%s:%s\n"%(i[0], i[1], ID, REF, i[3], INFO, GT, GL))
		ID += 1

	file.close()

VERSION="1.0"

USAGE="""\
	Non-reference TE calling and genotyping.

	Optional output format: .bed and .vcf
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="process.py call", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	# parser.add_argument("AlignmentFile", type=str, help="the bam format file generated by ngmlr, within a '.bai' index file")
	parser.add_argument("input", metavar="SAM", type=str, help="Input TE seqs to be mapped")
	parser.add_argument("Reference", metavar="REFERENCE", type=str, help="the reference genome(fasta format)")
	parser.add_argument("format", metavar="[BED,VCF]", type=str, help="The format of the output file. [%(default)s]", default = "bed")
	parser.add_argument('temp_dir', type=str, help = "temporary directory to use for distributed jobs")
	# parser.add_argument("--temp", type=str, default=tempfile.gettempdir(), help="Where to save temporary files")

	# parser.add_argument('-s', '--min_support', help = "Mininum number of reads that support a TE.[%(default)s]", default = 5, type = int)
	# parser.add_argument('-l', '--min_length', help = "Mininum length of TE to be reported.[%(default)s]", default = 50, type = int)
	# parser.add_argument('-d', '--min_distance', help = "Mininum distance of two TE clusters.[%(default)s]", default = 20, type = int)
	parser.add_argument('-hom', '--homozygous', help = "The mininum score of a genotyping reported as a homozygous.[%(default)s]", default = 0.8, type = float)
	parser.add_argument('-het','--heterozygous', help = "The mininum score of a genotyping reported as a heterozygous.[%(default)s]", default = 0.3, type = float)
	parser.add_argument('-q', '--min_mapq', help = "Mininum mapping quality.[%(default)s]", default = 20, type = int)
	parser.add_argument('--sample', help = "The name of the sample which be noted.", default = "None", type = str)
	# parser.add_argument('-t', '--threads', help = "Number of threads to use.[%(default)s]", default = 1, type = int)
	# parser.add_argument('-x', '--presets', help = "The sequencing type <pacbio,ont> of the reads.[%(default)s]", default = "pacbio", type = str)
	# parser.add_argument('--subread_length', help = "Length of fragments reads are split into [%(default)s]", default = 128, type = int)
	# parser.add_argument('--subread_corridor', help = "Length of corridor sub-reads are aligned with [%(default)s]", default = 20, type = int)
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

def run(argv):
    args = parseArgs(argv)
    setupLogging(False)
    starttime = time.time()
	
    if args.format == "bed":
    	call_bed(args)
    elif args.format == "vcf":
    	call_vcf(args)
    else:
    	logging.error("The format is available.")
    	exit(1)

    logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

if __name__ == '__main__':
	# path = sys.argv[2]
	# choice = sys.argv[1]
	# if choice == 'bed':
	# 	call_bed(path)
	# elif choice == 'vcf':
	# 	call_vcf(path)
	run(sys.argv[:1])
