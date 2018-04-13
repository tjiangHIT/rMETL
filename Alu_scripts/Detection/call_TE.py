import sys	
import pysam
from collections import Counter
from genotype import *

def acquire_count_max(_list_):
	c = Counter(_list_)
	return c.most_common(1)[0][0]

flag_dic = {0:1, 16:2, 2048:0, 2064:0, 4:0}

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
	seq = line.split('_')
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

def call_bed(p1):
	# samfile = pysam.AlignmentFile(p1)
	AlignmentFile = open(p1, 'r')
	for line in AlignmentFile:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '@':
			continue

		local_info = parse_name_tp(seq[0])
		Flag = int(seq[1])
		sub_type = seq[2]
		MAPQ = int(seq[4])
		if flag_dic[Flag] != 0 and MAPQ >= 20:
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
	for i in sort_list:
		# print("%s\t%s\t%s\t%s"%(i[0], breakpoint, insert_size, final_type))
		print "\t".join(i)

def print_vcf_head():
	import time
	Date = time.strftime("%Y%m%d")
	print("##fileformat=VCFv4.2")
	print("##fileDate=%s"%(Date))
	print("##source=tjiang_scripts")
	print("##reference=Grch37")
	print("##ALT=<ID=<DEL>,Description=\"Deletion relative to the reference\">")
	print("##ALT=<ID=<INS>,Description=\"Insertion of sequence relative to the reference\">")
	print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">")
	print("##INFO=<ID=SVLEN,Number=.,Type=String,Description=\"Difference in length between REF and ALT alleles\">")
	print("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count'\">")
	print("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele frequency'\">")
	print("##INFO=<ID=AN,Number=.,Type=String,Description=\"Allele name'\">")
	print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")


def parse_seq_head(line):
	# seq = line.split('_')
	# Type = seq[0]
	# chr = seq[1]
	# breakpoint = seq[2]
	# size = seq[3]
	# GT = seq[5]+':'+seq[6]
	# local_info = R_INFO(Type, chr, breakpoint, size, GT)
	# return local_info

	seq = line.split('_')
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

def call_vcf(p):
	AlignmentFile = open(p, 'r')
	for line in AlignmentFile:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '@':
			continue

		local_info = parse_seq_head(seq[0])
		Flag = int(seq[1])
		sub_type = seq[2]
		MAPQ = int(seq[4])
		if flag_dic[Flag] != 0 and MAPQ >= 20:
			# to do something
			key = "%s_%s_%s_%s"%(local_info.Chr, local_info.Pos, local_info.Len, local_info.GT)
			if key not in cluster_dic:
				cluster_dic[key] = list()
			cluster_dic[key].append("<%s:ME:%s>"%(local_info.Type, sub_type))
	AlignmentFile.close()

	sort_list = list()
	for i in cluster_dic:
		chr, breakpoint, insert_size, GT = parse_name(i)
		final_type = acquire_count_max(cluster_dic[i])
		# final_type = cluster_dic[i][1]
		# final_MAPQ = cluster_dic[i][0]
		# final_strand = final.split('&')[1]
		# final_strand = acquire_count_max(cluster_dic[i][1])
		# sort_list.append([chr, breakpoint, insert_size, final_type, strand_dic[int(final_strand)]])

		# if final_MAPQ >= 20:
		sort_list.append([chr, breakpoint, insert_size, final_type, GT])

		# sort_list.append([chr, breakpoint, insert_size, final_type, str(final_MAPQ)])
		# print("%s\t%s\t%s\t%s"%(chr, breakpoint, insert_size, final_type))
	sort_list = sorted(sort_list, key = lambda x:(x[0], int(x[1])))
	print_vcf_head()
	ID = 0
	for i in sort_list:
		# print("%s\t%s\t%s\t%s"%(i[0], breakpoint, insert_size, final_type))
		# print "\t".join(i)
		INFO = "SVTYPE=%s;SVLEN=%d;END=%d;SAMPLE=NA12878;STRAND=%s"%(i[3][1:4], int(i[2]), int(i[1])+int(i[2])-1, '+/-')
		concordant = int(i[4].split(':')[0])
		discordant = int(i[4].split(':')[1]) - int(i[4].split(':')[0])
		if discordant < 0:
			discordant = 0
		# print concordant, discordant
		# GT, GL = genotype_call_with_read_pair(concordant, discordant)
		# print GT, GL
		GT, GL = simple_call_genotype(concordant, concordant+discordant, 0.3, 0.8)
		print("%s\t%s\t%d\tN\t%s\t.\t.\t%s\tGT:DV:DR\t%s:%s"%(i[0], i[1], ID, i[3], INFO, GT, GL))
		ID += 1

def parseArgs(argv):
	pass

def run(argv):
    args = parseArgs(argv)


if __name__ == '__main__':
	path = sys.argv[2]
	choice = sys.argv[1]
	if choice == 'bed':
		call_bed(path)
	elif choice == 'vcf':
		call_vcf(path)
