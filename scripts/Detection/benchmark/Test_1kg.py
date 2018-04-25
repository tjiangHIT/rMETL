import sys
from data_collection import *

data_dic = dict()

def convert_data(_list_):
	for i in _list_:
		chr = i[0]
		pos = i[1]
		length = i[2][0]
		sub_type = i[2][1]
		strand = i[2][2]

		if chr not in data_dic:
			data_dic[chr] = list()
		data_dic[chr].append([pos, length, sub_type, strand, 0])
	# return data_dic

standard = 10
total_whole_right = 0
total_right = 0

def check_locus(_list_):
	chr = _list_[0]
	pos = int(_list_[1])
	length = int(_list_[2])
	sub_type = _list_[3]
	# strand = _list_[4]
	# print pos
	if chr in data_dic and length >= 50:
		for i in xrange(len(data_dic[chr])):
			truth_pos = data_dic[chr][i][0]
			# print pos, truth_pos - standard, truth_pos + standard
			if truth_pos - standard <= pos and pos <= truth_pos + standard:
				# print pos, truth_pos
				# print dic[chr][i][2], dic[chr][i][3], dic[chr][i][1], sub_type, strand, length
				# global total_right
				# total_right += 1
				data_dic[chr][i][4] = 1
				if sub_type[0] == data_dic[chr][i][2][0]:
					data_dic[chr][i][4] = 2
					# global total_whole_right
					# total_whole_right += 1
				# break


def evaluation(p):
	file = open(p, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		# chr = seq[0]
		# pos = int(seq[1])
		# length = int(seq[2])
		# sub_type = seq[3]
		# strand = seq[4]
		check_locus(seq)
	total_0 = 0
	total_1 = 0
	total_2 = 0
	for i in data_dic:
		for j in data_dic[i]:
			if j[4] == 1:
				total_1 += 1
			elif j[4] == 2:
				total_2 += 1
			else:
				total_0 += 1
	print "%d\t%d\t%d\t%0.5f\t%0.5f"%(total_0, total_1, total_2, (total_1 + total_2)*1.0/(total_0+total_1+total_2), total_2*1.0/(total_1+total_2))

	file.close()

def main():
	na1k_path = sys.argv[1]
	na12878_list = collect_1000G(na1k_path)
	convert_data(na12878_list)
	callset_path = sys.argv[2]
	evaluation(callset_path)
	# print total_right
	# print total_whole_right

if __name__ == '__main__':
	main()