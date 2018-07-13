import sys
import random
import time

def main():
	dic = dict()

	file = open(sys.argv[1], 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		if chr != 'chr1':
			continue
		subtype = seq[3]

		if subtype not in dic:
			# dic[subtype] = list()
			dic[subtype] = 0

		# dic[subtype].append()
		dic[subtype] += 1
	file.close()

	sumtotal = 0
	for key in dic:
		# print key, dic[key]
		if key[:3] in ['Alu', 'SVA'] or key[:2] == 'L1':
			continue
		sumtotal += dic[key]

	print sumtotal

if __name__ == '__main__':
	main()