'''
usage:
	python2 (dbRIP)ALU-L1-SVA > filter.txt
	filter.txt: 1. target sequence duplication
				2. element sequence
'''

import sys

def main():
	path_1 = sys.argv[1]
	path_2 = sys.argv[2]
	file = open(path_1, 'r')
	flag = dict()
	for line in file:
		seq = line.strip('\n').split('\t')[0]
		if seq not in flag:
			flag[seq] = 0
	file.close()

	file = open(path_2, 'r')
	num = 0
	for line in file:
		num += 1
		seq = line.strip('\n')
		if str(num) in flag:
			print(seq)		
	file.close()	

if __name__ == '__main__':
	main()