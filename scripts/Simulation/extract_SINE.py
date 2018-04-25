'''
usage:
	python2 (dbRIP)ALU-L1-SVA > filter.txt
	filter.txt: 1. target sequence duplication
				2. element sequence
'''

import sys

def judge_str(seq):
	if seq.find('n') == -1 and seq.find('N') == -1:
		return 1
	else:
		return 0

def acquire_TSD_realELE(seq):
	split_TSD_signal = "<font color=\"green\"><u>"
	split_realELE_signal = "<font color=\"red\">"
	seq = seq.split(split_TSD_signal)
	TSD = "".join(seq[2].split('</u>')[0].split("<br>"))
	realELE = "".join(seq[1].split(split_realELE_signal)[1].split("</font>")[0].split("<br>"))
	return TSD, realELE

def main():
	file = open(sys.argv[1], 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		TSD, realELE = acquire_TSD_realELE(seq[13])
		# seq_2 = seq[13].split(flag_1)
		# for i in seq_2:
		# 	print(i)
		# if TSD.find('n') == -1 and realELE.find('N') == -1:

		# for SVA and L1
		if judge_str(TSD) == 1 and judge_str(realELE) == 1:
			# if len(TSD) <= 25:
			print(TSD+"\t"+realELE)
				'''
		for ALU:
		if judge_str(TSD) == 1 and judge_str(realELE) == 1:
			if len(TSD) <= 25:
				print(TSD+"\t"+realELE)
		'''
		# break
	file.close()	

if __name__ == '__main__':
	main()