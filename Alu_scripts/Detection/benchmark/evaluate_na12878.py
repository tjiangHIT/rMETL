
import sys
from data_collection import *

# dic_tea = dict()
# dic_tangram = dict()
# dic_retroseq = dict()
DICT = list()
standard = 50

def check_locus(chr, pos):
	dic_tea = DICT[0]
	dic_tangram = DICT[1]
	dic_retroseq = DICT[2]
	result = [0, 0, 0]
	if chr in dic_tea:
		for i in xrange(len(dic_tea[chr])):
			truth_pos = dic_tea[chr][i][0]
			if truth_pos - standard <= pos and pos <= truth_pos + standard:
				dic_tea[chr][i][3] = 1
				result[0] = 1
				# break
	# for key in dic_tangram:
	# 	print key, type(key), chr, type(chr)
	if chr in dic_tangram:
		# print chr
		for i in xrange(len(dic_tangram[chr])):
			truth_pos = dic_tangram[chr][i][0]
			# print truth_pos
			if truth_pos - standard <= pos and pos <= truth_pos + standard:
				dic_tangram[chr][i][3] = 1
				result[1] = 1
				# break
	if chr in dic_retroseq:
		for i in xrange(len(dic_retroseq[chr])):
			truth_pos = dic_retroseq[chr][i][0]
			if truth_pos - standard <= pos and pos <= truth_pos + standard:
				dic_retroseq[chr][i][3] = 1
				result[2] = 1
				# break
	return result

def call_overlap(local_l, key1, key2):
	chr = local_l[0]
	pos = local_l[1]
	overlap = 0
	if chr in DICT[key1]:
		for i in xrange(len(DICT[key1][chr])):
			if DICT[key1][chr][i][3] == 0:
				local_pos = DICT[key1][chr][i][0]
				if local_pos - standard <= pos and pos <= local_pos + standard:
					overlap += 1
	if chr in DICT[key2]:
		for i in xrange(len(DICT[key2][chr])):
			if DICT[key2][chr][i][3] == 0:
				local_pos = DICT[key2][chr][i][0]
				if local_pos - standard <= pos and pos <= local_pos + standard:
					overlap += 1
	return overlap

def evaluation(p, op):
	# r = check_locus('1', 5934294)
	# print r
	nice_job_1 = 0
	nice_job_2 = 0
	nice_job_3 = 0
	file = open(p, 'r')
	temp_chr = ''
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		pos = int(seq[1])
		length = int(seq[2])
		read_conut = int(seq[3])

		if temp_chr != chr:
			temp_chr = chr
			# print("[INFO]: Testing chromsome %s."%(temp_chr))
		# print chr, pos
		result = check_locus(chr, pos)
		if sum(result) == 1:
			nice_job_1 += 1
		if sum(result) == 2:
			nice_job_2 +=1
		if sum(result) == 3:
			nice_job_3 +=1
	file.close()
	print("[INFO]: (1,2,3) %d,%d,%d signal locuses predictions are good!"%(nice_job_1, nice_job_2, nice_job_3))

	out_file = open(op, 'w')
	Tp = 0
	Np = 0
	for key in DICT[0]:
		for i in xrange(len(DICT[0][key])):
			if DICT[0][key][i][3] == 0:
				Np += 1
				overlap = call_overlap(DICT[0][key][i], 1, 2)
				out_file.write("Tea\t%s\t%d\t%s\t%d\n"%(key, DICT[0][key][i][0], DICT[0][key][i][2], overlap))
			else:
				Tp += 1
	print("[INFO]: For Tea callset: Y for %d, N for %d(%0.5f)."%(Tp, Np, Tp*1.0/(Tp+Np)))

	Tp = 0
	Np = 0
	for key in DICT[1]:
		for i in xrange(len(DICT[1][key])):
			if DICT[1][key][i][3] == 0:
				Np += 1
				overlap = call_overlap(DICT[1][key][i], 0, 2)
				out_file.write("Tangram\t%s\t%d\t%s\t%d\n"%(key, DICT[1][key][i][0], DICT[1][key][i][2], overlap))
			else:
				Tp += 1
	print("[INFO]: For Tangram callset: Y for %d, N for %d(%0.5f)."%(Tp, Np, Tp*1.0/(Tp+Np)))

	Tp = 0
	Np = 0
	for key in DICT[2]:
		for i in xrange(len(DICT[2][key])):
			if DICT[2][key][i][3] == 0:
				Np += 1
				overlap = call_overlap(DICT[2][key][i], 0, 1)
				out_file.write("RetroSeq\t%s\t%d\t%s\t%d\n"%(key, DICT[2][key][i][0], DICT[2][key][i][2], overlap))
			else:
				Tp += 1
	print("[INFO]: For RetroSeq callset: Y for %d, N for %d(%0.5f)."%(Tp, Np, Tp*1.0/(Tp+Np)))
	out_file.close()


def main():
	import time
	starttime = time.time()
	Tea_path = sys.argv[1]
	Tea_path_2 = sys.argv[2]
	Tangram_path = sys.argv[3]
	RetroSeq_path = sys.argv[4]
	print("[INFO]: The path of the Tea(Alu) callset: %s"%(Tea_path))
	print("[INFO]: The path of the Tea(L1) callset: %s"%(Tea_path_2))
	print("[INFO]: The path of the Tangram callset: %s"%(Tangram_path))
	print("[INFO]: The path of the RetroSeq callset: %s"%(RetroSeq_path))
	dic_tea = collect_Tea(Tea_path, Tea_path_2)
	dic_tangram = collect_Tangram(Tangram_path)
	dic_retroseq = collect_RetroSeq(RetroSeq_path)
	DICT.append(dic_tea)
	DICT.append(dic_tangram)
	DICT.append(dic_retroseq)
	# print len(dic_tea)
	# print len(dic_tangram)
	# print len(dic_retroseq)

	signal_path = sys.argv[5]
	out_path = sys.argv[6]
	print("[INFO]: The path of the Signal file: %s"%(signal_path))
	print("[INFO]: Starting to evaluate...")
	evaluation(signal_path, out_path)

	print("[INFO]: Finished in %0.2f seconds."%(time.time() - starttime))
	print("[INFO]: Have a nice day! ^.^ ")

if __name__ == '__main__':
# 	collect()
	main()