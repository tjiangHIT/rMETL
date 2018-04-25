
import sys
from data_collection import *

# dic_tea = dict()
# dic_tangram = dict()
# dic_retroseq = dict()
DICT = list()
standard = 20

final_statics = list()

def check_locus(chr, pos, flag):
	dic_tea = DICT[0]
	dic_tangram = DICT[1]
	dic_retroseq = DICT[2]
	result = [0, 0, 0]
	if chr in dic_tea:
		for i in xrange(len(dic_tea[chr])):
			truth_pos = dic_tea[chr][i][0]
			if truth_pos - standard <= pos and pos <= truth_pos + standard:
				DICT[0][chr][i][3] = flag
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
				DICT[1][chr][i][3] = flag
				result[1] = 1
				# break
	if chr in dic_retroseq:
		for i in xrange(len(dic_retroseq[chr])):
			truth_pos = dic_retroseq[chr][i][0]
			if truth_pos - standard <= pos and pos <= truth_pos + standard:
				DICT[2][chr][i][3] = flag
				result[2] = 1
				# break
	return result

def call_overlap(local_l, key1, key2, flag):
	chr = local_l[0]
	pos = local_l[1]
	overlap = 0
	if chr in DICT[key1]:
		for i in xrange(len(DICT[key1][chr])):
			if DICT[key1][chr][i][3] != flag:
				local_pos = DICT[key1][chr][i][0]
				if local_pos - standard <= pos and pos <= local_pos + standard:
					overlap += 1
	if chr in DICT[key2]:
		for i in xrange(len(DICT[key2][chr])):
			if DICT[key2][chr][i][3] != flag:
				local_pos = DICT[key2][chr][i][0]
				if local_pos - standard <= pos and pos <= local_pos + standard:
					overlap += 1
	return overlap

def check_1kg(chr, pos, flag):
	result = 0
	dic_1kg = DICT[3]
	if chr in dic_1kg:
		for i in xrange(len(dic_1kg[chr])):
			truth_pos = dic_1kg[chr][i][0]
			if truth_pos - standard <= pos and pos <= truth_pos + standard:
				DICT[3][chr][i][1] = flag
				result = 1
	return result

def build_statics():
	print("\n[INFO]: Threshold is %d"%(standard))
	# print("*"*50)
	print("[INFO]: Tea\tTangram\tRetroSeq\t1KG")
	print("[INFO]: %d/%d\t%d/%d\t%d/%d\t%d/%d"%(final_statics[0][0], final_statics[4][0], final_statics[1][0], final_statics[5][0], final_statics[2][0], final_statics[6][0], final_statics[3][0], final_statics[7][0]))
	print("[INFO]: %d/%d\t%d/%d\t%d/%d\t%d/%d"%(final_statics[0][1], final_statics[4][1], final_statics[1][1], final_statics[5][1], final_statics[2][1], final_statics[6][1], final_statics[3][1], final_statics[7][1]))
	print("[INFO]: %0.5f/%0.5f\t%0.5f/%0.5f\t%0.5f/%0.5f\t%0.5f/%0.5f\n"%(final_statics[0][2], final_statics[4][2], final_statics[1][2], final_statics[5][2], final_statics[2][2], final_statics[6][2], final_statics[3][2], final_statics[7][2]))
	# print("*")
	# print("*"*50)

def evaluation(p, op, p2):
	# r = check_locus('1', 5934294)
	# print r
	nice_job_1 = 0
	nice_job_2 = 0
	nice_job_3 = 0
	nice_job_4 = 0
	file = open(p, 'r')
	temp_chr = ''
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		pos = int(seq[1])
		length = int(seq[2])
		# read_conut = int(seq[3])

		if temp_chr != chr:
			temp_chr = chr
			# print("[INFO]: Testing chromsome %s."%(temp_chr))
		# print chr, pos
		result = check_locus(chr, pos, 1)
		if sum(result) == 1:
			nice_job_1 += 1
		if sum(result) == 2:
			nice_job_2 +=1
		if sum(result) == 3:
			nice_job_3 +=1

		result = check_1kg(chr, pos, 1)
		if result == 1:
			nice_job_4 += 1

	file.close()
	print("[INFO]: (1,2,3,4) %d,%d,%d,%d signal locuses predictions are good!"%(nice_job_1, nice_job_2, nice_job_3, nice_job_4))

	out_file = open(op, 'w')
	Tp = 0
	Np = 0
	for key in DICT[0]:
		for i in xrange(len(DICT[0][key])):
			if DICT[0][key][i][3] == 0:
				Np += 1
				overlap = call_overlap(DICT[0][key][i], 1, 2, 1)
				out_file.write("Tea\t%s\t%d\t%s\t%d\n"%(key, DICT[0][key][i][0], DICT[0][key][i][2], overlap))
			else:
				Tp += 1
	print("[INFO]: For Tea callset: Y for %d, N for %d(%0.5f)."%(Tp, Np, Tp*1.0/(Tp+Np)))
	final_statics.append([Tp, Np, Tp*1.0/(Tp+Np)])

	Tp = 0
	Np = 0
	for key in DICT[1]:
		for i in xrange(len(DICT[1][key])):
			if DICT[1][key][i][3] == 0:
				Np += 1
				overlap = call_overlap(DICT[1][key][i], 0, 2, 1)
				out_file.write("Tangram\t%s\t%d\t%s\t%d\n"%(key, DICT[1][key][i][0], DICT[1][key][i][2], overlap))
			else:
				Tp += 1
	print("[INFO]: For Tangram callset: Y for %d, N for %d(%0.5f)."%(Tp, Np, Tp*1.0/(Tp+Np)))
	final_statics.append([Tp, Np, Tp*1.0/(Tp+Np)])

	Tp = 0
	Np = 0
	for key in DICT[2]:
		for i in xrange(len(DICT[2][key])):
			if DICT[2][key][i][3] == 0:
				Np += 1
				overlap = call_overlap(DICT[2][key][i], 0, 1, 1)
				out_file.write("RetroSeq\t%s\t%d\t%s\t%d\n"%(key, DICT[2][key][i][0], DICT[2][key][i][2], overlap))
			else:
				Tp += 1
	print("[INFO]: For RetroSeq callset: Y for %d, N for %d(%0.5f)."%(Tp, Np, Tp*1.0/(Tp+Np)))
	final_statics.append([Tp, Np, Tp*1.0/(Tp+Np)])
	out_file.close()
	Tp = 0
	Np = 0
	for key in DICT[3]:
		for i in xrange(len(DICT[3][key])):
			if DICT[3][key][i][1] == 0:
				Np += 1
				# overlap = call_overlap(DICT[0][key][i], 1, 2, 1)
				# out_file.write("Tea\t%s\t%d\t%s\t%d\n"%(key, DICT[0][key][i][0], DICT[0][key][i][2], overlap))
			else:
				Tp += 1
	print("[INFO]: For 1kg callset: Y for %d, N for %d(%0.5f)."%(Tp, Np, Tp*1.0/(Tp+Np)))
	final_statics.append([Tp, Np, Tp*1.0/(Tp+Np)])

	nice_job_1 = 0
	nice_job_2 = 0
	nice_job_3 = 0
	nice_job_4 = 0
	file = open(p2, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		subtype = seq[10]
		start_pos = int(seq[1])
		# end_pos = int(seq[15])
		if subtype == 'INS' or subtype == "DUP":
			result = check_locus(chr, start_pos, 2)
			if sum(result) == 1:
				nice_job_1 += 1
			if sum(result) == 2:
				nice_job_2 +=1
			if sum(result) == 3:
				nice_job_3 +=1

			result = check_1kg(chr, start_pos, 2)
			if result == 1:
				nice_job_4 += 1
	file.close()
	print("[INFO]: (1,2,3,4) %d,%d,%d,%d signal locuses Sniffles are good!"%(nice_job_1, nice_job_2, nice_job_3, nice_job_4))

	Tp = 0
	Np = 0
	for key in DICT[0]:
		for i in xrange(len(DICT[0][key])):
			if DICT[0][key][i][3] != 2:
				Np += 1
				overlap = call_overlap(DICT[0][key][i], 1, 2, 2)
				# out_file.write("Tea\t%s\t%d\t%s\t%d\n"%(key, DICT[0][key][i][0], DICT[0][key][i][2], overlap))
			else:
				Tp += 1
	print("[INFO]: For Tea callset: Y for %d, N for %d(%0.5f)."%(Tp, Np, Tp*1.0/(Tp+Np)))
	final_statics.append([Tp, Np, Tp*1.0/(Tp+Np)])

	Tp = 0
	Np = 0
	for key in DICT[1]:
		for i in xrange(len(DICT[1][key])):
			if DICT[1][key][i][3] != 2:
				Np += 1
				overlap = call_overlap(DICT[1][key][i], 0, 2, 2)
				# out_file.write("Tangram\t%s\t%d\t%s\t%d\n"%(key, DICT[1][key][i][0], DICT[1][key][i][2], overlap))
			else:
				Tp += 1
	print("[INFO]: For Tangram callset: Y for %d, N for %d(%0.5f)."%(Tp, Np, Tp*1.0/(Tp+Np)))
	final_statics.append([Tp, Np, Tp*1.0/(Tp+Np)])

	Tp = 0
	Np = 0
	for key in DICT[2]:
		for i in xrange(len(DICT[2][key])):
			if DICT[2][key][i][3] != 2:
				Np += 1
				overlap = call_overlap(DICT[2][key][i], 0, 1, 2)
				# out_file.write("RetroSeq\t%s\t%d\t%s\t%d\n"%(key, DICT[2][key][i][0], DICT[2][key][i][2], overlap))
			else:
				Tp += 1
	print("[INFO]: For RetroSeq callset: Y for %d, N for %d(%0.5f)."%(Tp, Np, Tp*1.0/(Tp+Np)))
	final_statics.append([Tp, Np, Tp*1.0/(Tp+Np)])

	Tp = 0
	Np = 0
	for key in DICT[3]:
		for i in xrange(len(DICT[3][key])):
			if DICT[3][key][i][1] != 2:
				Np += 1
				# overlap = call_overlap(DICT[0][key][i], 1, 2, 1)
				# out_file.write("Tea\t%s\t%d\t%s\t%d\n"%(key, DICT[0][key][i][0], DICT[0][key][i][2], overlap))
			else:
				Tp += 1
	print("[INFO]: For 1kg callset: Y for %d, N for %d(%0.5f)."%(Tp, Np, Tp*1.0/(Tp+Np)))
	final_statics.append([Tp, Np, Tp*1.0/(Tp+Np)])
	build_statics()


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
	out_path = sys.argv[7]
	sniffles_path = sys.argv[6]
	print("[INFO]: The path of the Signal file: %s"%(signal_path))
	print("[INFO]: The path of the Sniffles file: %s"%(sniffles_path))

	na1k_path = sys.argv[8]
	print("[INFO]: The path of the 1KG file: %s"%(na1k_path))
	dic_1kg = collect_1kg(na1k_path)
	DICT.append(dic_1kg)

	print("[INFO]: Starting to evaluate...")
	# print("[INFO]: ====Tea/Tangram/RetroSeq Part====")
	evaluation(signal_path, out_path, sniffles_path)

	# print("[INFO]: ======== 1000 Genome Part========")

	print("[INFO]: Finished in %0.2f seconds."%(time.time() - starttime))
	print("[INFO]: Have a nice day! ^.^ ")

if __name__ == '__main__':
# 	collect()
	main()