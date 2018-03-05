import sys
import pysam
# import tools

CLIP_note = dict()

def store_clip_pos(locus, chr, length, alu_type):
	# about collecting breakpoint from clipping 
	hash_1 = int(locus /10000)
	mod = locus % 10000
	hash_2 = int(mod / 50)

	if hash_1 not in CLIP_note[chr]:
		CLIP_note[chr][hash_1] = dict()
		CLIP_note[chr][hash_1][hash_2] = list()
		CLIP_note[chr][hash_1][hash_2].append([locus, length, alu_type])
	else:
		if hash_2 not in CLIP_note[chr][hash_1]:
			CLIP_note[chr][hash_1][hash_2] = list()
			CLIP_note[chr][hash_1][hash_2].append([locus, length, alu_type])
		else:
			CLIP_note[chr][hash_1][hash_2].append([locus, length, alu_type])

def find_common(chr, start_pos, length):
	answer = list()
	up = start_pos + length
	down = start_pos
	# print down, up
	for k in xrange(int(up/10000) - int(down/10000) + 1):
		key_1 = int(down/10000) + k
		# print key_1
		if key_1 not in CLIP_note[chr]:
			continue
			# return answer
		for i in xrange(int((up%10000)/50)-int((down%10000)/50)+1):
			key_2 = int((down%10000)/50)+i
			# print key_2
			if key_2 not in CLIP_note[chr][key_1]:
				continue
				# return answer
			for ele in CLIP_note[chr][key_1][key_2]:
				this_s = ele[0]
				this_e = ele[0] + ele[1]
				# print this_s, this_e
				if (this_s <= down and this_e >= down) or (this_s <= up and this_e >= up) or (this_s >= down and this_e <= up):
					answer.append(ele)
	return answer

def compare(p1, p2, p3):
	count_1 = 0
	count_2 = 0
	total = 0
	truth_file = open(p1, 'r')
	print("[INFO]: Loading the grand truth...")
	for line in truth_file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		start_pos = int(seq[1])
		length = int(seq[2])
		alu_type = seq[3]
		if chr not in CLIP_note:
			CLIP_note[chr] = dict()
		store_clip_pos(start_pos, chr, length, alu_type)
		total += 1

	truth_file.close()

	predict_file = open(p2, 'r')
	print("[INFO]: Processing the scripts' prediction...")
	for line in predict_file:
		seq = line.strip('\n').split('\t')
		# print seq
		chr = seq[0]
		start_pos = int(seq[1])
		length = int(seq[2])
		read_count = int(seq[3])
		truth = find_common(chr, start_pos, length)
		if len(truth):
			count_1 += 1
			# print truth
		# break
	predict_file.close()

	sniffles_file = open(p3, 'r')
	print("[INFO]: Processing the Sniffles's prediction...")
	for line in sniffles_file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		start_pos = int(seq[1])
		end_pos = int(seq[2])
		if start_pos != end_pos:
			continue
		length = int(seq[-1])
		truth = find_common(chr, start_pos, length)
		if len(truth):
			count_2 += 1

	sniffles_file.close()
	return count_1, count_2, total


def main():
	import time
	starttime = time.time()
	truth_path = sys.argv[1]
	predict_path = sys.argv[2]
	sniffles_path = sys.argv[3]
	# pre_out = sys.argv[4]

	print("[INFO]: 1. The path of the grand truth: %s"%(truth_path))
	print("[INFO]: 2. The path of the scripts' predictions: %s"%(predict_path))
	print("[INFO]: 3. The path of the Sniffles' answer: %s"%(sniffles_path))
	# print("[INFO]: 4. The prefix path of output files: %s"%(pre_out))
	# print("[INFO]: This script is uesd for simulate Alus on the reference genome.")

	count_1, count_2, total = compare(truth_path, predict_path, sniffles_path)
	# print count_1, count_2
	print("[INFO]: The total number of the grand truth: %d"%(total))
	print("[INFO]: The total number of the scripts' predictions: %d"%(count_1))
	print("[INFO]: The recall rate is %0.3f"%(count_1*1.0/total))
	print("[INFO]: The total number of the Sniffles's predictions: %d"%(count_2))
	print("[INFO]: The recall rate is %0.3f"%(count_2*1.0/total))
	print("[INFO]: Finished in %0.2f seconds."%(time.time() - starttime))
	print("[INFO]: Have a nice day! ^.^ ")

if __name__ == '__main__':
	main()