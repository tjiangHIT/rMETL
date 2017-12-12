import sys

def init_dic():
	data = []
	for i in xrange(25):
		# dic = {}
		dic = []
		data.append(dic)
	return data

def get_chr_no(Chr_str):
	code = Chr_str[3:]
	if code.isdigit() == True:
		return int(code)-1
	else:
		if code == 'X':
			return 22
		elif code == 'Y':
			return 23
		else:
			return 24


def acquire_info(path):
	# data = []
	SV_TYPE = {}
	file = open(path)
	for line in file:
		if line[0] != '#':
			info = []
			SVtype = line.strip().split('\t')[10]
			best_chr1 = line.strip().split('\t')[12]
			best_start = line.strip().split('\t')[13]
			best_chr2 = line.strip().split('\t')[14]
			best_stop = line.strip().split('\t')[15]
			# info.append(SVtype)
			# info.append(best_chr1)
			# info.append(best_start)
			# info.append(best_chr2)
			# info.append(best_stop)
			# data.append(info)

			if abs(int(best_stop)-int(best_start)) < 50:
				continue

			if SVtype not in SV_TYPE:
				SV_TYPE[SVtype] = init_dic()

			if get_chr_no(best_chr1) == get_chr_no(best_chr2):
				info.append(best_start)
				info.append(best_stop)
				# info.append(int(best_stop)-int(best_start))
				SV_TYPE[SVtype][get_chr_no(best_chr1)].append(info)
			else:
				info.append(best_chr1)
				info.append(best_start)
				info.append(best_chr2)
				info.append(best_stop)
				SV_TYPE[SVtype][get_chr_no(best_chr1)].append(info)

	file.close()
	return SV_TYPE

# def complete_dic(List, data):
# 	for i in xrange(len(data)):
# 		CHR = 


def main():
	# LIST = init_dic()
	# print LIST
	data_1 = acquire_info(sys.argv[1])
	data_2 = acquire_info(sys.argv[2])

	for key in data_1:
		for i in xrange(25):
			# print key, "chr"+str(i+1)+":", len(data_1[key][i]), len(data_2[key][i])
			if key == "DUP":
				print len(data_2[key][i])
	# cmp scripts
	

if __name__ == '__main__':
	main()