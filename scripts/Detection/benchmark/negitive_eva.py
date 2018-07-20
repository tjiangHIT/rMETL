
import sys

def load_rMETL(path, dic):
	file = open(path, 'r')
	num = 0
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		pos = int(seq[1])
		lengtn = int(seq[2])
		subtype = seq[3][1]

		flag = 0
		if chr in dic:
			for i in dic[chr]:
				if i[2] == subtype:
					if (pos <= i[0] and pos+ lengtn >= i[0]) or (pos <= i[1] and pos+ lengtn >= i[1]) or (pos >= i[0] and pos+ lengtn <= i[1]):
						flag = 1
				# if i[2] == subtype:
				# 	if pos <= i[0] + 50 and pos >= i[0] - 50:
				# 		flag = 1
		if flag == 1:
			# print line.strip('\n')
			num	+= 1
	print num
	file.close()

def load_sniffles(path):
	dic = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		start = int(seq[13])
		end = int(seq[15])
		# end = int(seq[15])
		subtype = seq[10][:3]
		if subtype == "DEL":
			subtype = 'D'
		else:
			subtype = 'I'

		if chr not in dic:
			dic[chr] = list()
		dic[chr].append([start, end, subtype])
	file.close()
	return dic

def main():
	# print "ASD"
	sniffles = sys.argv[1]
	dic = load_sniffles(sniffles)

	# for i in dic['1']:
	# 	print i
	rMETL = sys.argv[2]
	load_rMETL(rMETL, dic)

if __name__ == '__main__':
	main()
