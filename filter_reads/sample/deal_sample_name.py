import sys

def deal_10X_data(path):
	dic = {}
	file = open(path)
	for line in file:
		seq = line.strip()
		flag_1 = seq.split('.')
		if flag_1[-1] == "xml":
			key = seq.split('/')[2].split('.')[0]+"."+"fastq"
			if key not in dic:
				dic[key] = 1
	file.close()
	return dic

def deal_54X_data(path, dic):
	file = open(path)
	for line in file:
		seq = line.strip()
		outline = seq.split(':')[1].split(' ')[1].split('/')[-1]
		if outline in dic:
			print seq
	file.close()

def main():
	data = deal_10X_data(sys.argv[1])
	deal_54X_data(sys.argv[2], data)

if __name__ == '__main__':
	main()