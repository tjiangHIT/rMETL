import argparse

def opt():
	# parser = argparse.ArgumentParser(description = 'Detect non-reference Transcription Elements.', usage = '%(prog)s [options]')
	parser = argparse.ArgumentParser(description = 'Detect non-reference Transcription Elements.')
	parser.add_argument('AlignmentFile', help = "The bam format file established by ngmlr, within a '.bai' index file.")
	parser.add_argument('Output_prefix', help = "The prefix of the output files.")
	parser.add_argument('Reference', help = "The reference genome(fasta format).")
	parser.add_argument('-s', '--min_support', help = "Mininum number of reads that support a TE.[5]", default = 5, type = int)
	parser.add_argument('-l', '--min_length', help = "Mininum length of TE to be reported.[50]", default = 50, type = int)
	parser.add_argument('-d', '--min_distance', help = "Mininum distance of two TE clusters.[20]", default = 20, type = int)
	parser.add_argument('-hom', '--homozygous', help = "The mininum score of a genotyping reported as a homozygous.[0.8]", default = 0.8, type = float)
	parser.add_argument('-het','--heterozygous', help = "The mininum score of a genotyping reported as a heterozygous.[0.3]", default = 0.3, type = float)
	parser.add_argument('-q', '--min_mapq', help = "Mininum mapping quality.[20]", default = 20, type = int)

	parser.add_argument('-t', '--threads', help = "Number of threads to use.[1]", default = 1, type = int)

	args = parser.parse_args()
	return args

# if __name__ == '__main__':
# 	args = opt()
# 	print args.heterozygous