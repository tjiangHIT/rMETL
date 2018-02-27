import sys
import pysam
import tools

def main():
	sam_path = sys.argv[1]
	tools.load_sam(sam_path)

if __name__ == '__main__':
	main()