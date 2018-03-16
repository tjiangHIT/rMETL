import sys
from data_collection import *

def main():
	import time
	starttime = time.time()
	path_in = sys.argv[1]
	print("[INFO]: The path of the 1000G File: %s"%(path_in))
	path_out = sys.argv[2]
	print("[INFO]: The path of the output File: %s"%(path_out))
	data = collect_1000G(path_in)
	file = open(path_out, 'w')
	for i in data:
		line = "%s\t%d\t%d\t%s\t%s\n"%(i[0], i[1], i[2][0], i[2][1], i[2][2])
		file.write(line)
	file.close()
	print("[INFO]: Finished in %0.2f seconds."%(time.time() - starttime))
	print("[INFO]: Have a nice day! ^.^ ")

if __name__ == '__main__':
	main()