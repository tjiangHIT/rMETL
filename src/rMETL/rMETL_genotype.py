# * @author: Jiang Tao (tjiang@hit.edu.cn)

GL_TAG = ['1/0', '0/1', '1/1']

def simple_call_genotype(Nalt, Ntotal, P_heterozygous, P_homozygous):
	bound_low = Ntotal * P_heterozygous
	bound_up = Ntotal * P_homozygous
	if Nalt < bound_low:
        # reliability = 0
		return GL_TAG[0], "%d:%d"%(Nalt, Ntotal - Nalt), 0
	elif bound_low <= Nalt and Nalt < bound_up:
        # reliability = 1
		return GL_TAG[1], "%d:%d"%(Nalt, Ntotal - Nalt), 1
	else:
        # reliability = 1
		return GL_TAG[2], "%d:%d"%(Nalt, Ntotal - Nalt), 1

def simple_filter_genotype(Nalt, Ntotal, P_heterozygous):
    bound_low = Ntotal * P_heterozygous
    if Nalt < bound_low:
        return 0
    else:
        return 1

def count_coverage(chr, s, e, f):
    total = 0
    for i in f.fetch(chr, s, e):
        total += 1
    return total

def add_genotype(info_list, file, low_bandary):
    '''
    allocate genotype for each MEI/MED
    '''
    for i in xrange(len(info_list)):
        if info_list[i][0][0] == 'INS':
            chr = info_list[i][0][1]
            start = info_list[i][0][2]-low_bandary
            end = info_list[i][0][2] + low_bandary
            locus_cov = count_coverage(chr, start, end, file)
            for j in xrange(len(info_list[i])):
                info_list[i][j].append(locus_cov)
        else:
            for j in xrange(len(info_list[i])):
                chr = info_list[i][j][1]
                start = info_list[i][j][2]
                end = info_list[i][j][2]+info_list[i][j][3]
                locus_cov = count_coverage(chr, start, end, file)
                info_list[i][j].append(locus_cov)
    return info_list

if __name__ == '__main__':
	pass
