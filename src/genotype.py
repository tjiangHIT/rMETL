# * @author: Jiang Tao (tjiang@hit.edu.cn)


from scipy.special import comb
import numpy as np

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
        # return    0
        return 0
    else:
        return 1

def genotype_call_with_read_pair(concordant, discordant, std_depth=5):
    """
    Genotype call based on total depth of concordant and discordant reads a la
    Hormozdiari et al. 2010 (Genome Research).
    """
    homozygous_deletion_threshold = 5
    if concordant < homozygous_deletion_threshold and discordant < homozygous_deletion_threshold:
        genotype = "./."
        discordant_genotype_likelihood = np.power(2, discordant) / np.power(2, homozygous_deletion_threshold)
        concordant_genotype_likelihood = np.power(2, concordant) / np.power(2, homozygous_deletion_threshold)
        shared_likelihood = np.sqrt(np.square(concordant_genotype_likelihood) + np.square(discordant_genotype_likelihood))
        genotype_likelihood = np.floor(-10 * np.log10(shared_likelihood))
    else:
        expected_discordant_lower_bound = concordant * 0.25
        expected_discordant_upper_bound = concordant * 4
        if discordant < expected_discordant_lower_bound:
            genotype = "1/1"
            genotype_likelihood = np.floor(-10*np.log10(1.0/np.power(2, expected_discordant_lower_bound - discordant)))
        elif expected_discordant_lower_bound <= discordant < expected_discordant_upper_bound:
            if np.abs(discordant - expected_discordant_upper_bound) - expected_discordant_upper_bound >= 0:
                genotype_ratio = np.power(2, discordant)
            else:
                print np.power(2, discordant)
                genotype_ratio = 1.0 / np.power(2, discordant)
            genotype_likelihood = np.floor(-10 * np.log10(1 - min(genotype_ratio, 1)))
            genotype = "1/0"
        else:
            if np.abs(discordant - expected_discordant_upper_bound) - expected_discordant_upper_bound >= 0:
            	genotype_ratio = np.power(2, np.abs(discordant - expected_discordant_upper_bound) - expected_discordant_upper_bound)
            else:
            	genotype_ratio = 1.0/np.power(2, expected_discordant_upper_bound - np.abs(discordant - expected_discordant_upper_bound))
            genotype_likelihood = np.floor(-10 * np.log10(1 - min(genotype_ratio, 1)))
            genotype = "0/0"
    return genotype, genotype_likelihood

def caculate_genotype_likelyhood(Nalt, Ntotal):
	ref_ref = 1.0 / 3
	ref_ME = 1.0 / 3
	ME_ME = 1.0 / 3
	P_ref = ref_ref * comb(Ntotal, Nalt) * P_homozygous_ref ** Nalt * (1 - P_homozygous_ref) ** (Ntotal - Nalt)
	P_het = ref_ME * comb(Ntotal, Nalt) * P_heterozygous ** Nalt * (1 - P_heterozygous) ** (Ntotal - Nalt)
	P_hom = ME_ME * comb(Ntotal, Nalt) * P_homozygous_ME ** Nalt * (1 - P_homozygous_ME) ** (Ntotal - Nalt)
	P_total = P_ref + P_het + P_hom
	P_ref = P_ref / P_total
	P_het = P_het / P_total
	P_hom = P_hom / P_total
	genotype_likelyhood = [P_ref, P_het, P_hom]
	final_gl = max(genotype_likelyhood)
	tag = genotype_likelyhood.index(final_gl) 
	return GL_TAG[tag], final_gl

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
	genotype_call_with_read_pair(29, 5)
