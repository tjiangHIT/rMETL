
from scipy.special import comb
# from math import log10
import numpy as np


P_homozygous_ref = 0.02
P_heterozygous = 0.5
P_homozygous_ME = 0.98

GL_TAG = ['0/0', '0/1', '1/1']

def caculate_genotype_likelyhood(Nalt, Ntotal):
	ref_ref = 1.0 / 3
	ref_ME = 1.0 / 3
	ME_ME = 1.0 / 3

	P_ref = ref_ref * comb(Ntotal, Nalt) * P_homozygous_ref ** Nalt * (1 - P_homozygous_ref) ** (Ntotal - Nalt)
	P_het = ref_ME * comb(Ntotal, Nalt) * P_heterozygous ** Nalt * (1 - P_heterozygous) ** (Ntotal - Nalt)
	P_hom = ME_ME * comb(Ntotal, Nalt) * P_homozygous_ME ** Nalt * (1 - P_homozygous_ME) ** (Ntotal - Nalt)
	P_total = P_ref + P_het + P_hom
	# P_ref = -10 * log10(P_ref / P_total)
	# P_het = -10 * log10(P_het / P_total)
	# P_hom = -10 * log10(P_hom / P_total)
	P_ref = P_ref / P_total
	P_het = P_het / P_total
	P_hom = P_hom / P_total
	genotype_likelyhood = [P_ref, P_het, P_hom]
	final_gl = max(genotype_likelyhood)
	tag = genotype_likelyhood.index(final_gl) 
	# P_ref = np.floor(-10 * np.log10(P_ref / P_total))
	# P_het = np.floor(-10 * np.log10(P_het / P_total))
	# P_hom = np.floor(-10 * np.log10(P_hom / P_total))
	# np.floor(-10 * np.log10(np.power(2, discordant) / np.power(2, expected_discordant_lower_bound)))
	return GL_TAG[tag], final_gl

# if __name__ == '__main__':
# 	caculate_genotype_likelyhood(20, 24)