
from scipy.special import comb


P_homozygous_ref = 0.02
P_heterozygous = 0.5
P_homozygous_ME = 0.98

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
	print P_ref, P_het, P_hom

if __name__ == '__main__':
	caculate_genotype_likelyhood(10, 40)