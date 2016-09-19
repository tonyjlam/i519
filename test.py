from itertools import product

def permutation_gen(kmer):
	l = list(product('ATGC', repeat = kmer))
	perm = []
	for i in l:
		perm.append(''.join(i))
	return perm	

def count_kmers():
	kmer_count = dict()
	
	for i in permutation_gen(2):
		kmer_count[i] = 1
	
	print kmer_count

count_kmers()
