import numpy as np
import sys

uniprotid = sys.argv[1]

def get_depth(uniprotid):
	with open('/home/adr66/eclair/data/conservation/psiblast/%s_aligned.msa'%uniprotid) as inpfile:
		n = 0 # The length of the target protein sequence
		m = 0 # The total number of PSIBLAST matches
		seq = '' # The target protein sequence

		# FASTA format, one line is a header and the entire sequence is on the following line
		for i, l in enumerate(inpfile):
			if i == 1:
				# l is the sequence of the target protein
				seq = list(l.strip())
				n = len(seq)
				result = np.zeros(n)
			if i > 2 and i % 2 == 1:
				# OLD code: result counts the number of non-gaps at each location
				# result = result + (list(l.strip()) != np.repeat('-',n))

				# NEW code: result counts the number of sequence matches at each residue
				result = result + (np.array(list(l.strip())) == np.array(seq))

				m += 1
	return (result, m)

print get_depth(uniprotid)
