#!/usr/bin/env python

from mjm_tools import get_aa_single
from mjm_parsers import parse_fasta
import os, sys

import numpy as np
from scipy.stats import pearsonr

def nanpearsonr(x, y, cutoff=2):
    '''cutoff: number of remaining elements after nan remove to return non-nan value'''
    finite_indices = np.where(np.isfinite(x) & np.isfinite(y))[0]
    if len(finite_indices) < max(cutoff, 2):
        return (np.nan, np.nan)
    return pearsonr(x[finite_indices], y[finite_indices])

'''
http://web.expasy.org/protscale/
1. Hphob. / Kyte & Doolittle
2. Polarity / Grantham
3. Average area buried
4. % accessible residues
5. Transmembrane tendency
6. Bulkiness
7. Amino acid composition
'''

expasy = {

'HPHO': {
'Ala':  1.800, 'Arg': -4.500, 'Asn': -3.500, 'Asp': -3.500, 'Cys':  2.500, 'Gln': -3.500, 'Glu': -3.500, 'Gly': -0.400, 'His': -3.200, 'Ile':  4.500, 'Leu':  3.800, 'Lys': -3.900, 'Met':  1.900, 'Phe':  2.800, 'Pro': -1.600, 'Ser': -0.800, 'Thr': -0.700, 'Trp': -0.900, 'Tyr': -1.300, 'Val':  4.200,
},

'POLA': {
'Ala':  8.100, 'Arg': 10.500, 'Asn': 11.600, 'Asp': 13.000, 'Cys':  5.500, 'Gln': 10.500, 'Glu': 12.300, 'Gly':  9.000, 'His': 10.400, 'Ile':  5.200, 'Leu':  4.900, 'Lys': 11.300, 'Met':  5.700, 'Phe':  5.200, 'Pro':  8.000, 'Ser':  9.200, 'Thr':  8.600, 'Trp':  5.400, 'Tyr':  6.200, 'Val':  5.900,
},

'AREA': {
'Ala': 86.600, 'Arg': 162.200, 'Asn': 103.300, 'Asp': 97.800, 'Cys': 132.300, 'Gln': 119.200, 'Glu': 113.900, 'Gly': 62.900, 'His': 155.800, 'Ile': 158.000, 'Leu': 164.100, 'Lys': 115.500, 'Met': 172.900, 'Phe': 194.100, 'Pro': 92.900, 'Ser': 85.600, 'Thr': 106.500, 'Trp': 224.600, 'Tyr': 177.700, 'Val': 141.000,
},

'ACCE': {
'Ala':  6.600, 'Arg':  4.500, 'Asn':  6.700, 'Asp':  7.700, 'Cys':  0.900, 'Gln':  5.200, 'Glu':  5.700, 'Gly':  6.700, 'His':  2.500, 'Ile':  2.800, 'Leu':  4.800, 'Lys': 10.300, 'Met':  1.000, 'Phe':  2.400, 'Pro':  4.800, 'Ser':  9.400, 'Thr':  7.000, 'Trp':  1.400, 'Tyr':  5.100, 'Val':  4.500, 
},

'TRAN': {
'Ala':  0.380, 'Arg': -2.570, 'Asn': -1.620, 'Asp': -3.270, 'Cys': -0.300, 'Gln': -1.840, 'Glu': -2.900, 'Gly': -0.190, 'His': -1.440, 'Ile':  1.970, 'Leu':  1.820, 'Lys': -3.460, 'Met':  1.400, 'Phe':  1.980, 'Pro': -1.440, 'Ser': -0.530, 'Thr': -0.320, 'Trp':  1.530, 'Tyr':  0.490, 'Val':  1.460,
},

'BULK': {
'Ala': 11.500, 'Arg': 14.280, 'Asn': 12.820, 'Asp': 11.680, 'Cys': 13.460, 'Gln': 14.450, 'Glu': 13.570, 'Gly': 3.400, 'His': 13.690, 'Ile': 21.400, 'Leu': 21.400, 'Lys': 15.710, 'Met': 16.250, 'Phe': 19.800, 'Pro': 17.430, 'Ser': 9.470, 'Thr': 15.770, 'Trp': 21.670, 'Tyr': 18.030, 'Val': 21.570
},

'COMP': {
'Ala':8.25, 'Arg':5.53, 'Asn':4.06, 'Asp':5.45, 'Cys':1.37, 'Gln':3.93, 'Glu':6.75, 'Gly':7.07, 'His':2.27, 'Ile':5.96, 'Leu':9.66, 'Lys':5.84, 'Met':2.42, 'Phe':3.86, 'Pro':4.70, 'Ser':6.56, 'Thr':5.34, 'Trp':1.08, 'Tyr':2.92, 'Val':6.87
}

}


#add single-letter AA keys in each expasy dict:
for prop in expasy.keys():
	for aa in expasy[prop].keys():
		expasy[prop][get_aa_single(aa)] = expasy[prop][aa]
		
		for code in ['-', 'X', 'B', 'Z', 'J', 'O', 'U', '?', '.', ':', '|', '*', ' ']:
			expasy[prop][code] = np.nan



#------------------------------------------------------------------------------------

input_msa = sys.argv[1]
output_file = sys.argv[2]
prop = sys.argv[3]  #i.e. POLA, HPHO, etc.

sequences = parse_fasta(input_msa)[1].values()

sequence_arrays = []
for seq in sequences:
	sequence_arrays.append(np.array([expasy[prop][aa] for aa in seq]))

sequence_matrix = np.vstack(sequence_arrays)

output = open(output_file, 'w')
for i in range(sequence_matrix.shape[1]):
	for j in range(i+1, sequence_matrix.shape[1]):
		
		corr = nanpearsonr(sequence_matrix[:,i], sequence_matrix[:,j], cutoff=2)[0]
		
		output.write('%i\t%i\t%.6f\n' %(i+1, j+1, corr))

output.close()
