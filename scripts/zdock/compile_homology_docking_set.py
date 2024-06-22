import time, subprocess, os, glob, sys
from collections import defaultdict
from mjm_parsers import parse_dictionary_list, unzip_res_range

homology_file = sys.argv[1]

sifts_file = '../pdb/pdbresiduemapping.txt'
excluded_pdbs_file = '../pdb/excluded_pdbs.txt'
output_file = 'homology_docking_set.txt'

sifts_data = dict()
for e in parse_dictionary_list(sifts_file):
	sifts_data['%s_%s' % (e['PDB'],e['Chain'])] = e
homology_pairs = parse_dictionary_list(homology_file)

#----------------------------------------------------------

if not os.path.exists(output_file):
	with open(output_file,'w') as output:
		output.write('\t'.join(['ProtA', 'ProtB', 'SubA', 'SubB', 'CovA', 'CovB']) + '\n')

output = open(output_file,'a')
for hp in homology_pairs:
	
	p1 = hp['UniProt_1']
	p2 = hp['UniProt_2']
	
	try:
		covered_res1 = set(unzip_res_range(sifts_data[hp['PDB_1']]['MappableResInPDBChainOnUniprotBasis']))
	except:
		covered_res1 = set()

	try:
		covered_res2 = set(unzip_res_range(sifts_data[hp['PDB_2']]['MappableResInPDBChainOnUniprotBasis']))
	except:
		covered_res2 = set()
	coverage1 = len(covered_res1)
	coverage2 = len(covered_res2)
	
	output.write('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (p1,p2,hp['PDB_1'],hp['PDB_2'],coverage1,coverage2))

output.close()
