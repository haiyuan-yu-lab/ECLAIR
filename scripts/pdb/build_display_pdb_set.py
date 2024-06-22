from mjm_parsers import parse_dictionary_list
from mjm_tools import unzip_res_range
from collections import defaultdict

input_file = 'ires_perpdb_alltax.txt'
output_file = 'pdb_display_structures.txt'
max_pdbs = 3

interaction2entries = defaultdict(list)

for e in parse_dictionary_list(input_file):
	interaction = e['UniProtA'], e['UniProtB']
	all_ires = [r+'A' for r in unzip_res_range(e['UniProtIresA'])] + [r+'B' for r in unzip_res_range(e['UniProtIresB'])]
	if len(all_ires) < 2: continue
	
	interaction2entries[interaction].append(e)
	interaction2entries[interaction][-1].update({'all_ires': set(all_ires[:])})


all_keep_pdbs = []

for i in interaction2entries:
	pdbs = interaction2entries[i]
	if len(pdbs)==0: continue
	
	pdbs.sort(reverse=True, key=lambda d: len(d['all_ires']))	
	
	keep_pdbs = [pdbs[0],]
	seen_ires = set(pdbs[0]['all_ires'])
	
	pdbs = pdbs[1:]
	
	while len(pdbs) >= 1 and len(keep_pdbs) < max_pdbs:
		ires_diffs = []
		for p in pdbs:
			new_ires = p['all_ires'] - seen_ires
			ires_diffs.append((len(new_ires), p))
		ires_diffs.sort(reverse=True)
		
		if ires_diffs[0][0]==0:
			break
		
		keep_pdbs.append(ires_diffs[0][1])
		seen_ires.update(keep_pdbs[-1]['all_ires'])
		pdbs = [i[1] for i in ires_diffs[1:]]
	
	all_keep_pdbs += keep_pdbs


header = ['UniProtA', 'UniProtB', 'PDB', 'ChainA', 'ChainB', 'TaxIDA', 'TaxIDB', 'NumIresA', 'NumIresB', 'UniProtIresA', 'UniProtIresB', 'PDBIresA', 'PDBIresB', 'UniProtAllResA', 'UniProtAllResB', 'PDBAllResA', 'PDBAllResB']
output = open(output_file, 'w')
output.write('\t'.join(header)+'\n')
for e in sorted(all_keep_pdbs, key=lambda d: (d['UniProtA'], d['UniProtB'])):
	output.write('\t'.join([e[h] for h in header]) + '\n')
output.close()
		
		







