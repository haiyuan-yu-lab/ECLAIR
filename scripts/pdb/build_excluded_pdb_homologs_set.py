import os, glob
from collections import defaultdict
from mjm_parsers import parse_dictionary_list

interactomes_dir = '../interactomes'
sifts_file = '../pdb/pdbresiduemapping.txt'
homologs_file = '../uniprot/homologs/pdb/all_homologs.txt'
output_file = 'excluded_pdbs.txt'

homologs = defaultdict(set)
for l in open(homologs_file):
	u = l.split('\t')[0]
	h = set(l.strip().split('\t')[1].split(';'))
	homologs[u].update(h)


interactomes = set()
interactome_files = glob.glob(os.path.join(interactomes_dir, '*_binary_hq.txt'))
for f in interactome_files:
	print 'parsing', f, '...'
	interactomes.update(  set([tuple(sorted(l.strip().split('\t')[:2])) for l in open(f)][1:])  )


#ensure homologs file includes proteins as their own homolog (important later)
for i in interactomes:
	p1, p2 = i
	homologs[p1].add(p1)
	homologs[p2].add(p2)

	
pdb2uniprots = defaultdict(set)  #store all uniprots seen in each PDB
pdbuniprot2count = defaultdict(int)  #store number of times a uniprot is seen in each PDB
uniprot2pdb = defaultdict(set)  #all pdbs associted with uniprot and its homologs (reduce the set of uniprots to check for each interaction)
for e in parse_dictionary_list(sifts_file):
	pdb2uniprots[e['PDB']].add(e['UniProt'])
	pdbuniprot2count[(e['PDB'], e['UniProt'])] += 1
	
	homologs[e['UniProt']].add(e['UniProt'])
	
	for u in homologs[e['UniProt']]:
		uniprot2pdb[u].add(e['PDB'])


output = open(output_file, 'w')
output.write('\t'.join(['UniProtA', 'UniProtB', 'hasCC', 'excludedPDBs'])+'\n')

print 'Calculating excluded PDBs for %i interactions in %i interactomes...' %(len(interactomes), len(interactome_files))
for i_i, i in enumerate(sorted(interactomes)):
	if i_i%1000==0: print i_i,
	
	p1, p2 = i
	excluded_pdbs = set()	
	has_CC = 'N'
	
	for pdb in uniprot2pdb[p1].union(uniprot2pdb[p2]):
		
		#homodimers
		if p1==p2:
			if pdbuniprot2count[(pdb, p1)] > 1:
				excluded_pdbs.add(pdb)
				has_CC = 'Y'
			
			num_homologs_in_pdb = sum([pdbuniprot2count[(pdb, h)] for h in homologs[p1]])
			if num_homologs_in_pdb > 1:
				excluded_pdbs.add(pdb)
		
		#heterodimers
		else:
			if p1 in pdb2uniprots[pdb] and p2 in pdb2uniprots[pdb]:
				excluded_pdbs.add(pdb)
				has_CC = 'Y'
			
			if len(homologs[p1].intersection(pdb2uniprots[pdb])) > 0 and len(homologs[p2].intersection(pdb2uniprots[pdb])) > 0:
				excluded_pdbs.add(pdb)
			
	
	output.write('%s\t%s\t%s\t%s\n' %(p1, p2, has_CC, ';'.join(sorted(excluded_pdbs))))

output.close()

			
