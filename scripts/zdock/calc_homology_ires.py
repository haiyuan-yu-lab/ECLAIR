
# Batch calculate all IRES for all zdock models

import time, sys, os, glob, re, sys
from mjm_parsers import unzip_res_range, parse_dictionary_list
from mjm_tools import fetch_uniprot, xtool, time_elapsed, zip_res_range
start_time = time.time()

homology_data = parse_dictionary_list(sys.argv[1])
interlist = map(lambda e: [e['UniProt_1'],e['UniProt_2']], homology_data)

docking_results_dir = 'pdb_docked_models'
output_file = 'ires/IRES_homology_docking.txt'
uniprotInfoFile = sys.argv[2]
sifts_file = '../pdb/pdbresiduemapping.txt'

docked_model_files = sorted(glob.glob(os.path.join(docking_results_dir, '*.pdb')))
uniprot2info = dict((e['id'], e) for e in parse_dictionary_list(uniprotInfoFile))
sifts_data = parse_dictionary_list(sifts_file)

#don't recalculate IRES for already processed docking results
already_calculated = {}
if os.path.exists(output_file):
	for l in open(output_file).readlines()[1:]:
		already_calculated[l.split('\t')[5]] = l


pdbChain2siftsEntry = {}
for e in sifts_data:
	pdb_chain = e['PDB']+'_'+e['Chain']
	#if chain has no entry yet, or the currently-parsed entry has less uniprot coverage
	if pdb_chain not in pdbChain2siftsEntry or len(unzip_res_range(e['MappableResInPDBChainOnUniprotBasis'])) > len(unzip_res_range(pdbChain2siftsEntry[pdb_chain]['MappableResInPDBChainOnUniprotBasis'])):
		pdbChain2siftsEntry[pdb_chain] = e

#------------------------------------------------------------------------

print 'Found %i docked models in %s...' %(len(docked_model_files), docking_results_dir)

output = open(output_file, 'w')
output.write('\t'.join(['UniProtA', 'UniProtB', 'Rank', 'SubunitA', 'SubunitB', 'File', 'ZDOCK_Score', 'UniProtA_IRES', 'UniProtB_IRES']) + '\n')

#for fi, f in list(enumerate(docked_model_files)):
#	if fi%100==0: print '%s/%s complete' %(fi, len(docked_model_files))

for e in homology_data:

	fs = glob.glob(os.path.join(docking_results_dir,'%s--%s*.pdb' % (e['PDB_1'],e['PDB_2'])))
	fs += glob.glob(os.path.join(docking_results_dir,'%s--%s*.pdb' % (e['PDB_2'],e['PDB_1'])))

	for f in fs:
		if os.path.basename(f) in already_calculated:
			output.write(already_calculated[os.path.basename(f)])
			continue
	
		subA, subB, _, num, _ = re.split('--|-|\.', os.path.basename(f))
	
		if subA not in pdbChain2siftsEntry or subB not in pdbChain2siftsEntry:
			continue
		
		uniprotA = pdbChain2siftsEntry[subA]['UniProt']
		uniprotB = pdbChain2siftsEntry[subB]['UniProt']
	
	
		#subunit A sifts residue mapping
		subA2uniprot = {}
		pdbres = unzip_res_range(pdbChain2siftsEntry[subA]['MappableResInPDBChainOnPDBBasis'])
		uniprotres = unzip_res_range(pdbChain2siftsEntry[subA]['MappableResInPDBChainOnUniprotBasis'])
	
		for i in range(len(pdbres)):
			subA2uniprot[pdbres[i]] = uniprotres[i]
	
		#subunit B sifts residue mapping
		subB2uniprot = {}
		pdbres = unzip_res_range(pdbChain2siftsEntry[subB]['MappableResInPDBChainOnPDBBasis'])
		uniprotres = unzip_res_range(pdbChain2siftsEntry[subB]['MappableResInPDBChainOnUniprotBasis'])
	
		for i in range(len(pdbres)):
			subB2uniprot[pdbres[i]] = uniprotres[i]
	
		#--------get interface residues----------
		out = [x.strip().split('\t') for x in xtool('irescalc.py', f, c1='A', c2='B', o='residue_stats').strip().split('\n')]
	
		if out == [['']]:
			A_ires, B_ires = 'N/A', 'N/A'
		else:
			A_ires = zip_res_range([subA2uniprot[q[1]] for q in out if q[1] in subA2uniprot and q[0]=='A'])
			B_ires = zip_res_range([subB2uniprot[q[1]] for q in out if q[1] in subB2uniprot and q[0]=='B'])
	
		try:
			zdock_score = open(f[:-7]+'.out').readlines()[int(num)+3].strip().split()[-1]
		except:
			zdock_score = 'n/a'
	
		if uniprotA > uniprotB:
			output.write('\t'.join([uniprotB, uniprotA, num, subB, subA, os.path.basename(f), zdock_score, B_ires, A_ires])+'\n')
		else:
			output.write('\t'.join([uniprotA, uniprotB, num, subA, subB, os.path.basename(f), zdock_score, A_ires, B_ires])+'\n')
	
	
output.close()

print 'Finished calculating IRES for %s docked models: %s' %(len(docked_model_files), time_elapsed(start_time))
