
# Batch calculate all dSASA for all zdock models

import time, sys, os, glob, re
from mjm_parsers import unzip_res_range, parse_dictionary_list
from mjm_tools import fetch_uniprot, xtool, time_elapsed, natural_keys
start_time = time.time()

homology_data = parse_dictionary_list(sys.argv[1])
interlist = map(lambda e: [e['UniProt_1'],e['UniProt_2']], homology_data)

docking_results_dir = 'pdb_docked_models'
output_file = 'ires/dist3d_homology_docking.txt'
uniprotInfoFile = sys.argv[2]
sifts_file = '../pdb/pdbresiduemapping.txt'

docked_model_files = sorted(glob.glob(os.path.join(docking_results_dir, '*.pdb')) + glob.glob(os.path.join(docking_results_dir, '*.pdb.gz')), key=natural_keys)
uniprot2info = dict((e['id'], e) for e in parse_dictionary_list(uniprotInfoFile))
sifts_data = parse_dictionary_list(sifts_file)

#don't recalculate dSASA for already processed docking results
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

if not os.path.exists(output_file):
	output = open(output_file, 'w')
	output.write('\t'.join(['UniProtA', 'UniProtB', 'Rank', 'SubunitA', 'SubunitB', 'File', 'ZDOCK_Score', 'UniProtA_dSASA', 'UniProtB_dSASA']) + '\n')
	output.close()

output = open(output_file,'a')
for e in homology_data:

	fs = glob.glob(os.path.join(docking_results_dir,'%s--%s*.pdb' % (e['PDB_1'],e['PDB_2'])))
	fs+= glob.glob(os.path.join(docking_results_dir,'%s--%s*.pdb' % (e['PDB_2'],e['PDB_1'])))

	for f in fs:
		if os.path.basename(f) in already_calculated:
			continue

		subA, subB, _, num = re.split('--|-|\.', os.path.basename(f))[:4]

		if subA not in pdbChain2siftsEntry or subB not in pdbChain2siftsEntry:
			continue

		uniprotA = pdbChain2siftsEntry[subA]['UniProt']
		uniprotB = pdbChain2siftsEntry[subB]['UniProt']

		if [uniprotA,uniprotB] not in interlist and [uniprotB,uniprotA] not in interlist:
			continue

		if uniprotA not in uniprot2info or uniprotB not in uniprot2info:
			continue

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
		out = [x.strip().split('\t') for x in xtool('chaindistcalc.py', f, c1='A', c2='B', n=10).strip().split('\n')]

		if out == [['']]:
			print 'chaindistcalc failed: %s' %(f)
			continue

		A_dSASA = dict([(int(subA2uniprot[q[1]]), q[3]) for q in out if q[1] in subA2uniprot and q[0]=='A'])  #map UniProt residues to SASA values
		B_dSASA = dict([(int(subB2uniprot[q[1]]), q[3]) for q in out if q[1] in subB2uniprot and q[0]=='B'])  #map UniProt residues to SASA values

		A_dSASA = [str(A_dSASA[r]) if r in A_dSASA else 'nan' for r in range(1, int(uniprot2info[uniprotA]['length'])+1)]
		B_dSASA = [str(B_dSASA[r]) if r in B_dSASA else 'nan' for r in range(1, int(uniprot2info[uniprotB]['length'])+1)]

		try:
			if f[-4:] == '.pdb':
				zdock_score = open(re.sub('ZDOCK.*pdb$', 'ZDOCK.out', f)).readlines()[int(num)+3].strip().split()[-1]
			else:
				zdock_score = open(re.sub('ZDOCK.*pdb.gz$', 'ZDOCK.out', f)).readlines()[int(num)+3].strip().split()[-1]
		except:
			zdock_score = 'n/a'

		if uniprotA > uniprotB:
			output.write('\t'.join([uniprotB, uniprotA, num, subB, subA, os.path.basename(f), zdock_score, ';'.join(B_dSASA), ';'.join(A_dSASA)]) + '\n')
		else:
			output.write('\t'.join([uniprotA, uniprotB, num, subA, subB, os.path.basename(f), zdock_score, ';'.join(A_dSASA), ';'.join(B_dSASA)]) + '\n')

print 'Finished calculating dist3d values for %s docked models: %s' %(len(docked_model_files), time_elapsed(start_time))
