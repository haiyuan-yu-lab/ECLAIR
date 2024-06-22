# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script calculates the average atomic distances
# between the alpha carbon (CA) atoms for each residue
# in protein 1 to any of the CA atoms in protein 2 for
# all of the docking results in the PDB docking
# set. dist3d values are calculated using mjms
# chaindistcalc.py. Wherever possible the script avoids
# recalculating dist3d values for docking results
# that have already had their dist3d values
# calculated.

# Expected Outcomes:
# The output file should be updated to include entries
# for dist3d values of all residues for all of the docking
# results that were included in the PDB docking
# set.

# Known Bugs:

# Imports
import time, sys, os, glob, re
from mjm_parsers import unzip_res_range, parse_dictionary_list
from mjm_tools import fetch_uniprot, xtool, time_elapsed, natural_keys

# Keep track of runtime
start_time = time.time()

# Obtain list of all interactions in input file
interlist = [sorted(l.strip().split('\t')[:2]) for l in open(sys.argv[1])]

# Path to docking results Directory
docking_results_dir = 'pdb_docked_models'

# Fixed path to output file
output_file = 'ires/dist3d_pdb_docking.txt'

# UniProt Info File from Input
uniprotInfoFile = sys.argv[2]

# SIFTs Mapping file
# This is a static copy of our resources
sifts_file = '../pdb/pdbresiduemapping.txt'

# Obtain list of all docking results
docked_model_files = sorted(glob.glob(os.path.join(docking_results_dir, '*.pdb')) + glob.glob(os.path.join(docking_results_dir, '*.pdb.gz')), key=natural_keys)

# Read UniProt Info file as dictionary
uniprot2info = dict((e['id'], e) for e in parse_dictionary_list(uniprotInfoFile))

# Read SIFTs Mapping file as dictionary
sifts_data = parse_dictionary_list(sifts_file)

# Obtain list of all docking results that have
# already had their delta SASAs calculated
# Here we save the actual line (but we don't
# actually use it?)
already_calculated = {}
if os.path.exists(output_file):
	with open(output_file) as f:
		f.read()
		for l in f:
			already_calculated[l.split('\t')[5]] = l

# Obtain mapping of PDB ID + Chain to UniProt
# Only map each chain to the UniProt ID with the
# highest coverage?
# Something seems off here?
pdbChain2siftsEntry = {}
for e in sifts_data:
	pdb_chain = e['PDB']+'_'+e['Chain']
	
	# Save each PDB_Chain pair to map to the UniProt
	# ID that has the highest coverage (most mappable residues)
	if pdb_chain not in pdbChain2siftsEntry or len(unzip_res_range(e['MappableResInPDBChainOnUniprotBasis'])) > len(unzip_res_range(pdbChain2siftsEntry[pdb_chain]['MappableResInPDBChainOnUniprotBasis'])):
		pdbChain2siftsEntry[pdb_chain] = e

print 'Found %i docked models in %s...' %(len(docked_model_files), docking_results_dir)

# If the output files does not exist,
# open it up for writting and add the
# header lines
if not os.path.exists(output_file):
	output = open(output_file, 'w')
	output.write('\t'.join(['UniProtA', 'UniProtB', 'Rank', 'SubunitA', 'SubunitB', 'File', 'ZDOCK_Score', 'UniProtA_dSASA', 'UniProtB_dSASA']) + '\n')
	output.close()

# Break the docking results into batches of 1000
# I'm not sure why this is suddently done using batches
# rather than just one for loop?
# I think this is optimizing for a potentially non-existent
# trade off between write operation time vs. cost of storing
# all of the lines in memory at once.
# Technically the fastest should just be to write ALL of the
# output in one write operation (I don't think the memory
# cost should be too large). Need to do testing / standardize
# the decision in all of the scripts.
batch_size = 1000
batches = [docked_model_files[i*batch_size:(i+1)*batch_size] for i in range(0, len(docked_model_files)/batch_size+1)]

# Iterate over each batch
for batch_num, docked_batch in enumerate(batches):
	# Status Update
	print '%d/%d complete' % (batch_num*batch_size, len(docked_model_files))
	
	# Store the output lines for this batch
	batch_output = []
	
	# Iterate over each docking result
	for num, f in enumerate(docked_batch):
		# Skip entries that have already had their delta
		# SASA values calculated
		if os.path.basename(f) in already_calculated:
			continue
		
		# Obtain the two Models that made up the dock and
		# the rank of the docking result
		subA, subB, _, num = re.split('--|-|\.', os.path.basename(f))[:4]
		
		# Skip entries where either structure is not included
		# in the SIFTs mapping
		if subA not in pdbChain2siftsEntry or subB not in pdbChain2siftsEntry:
			continue
		
		# Obtain the UniProt IDs corresponding to these
		# structures
		uniprotA = pdbChain2siftsEntry[subA]['UniProt']
		uniprotB = pdbChain2siftsEntry[subB]['UniProt']
		
		# Skip entries where either UniProt ID is not represented
		# in the UniProt Info File
		if uniprotA not in uniprot2info or uniprotB not in uniprot2info:
			continue
		
		# Skip entries where the docking results do not
		# represent one of the interactions included in
		# the current input file
		if [uniprotA,uniprotB] not in interlist and [uniprotB,uniprotA] not in interlist:
			continue
		
		# Generate dicitonary mapping PDB positions to
		# UniProt positions
		subA2uniprot = {}
		pdbres = unzip_res_range(pdbChain2siftsEntry[subA]['MappableResInPDBChainOnPDBBasis'])
		uniprotres = unzip_res_range(pdbChain2siftsEntry[subA]['MappableResInPDBChainOnUniprotBasis'])
		
		# Update the dictionary
		for i in range(len(pdbres)):
			subA2uniprot[pdbres[i]] = uniprotres[i]
		
		# Generate dicitonary mapping PDB positions to
		# UniProt positions
		subB2uniprot = {}
		pdbres = unzip_res_range(pdbChain2siftsEntry[subB]['MappableResInPDBChainOnPDBBasis'])
		uniprotres = unzip_res_range(pdbChain2siftsEntry[subB]['MappableResInPDBChainOnUniprotBasis'])
		
		# Update the dictionary
		for i in range(len(pdbres)):
			subB2uniprot[pdbres[i]] = uniprotres[i]
		
		# Calculate average distance to the other chain
		# by calling mjm chaindistcalc.py. The output looks like...
		#
		# HEAD:
		#
		# Chain	Pos	Resn	avgDist
		# A	7	PHE	63.6398
		# A	8	GLU	62.4268
		# A	9	THR	65.2177
		# A	10	ASP	63.1981
		# A	11	MET	60.1369
		#
		out = [x.strip().split('\t') for x in xtool('chaindistcalc.py', f, c1='A', c2='B', n=10).strip().split('\n')]
		
		# If there is no output, then there the calculation
		# is considered to have failed
		if out == [['']]:
			print 'chaindistcalc failed: %s' %(f)
			continue
		
		# Split the dSASA values by chain / generate
		# dictionary (pos --> dSASA) (After mapping PDB Positions
		# back to UniProt positions
		A_dSASA = dict([(int(subA2uniprot[q[1]]), q[3]) for q in out if q[1] in subA2uniprot and q[0]=='A'])  #map UniProt residues to SASA values
		B_dSASA = dict([(int(subB2uniprot[q[1]]), q[3]) for q in out if q[1] in subB2uniprot and q[0]=='B'])  #map UniProt residues to SASA values
		
		# Update lists to include NaN values for all of the positions
		# that were not covered by the PDB structure
		A_dSASA = [str(A_dSASA[r]) if r in A_dSASA else 'nan' for r in range(1, int(uniprot2info[uniprotA]['length'])+1)]
		B_dSASA = [str(B_dSASA[r]) if r in B_dSASA else 'nan' for r in range(1, int(uniprot2info[uniprotB]['length'])+1)]
		
		# Attempt to parse the Zdock output file to obtain the
		# score for this docking result. I am not entirely sure
		# how to interpret these outputs, but I have made a guess.
		#
		# 1) The first four lines seem to be setup / input parameters
		#    - For instance, we can tell the third and fourth line
		#      are clearly specific to the two input models that were
		#      used in the docking (receptor and ligand) respectively
		# 2) The remaining lines seem to be diagnostic output for each
		#    of the 2000 docks that there generated, ranked in order of
		#    the docking score (which is the last field in these columns).
		#
		# In this script we just want to grab the score (last column) for
		# the Nth ranked dock, which can be found on row N + 3
		#
		# HEAD:
		#
		# 140   1.2
		# 0.000000 0.000000 0.000000
		# 0006F41AE20CAED003C84CF1841CE2B7.R  17.880   -3.594   20.657
		# 35BFCA6537FC02352701FDE04AB07347.L  16.047   8.447 4.758
		# -1.308997   0.419640 2.974122 38 5  111   4.808
		# -2.094395   0.419640 2.974122 42 128   112   4.647
		# -3.141593   1.998705 2.370472 118   119   36 4.090
		# -1.308997   0.297065 2.160177 25 126   99 3.738
		#
		try:
			if f[-4:] == '.pdb':
				zdock_score = open(re.sub('ZDOCK.*pdb$', 'ZDOCK.out', f)).readlines()[int(num)+3].strip().split()[-1]
			else:
				zdock_score = open(re.sub('ZDOCK.*pdb.gz$', 'ZDOCK.out', f)).readlines()[int(num)+3].strip().split()[-1]
		# If we can't read the output file or the ouput file is
		# missing, just set the score to n/a
		except:
			zdock_score = 'n/a'
		
		# Add lines to batch_output. If the UniProt IDs are not
		# already sorted, flip all of the output columns so
		# that the output corresponds to the sorted interaction
		# pair. The output format is as follows...
		#
		# HEAD:
		#
		# UniProtA UniProtB Rank  SubunitA SubunitB File  ZDOCK_Score UniProtA_IRES  UniProtB_IRES
		# O00757   P01111   01 0006F41AE20CAED003C84CF1841CE2B7 35BFCA6537FC02352701FDE04AB07347 0006F41AE20CAED003C84CF1841CE2B7--35BFCA6537FC02352701FDE04AB07347--ZDOCK-01.pdb 4.808 [233] [123]
		# O00757   P01111   02 0006F41AE20CAED003C84CF1841CE2B7 35BFCA6537FC02352701FDE04AB07347 0006F41AE20CAED003C84CF1841CE2B7--35BFCA6537FC02352701FDE04AB07347--ZDOCK-02.pdb 4.647 [233] [123]
		# O00757   P01111   03 0006F41AE20CAED003C84CF1841CE2B7 35BFCA6537FC02352701FDE04AB07347 0006F41AE20CAED003C84CF1841CE2B7--35BFCA6537FC02352701FDE04AB07347--ZDOCK-03.pdb 4.090 [23]  [123]
		# O00757   P01111   04 0006F41AE20CAED003C84CF1841CE2B7 35BFCA6537FC02352701FDE04AB07347 0006F41AE20CAED003C84CF1841CE2B7--35BFCA6537FC02352701FDE04AB07347--ZDOCK-04.pdb 3.738 [230,235]   [128]
		#
		if uniprotA > uniprotB:
			batch_output.append('\t'.join([uniprotB, uniprotA, num, subB, subA, os.path.basename(f), zdock_score, ';'.join(B_dSASA), ';'.join(A_dSASA)]))
		else:
			batch_output.append('\t'.join([uniprotA, uniprotB, num, subA, subB, os.path.basename(f), zdock_score, ';'.join(A_dSASA), ';'.join(B_dSASA)]))
	
	# Actually write the batch_output
	with open(output_file,'a') as output:
		output.write('\n'.join(batch_output) + '\n')

# Close the Output
output.close()

# Final Status Update / Total runtime
print 'Finished calculating dist3d values for %s docked models: %s' %(len(docked_model_files), time_elapsed(start_time))

# Implementation for single for loop, many
# write operation output. I'm just leaving this
# here for posterity
'''
for fi, f in list(enumerate(docked_model_files)):
	if fi%100==0: print '%s/%s complete' %(fi, len(docked_model_files))

	#rewrite already calculated entries:
	if os.path.basename(f) in already_calculated:
		output.write(already_calculated[os.path.basename(f)])
		continue
	
	subA, subB, _, num = re.split('--|-|\.', os.path.basename(f))[:4]
	
	if subA not in pdbChain2siftsEntry or subB not in pdbChain2siftsEntry:
		continue
		
	#~ if subA == subB:  #parse only heterodimers
		#~ continue
	
	uniprotA = pdbChain2siftsEntry[subA]['UniProt']
	uniprotB = pdbChain2siftsEntry[subB]['UniProt']
	
	if [uniprotA,uniprotB] not in interlist and [uniprotB,uniprotA] not in interlist:
		continue

	if uniprotA not in uniprot2info or uniprotB not in uniprot2info:
		continue
	
	#~ try:
	
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
		output.write('\t'.join([uniprotB, uniprotA, num, subB, subA, os.path.basename(f), zdock_score, ';'.join(B_dSASA), ';'.join(A_dSASA)])+'\n')
	else:
		output.write('\t'.join([uniprotA, uniprotB, num, subA, subB, os.path.basename(f), zdock_score, ';'.join(A_dSASA), ';'.join(B_dSASA)])+'\n')
	
	
	#~ except:
		#~ continue
		
output.close()
'''