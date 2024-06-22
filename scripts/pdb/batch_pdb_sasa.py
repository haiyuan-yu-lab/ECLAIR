# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script combines the precomputed PDB SASA
# information to generate the SASA PDB feature used
# in Eclair This feature is generated by reporting
# the average SASA per residue over all PDB structures.
# that could be mapped to a UniProt using Sifts. For some
# reason this feature is reported per interaction
# rather than per protein and I am not sure why.

# Expected Outcomes:
# The output file should be updated so that it includes
# entries for the calculation of the average PDB
# SASA at each position for each of the interactions
# included in the input file.

# Known Bugs:
# - Dependence on static copies of the resources folder
# - Interaction order is sorted. I have confirmed there
#   was originally an assumption that the input
#   interactions should be sorted. I suspsect that
#   it is possible some features may be unavailable
#   when the interaction is not properly sorted in
#   the input. Will do testing to confirm.
# - Confused why the output here is per interaction
# - Add warnings

# Imports
import glob, os, sys
from mjm_parsers import parse_dictionary_list
from collections import defaultdict
import numpy as np

# Identify input files form parameters, otherwise
# use defaults
if len(sys.argv) > 1:
	uniprot_info_file = sys.argv[1]
	interactomes_regex = sys.argv[2]
else:
	uniprot_info_file = '../uniprot/uniprot_info.txt'
	interactomes_regex = '../interactomes/*binary_hq.txt'

# Sifts file for mapping UniProt to PDB
# Static copy from resources
sifts_file = 'pdbresiduemapping.txt'

# Precomputed SASA per PDB
# Static copy from resources
sasa_file = 'SASA_perpdb_alltax.txt'

# File for PDBs to exclude?
# This seems to be a mapping of PDBs that should
# be excluded for a given interaction pair
# as well as a boolean indicating whether a
# co-crystal structure for that pair exists?
#
# HEAD:
#
# UniProtA  UniProtB hasCC excludedPDBs
# A0A023PXA5   P47068   N
# A0A023PXP4   P39743   N
# A0A023PYF7   P39743   N
# A0A023PZD0   P53281   N  1KI1;2KGR;2KHN;3FIA;3HS8;3HS9;3JV3;3QBV;4IIM
# A0A023PZD0   P80667   N  1KI1;2KGR;2KHN;3FIA;3HS8;3HS9;3JV3;3QBV;4IIM
#
excluded_pdbs_file = 'excluded_pdbs.txt'

# Fixed output files
output_file = '/home/adr66/eclair/features/per_feature/SASA_PDB.txt'

# Read UniProt Info file as dictionary
uniprot2info = dict((e['id'], e) for e in parse_dictionary_list(uniprot_info_file))

# Obtain a list of all the interactions
interactomes = set()
for f in glob.glob(interactomes_regex):
        print 'parsing', f, '...'
        interactomes.update(  set([tuple(sorted(l.strip().split('\t')[:2])) for l in open(f)])  )


excluded_pdbs = defaultdict(set)
for e in parse_dictionary_list(excluded_pdbs_file):
	
	# If the interaction does not have a crystal structure
	# then we don't have to worry about excluding their
	# structures.
	#
	# This has something to do with "restricting docking
	# subunits" and is linked to things that "would be used
	# to train / test classifier"
	if e['hasCC']=='N':
		continue
	
	# If there are no excluded PDBs for this pair, we can continue
	if e['excludedPDBs'] == '':
		continue
	
	# Otherwise, we sort the interaction and look up the excluded
	# PDBs from our dictionary
	#
	# Note: This is the second place I have seen any mention of
	#       sorting the interaction order. This piece of code
	#       was previously commented to imply that the interaction
	#       should already be sorted. Am afraid there may be an
	#       unstated assumption that the user has already sorted
	#       his / her input. In either case, the final version
	#       of this should be modified so that the interaction
	#       input is automatically sorted as the first step.
	interaction = tuple(sorted(tuple([e['UniProtA'], e['UniProtB']])))
	excluded_pdbs[interaction] = set(e['excludedPDBs'].split(';'))

# Generate dictionary for mapping UniProt --> PDB --> list of Chains + SASA values?
uniprot2chains = defaultdict(lambda: defaultdict(list))

# Iterate over Previously Generated
# PDB SASA File
for e in parse_dictionary_list(sasa_file):
	#Parse SASA list for each PDB structure
	sasas = np.array([ float(r) if r != 'NaN' else np.nan for r in e['UniProt_SASA'].split(';') ])
	
	# Add to dictionary
	# UniProt --> PDB --> SASA
	uniprot2chains[e['UniProt']][e['PDB']].append(sasas)


# Figure out which interactions have aleady had
# their SASA reported in the output
already_fetched = set()
if os.path.exists(output_file):
	for l in open(output_file):
		already_fetched.add(''.join(l.split('\t')[:2]))

# Open output for appending
output = open(output_file, 'a')

# Iterate over each interaction
for i in interactomes:
	# Obtain UniProt IDs
	p1, p2 = i
	
	# Skip entries where either entry is not in
	# the UniProt Info file
	if p1 not in uniprot2info or p2 not in uniprot2info:
		continue
	
	# Skip entries that have already had their
	# SASA reported
	if ''.join((p1,p2)) in already_fetched:
		continue

	# Obtain SASAs for P1
	p1_sasas = []
	for pdb in uniprot2chains[p1]:
		# Skip excluded
		if pdb not in excluded_pdbs[i]:
			p1_sasas += uniprot2chains[p1][pdb]
	
	# Obtain SASAs for P2
	p2_sasas = []
	for pdb in uniprot2chains[p2]:
		# Skip excluded
		if pdb not in excluded_pdbs[i]:
			p2_sasas += uniprot2chains[p2][pdb]
			
	# Skip entries with no remaining SASA values in
	# either protein
	if len(p1_sasas)==0 and len(p2_sasas)==0:
		continue
	
	# If P1 has no values, set the full result to NaN
	if len(p1_sasas) == 0:
		p1_means = [np.nan for r in range(int(uniprot2info[p1]['length']))]
	# Otherwise, average all of the SASAs reported per structure
	else:
		# Create a masked array from all of the reported SASAs
		# The masked array just lets us easilly ignore the NaN
		# values where a given structure does not cover a position
		# when we go to calculate the average per position
		mdat = np.ma.masked_array(p1_sasas, np.isnan(p1_sasas))
		
		# Generate the mean SASA per position over all structures
		p1_means = np.mean(mdat, axis=0)
		
		# Unmask data filling back in NaN for any positions
		# that had no SASA information at all
		p1_means = [p1_means.data[r] if p1_means.mask[r]==False else np.nan for r in range(len(p1_means.data))]
	
	
	# If P2 has no values, set the full result to NaN
	if len(p2_sasas) == 0:
		p2_means = [np.nan for r in range(int(uniprot2info[p2]['length']))]
	# Otherwise, average all of the SASAs reported per structure
	else:
		# Create a masked array from all of the reported SASAs
		# The masked array just lets us easilly ignore the NaN
		# values where a given structure does not cover a position
		# when we go to calculate the average per position
		mdat = np.ma.masked_array(p2_sasas, np.isnan(p2_sasas))
		
		# Generate the mean SASA per position over all structures
		p2_means = np.mean(mdat, axis=0)
		
		# Unmask data filling back in NaN for any positions
		# that had no SASA information at all
		p2_means = [p2_means.data[r] if p2_means.mask[r]==False else np.nan for r in range(len(p2_means.data))]
		
	# Skip entries where the length of either of the UniProts
	# based on SASA features does not match the expected
	# length of the UniProt based on the UniProt Info
	# Should throw a warning
	if len(p1_means) != int(uniprot2info[p1]['length']) or len(p2_means) != int(uniprot2info[p2]['length']):
		print 'PDB SASA: Unexpected lengths for pair %s,%s' % (p1,p2)
		continue
	
	# Write output file
	# Not sure why this has to be a pair-wise
	# feature rather than a protein specific
	# feature?
	#
	# HEAD:
	#
	# P1  P2 P1AvgSASA   P2AvgSASA
	output.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in p1_means]), ';'.join(['%.2f' %(res) for res in p2_means])))

# Close Output
output.close()
