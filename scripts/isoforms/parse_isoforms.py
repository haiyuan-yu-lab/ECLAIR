# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script fetches processes all of the UniProt IDs in 
# a given uniprot_info_file to create a feature mask that
# identifies all positions in the UniProt sequence that
# are affected by isoforms / alternate sequences. Wherever
# possible, the script identifies UniProt IDs that have
# previously had this isoform mask written, and only
# writes information for new IDs. The output format
# is two columns, where the first column is the UniProt ID
# and the second column is the isoform mask, where 0
# represents no isoform, 1 represents at least one isoform
# and all positions are joined by semicolons.
#
# e.g. Q9BQI0	0;0;1;1;1;1;1;0;0;0;

# Expected Outcomes:
# All UniProt IDs in the input file should have their
# isoform mask written to the output file.

# Known Bugs:
# - Potential concerns about the caching process and
#   possible loss of information when UniProt entries
#   that we have previously fetched are changed.
# - Isoforms masks are generated for every UniProt ID
#   including the previously fetched. This could be
#   resolved by moving the already_fetched assignment
#   above the main for loop, and performing the iteration
#   on the difference between these two sets.

# Imports
import sys, os
from mjm_tools import parse_dictionary_list

# Set up input file based on input parameter
# If no input is provided, default to fixed
# input file
if len(sys.argv) > 1:
	uniprot_info_file = sys.argv[1]
else:
	uniprot_info_file = '../uniprot/uniprot_info_feat.txt'

# Output to fixed output file
output_file = '/home/adr66/eclair/features/per_feature/isoforms.txt'

# Read UniProt info file as a dictionary, uniprot2info, formatted as...
# Keys: UniProt IDs
# Values: Dictionary representation of the row of attributes
#         in the uniprot_info_file
# e.g: {UniID : {att1 : val1, att2 : val2, ...}}
uniprot2info = dict((e['id'], e) for e in parse_dictionary_list(uniprot_info_file))

# Iterate over each UniProt ID to generate an isoform mask
# where 0s indicate regions of the protein with no isoform
# and 1s indicate regions with at least one isoform
# NOTE: This calculates for ALL UniProts in the info file
#       even the ones that have already been fetched. This
#       code should be modified so that it only generates
#       the mask on new UniProts.
feats = []
for p in uniprot2info:
	# Save bounds for all alternate sequences associated with the protein
	#
	# Alternate sequences are (to my understanding) represented roughly as...
	#
	# TYPE START END ORIG_SEQ -> ALT_SEQ ISOFORM PUBLICATION ID?
	#
	# e.g. 
	# 'VAR_SEQ 397 432 NVQEAQKILNNSGLPITSAIDLEDAAKKAVASVAKK -> FMEKKGSYMHIKQETGNSNENITGIQENVAHQNLSKCAISIFLC
	#  (in isoform 2). {ECO:0000303|PubMed:14702039}. /FTId=VSP_042013.'
	#
	varseqs = []
	for var in uniprot2info[p]['feature(ALTERNATIVE SEQUENCE)'].split('; '):
		# If there are no alternate sequences, do nothing
		if var == '':
			continue
		# Save the start and end position of the alternate sequence
		varseqs.append(var.split(' ')[1:3])
	
	# Create binary mask for all isoforms
	#
	# e.g. if you had a protein of len 10, with an isoform
	#      occurring between positions 2 and 6 and another
	#      between positions 8 and 9 you would have...
	#
	#      ['0', '1', '1', '1', '1', '1', '0', '1', '1', '0']
	#
	result = ['0' for _ in range(int(uniprot2info[p]['length']))]
	for se in varseqs:
		s, e = se
		# Set isoform affected positions to 1
		result[int(s)-1:int(e)] = ['1'] * (int(e) - int(s) + 1)
	
	# Save mask for each protein
	feats.append(result)

# Read in output file to determine which UniProts have
# already been fetched / processed for isoforms.
already_fetched = {}
if os.path.exists(output_file):
	for l in open(output_file):
		p, d = l.split('\t')
		already_fetched[p] = d

# Write to output file by appending all new UniProt
# IDs / isoforms to the file
with open(output_file, 'a') as of:
	for i, p in enumerate(uniprot2info):
		# Skip previously fetched UniProt IDs
		if p not in already_fetched:
			of.write('%s\t%s\n' % (p, ';'.join(feats[i])))
