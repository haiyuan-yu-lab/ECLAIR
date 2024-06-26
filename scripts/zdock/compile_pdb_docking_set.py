# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script compiles a docking set to be used
# for performing docking simulation using PDB
# structures for each interaction provided in the input
# files. The docking set is generated by parsing
# the SIFTs mapped structures and selecting for each
# interaction, the single pair of PDB structures
# that provides the greatest coverage for each
# individual UniProt protein in the interaction.

# Expected Outcomes:
# All interactions in the input file that have
# valid SIFTs mapped structures should be represented
# in the output PDB docking set.

# Known Bugs:
# - Dependence on static resource copies.
# - Unsure whether using a single docking pair
#   for every interaction is the best practice
#   or if it would be possible to encorporate
#   multiple models that cover unique regions
#   of the protein / combining the results.
# - Questions about restricting structures to
#   exclude large complexes.
# - Should add warnings.

# Imports
import time, subprocess, os, glob, sys
from collections import defaultdict
from mjm_parsers import parse_dictionary_list, unzip_res_range

# Identify input UniProt info / Interactions
# files from input parameters, otherwise
# use defaults
if len(sys.argv) > 1:
	uniprot_info_file = sys.argv[1]
	interactomes_regex = sys.argv[2]
else:
	interactomes_regex = '../interactomes/*_binary_hq.txt'
	uniprot_info_file = '../uniprot/uniprot_info.txt'


# Sifts Mapping Files
# Static copy from resources
sifts_file = '../pdb/pdbresiduemapping.txt'

# List of excluded PDB structures
excluded_pdbs_file = '../pdb/excluded_pdbs.txt'

# Fixed Output file
output_file = 'pdb_docking_set.txt'

# Obtain list of all the interactions
# in the input file
interactomes = set()
for f in glob.glob(interactomes_regex):
	print 'parsing', f, '...'
	interactomes.update(  set([tuple(sorted(l.strip().split('\t')[:2])) for l in open(f)][1:])  )

# Obtain list of PDB structures
# to be excluded
excluded_pdbs = defaultdict(set)
for e in parse_dictionary_list(excluded_pdbs_file):
	
	# Refer to previous scripts for comments
	if e['hasCC']=='N':
		continue
	if e['excludedPDBs'] == '':
		continue
	
	interaction = tuple(sorted(tuple([e['UniProtA'], e['UniProtB']])))  #should already be sorted, but doesn't hurt to do it again
	excluded_pdbs[interaction] = set(e['excludedPDBs'].split(';'))

# Read UniProt Info file as dictionary
uniprot2info = dict((e['id'], e) for e in parse_dictionary_list(uniprot_info_file))

# Read SIFTS Mapping as dicitonary
sifts_data = parse_dictionary_list(sifts_file)

# Iterate over every Sifts mapping to generate
# a dictionary mapping UniProt ID to list
# of valid PDB structures to use for docking
uniprot2pdb = defaultdict(list)
for e in sifts_data:	
	# Skip entries that do not map to any of the
	# UniProts in the UniProt Info file
	# Warning?
	if e['UniProt'] not in uniprot2info:
		continue
	
	# Skip entries that contain two-letter chains
	# These PDB structures are generally very large
	# complexes (>26 subunits) and are not good
	# for docking
	if len(e['Chain']) > 1:
		continue
	
	# Obtain set of residues covered by this model
	covered_res = set(unzip_res_range(e['MappableResInPDBChainOnUniprotBasis']))
	coverage = len(covered_res)
	
	# Skip entries with insufficiently low coverage
	# (raw counts)
	if coverage < 50:
		continue
	
	# Skip entries with insufficiently low coverage
	# (as a proportion of total UniProt length)
	if coverage / float(uniprot2info[e['UniProt']]['length']) < 0.50:
		continue
	
	# Add structure to list for docking
	uniprot2pdb[e['UniProt']].append((coverage, covered_res, e['PDB'], e['Chain']))


# Prioritize Structure option by coverage
for u in uniprot2pdb:
	uniprot2pdb[u].sort(reverse=True)

# Figure out which interactions were
# already in the previous docking set
already_fetched = {}
if os.path.exists(output_file):
	with open(output_file) as outf:
		_ = outf.readline()
		for l in outf:
			already_fetched[','.join(l.split('\t')[:2])] = l

# Open output for writting / re-write all of
# the previously fetched docking sets
# Why not just use append?
output = open(output_file, 'w')
output.write('\t'.join(['ProtA', 'ProtB', 'SubA', 'SubB', 'CovA', 'CovB', 'pCovA', 'pCovB']) + '\n')
for interaction in already_fetched:
	output.write(already_fetched[interaction])

# Iterate over each interaction to find the
# best structure / chain to use for each protein
for interaction in interactomes:
	# Obtain UniProt IDs
	p1, p2 = interaction
	
	# Skip entries that were already fetched
	if ','.join(interaction) in already_fetched:
		continue
	
	# Restrict valid ModBase models to use for each
	# protein by removing any models that use excluded
	# PDB structures as their templates
	p1_subunits = [pdb for pdb in uniprot2pdb[p1] if pdb[2] not in excluded_pdbs[interaction]]
	p2_subunits = [pdb for pdb in uniprot2pdb[p2] if pdb[2] not in excluded_pdbs[interaction]]
	
	# Skip entries that have no valid models left
	if len(p1_subunits) < 1 or len(p2_subunits) < 1:
		continue
	
	# Use the highest coverage model for each protein
	subA = p1_subunits[0]
	subB = p2_subunits[0]
	
	# Write the optimal (highest coverage) docking
	# set fo the output file
	#
	# HEAD:
	#
	# ProtA ProtB SubA  SubB  CovA  CovB  pCovA pCovB
	# P10144   P34932   db937cd1d7682ac28e3ff2cacdf9c6ae dbeb6e05488502eb0776e94cd1a07745 238   697   0.9636   0.8298
	# O60678   P15880   3a8423ea01114f0cfdc0098dcca59f65 4e451b7f4492d36062a2bc5a98d1b2d7 312   226   0.5876   0.7713
	# P45880   Q53ES0   df3aebbc41c449627aedb7c081e35ea2 6a64878d8b92022fc7ec7d621bc9b728 285   998   0.9694   0.9852
	# P32969   V9HW01   cf79feb84ec1c42eff84dd266acccad2 f3fc1cc36bfdf9654ad339fbce35b011 189   133   0.9844   0.8471
	# Q13404   Q96FA3   b29d1054065d33b51fe1e48502be682e 8c181e092f5ed4313aba901162b8cb39 147   266   1.0000   0.6364
	#
	output.write('%s\t%s\t%s\t%s\t%s\t%s\t%.4f\t%.4f\n' %(p1, p2, subA[2]+'_'+subA[3], subB[2]+'_'+subB[3], subA[0], subB[0], subA[0]/float(uniprot2info[p1]['length']), subB[0]/float(uniprot2info[p2]['length'])))

# Close Output
output.close()
