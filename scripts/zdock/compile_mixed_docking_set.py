# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script compiles a docking set to be used
# for performing docking simulation using a combination
# of one ModBase model and one PDB structure for each 
# interaction provided in the input files. The docking
# set is generated by parsing the valid combinations
# of ModBase and PDB sources and selecting for each
# interaction, the single pair of structure
# that provides the greatest coverage for each
# individual UniProt protein in the interaction.

# Expected Outcomes:
# All interactions in the input file that have
# valid ModBase models and PDB structure should
# be represented in the output ModBase docking set.

# Known Bugs:
# - It looks like the select_modbase_models file
#   that it uses is not actually the file that is
#   actively being updated by Eclair but rather
#   an old one in mjm659.
# - Dependence on static resource copies.
# - Unsure whether using a single docking pair
#   for every interaction is the best practice
#   or if it would be possible to encorporate
#   multiple models that cover unique regions
#   of the protein / combining the results.
# - Should add warnings.

# Imports
import time, subprocess, os, glob, sys
from collections import defaultdict
from mjm_parsers import parse_dictionary_list, unzip_res_range

# Identify input UniProt info / Interactions
# files from input parameters, otherwise
# use defaults
interactomes_regex = sys.argv[2]
uniprot_info_file = sys.argv[1]

# Sifts Mapping Files
# Static copy from resources
sifts_file = '../pdb/pdbresiduemapping.txt'

# List of non-redundant ModBase Models
# It looks like this is using a static version
# from earlier, and that this should be updated
# to include the most recent Eclair results
# as follows...
mb_models_file = '../modbase/models/parsed/select_modbase_models.txt'
#mb_models_file = '/home/mjm659/ires_ml/site_data/models_interaction/zdock/select_modbase_models.txt'

# List of excluded PDB structures
excluded_pdbs_file = '../pdb/excluded_pdbs.txt'

# Fixed Output file
output_file = 'mixed_docking_set.txt'

# Obtain list of all the interactions
# in the input file
interactomes = set()
for f in glob.glob(interactomes_regex):
	print 'parsing', f, '...'
	interactomes.update(set([tuple(sorted(l.strip().split('\t')[:2])) for l in open(f)][1:])  )

# Obtain list of PDB structures
# to be excluded
excluded_pdbs = defaultdict(set)
for e in parse_dictionary_list(excluded_pdbs_file):
	
	# Refer to previous scripts for comments
	if e['hasCC']=='N':
		continue
	if e['excludedPDBs'] == '':
		continue
	
	interaction = tuple(sorted([e['UniProtA'], e['UniProtB']]))  #should already be sorted, but doesn't hurt to do it again
	excluded_pdbs[interaction] = set(e['excludedPDBs'].split(';'))


# Read UniProt Info file as dictionary
uniprot2info = dict((e['id'], e) for e in parse_dictionary_list(uniprot_info_file))

# Read SIFTS Mapping as dicitonary
sifts_data = parse_dictionary_list(sifts_file)

# Read ModBase Mapping as dictionary
modbase_data = parse_dictionary_list(mb_models_file)


# Iterate over every ModBase Model to generate
# a dictionary mapping UniProt ID to list
# of valid ModBase models to use for docking
uniprot2models = defaultdict(list)
for e in modbase_data:
	# Skip entries that do not map to any of the
	# UniProts in the UniProt Info file
	# Warning?
	if e['uniprot'] not in uniprot2info:
		continue
	
	# Skip entries with low quality
	if float(e['modpipe_quality_score']) < 1.1:
		continue
	
	# Obtain set of residues covered by this model
	covered_res = set(range(int(e['target_begin']), int(e['target_end'])+1))
	coverage = len(covered_res)
	
	# Skip entries with insufficiently low coverage
	# (raw counts)
	if coverage < 50:
		continue
	
	coverage_perc = coverage / float(uniprot2info[e['uniprot']]['length'])
	    
	# Skip entries with insufficiently low coverage
	# (as a proportion of total UniProt length)
	if coverage_perc < 0.50:
		continue
	
	# Add model to list for docking
	uniprot2models[e['uniprot']].append((coverage, coverage_perc, covered_res, e['modbase_modelID'], '', e['template_pdb'].upper(), 'MB'))


# Iterate over every Sifts mapping to generate
# a dictionary mapping UniProt ID to list
# of valid PDB structures to use for docking
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
	
	coverage_perc = coverage / float(uniprot2info[e['UniProt']]['length'])
	
	# Skip entries with insufficiently low coverage
	# (as a proportion of total UniProt length)
	if coverage_perc < 0.50:
		continue
	
	# Add structure to list for docking
	uniprot2models[e['UniProt']].append((coverage, coverage_perc, covered_res, e['PDB'], e['Chain'], e['PDB'], 'PDB'))


# Prioritize Structure option by coverage
for u in uniprot2models:
	uniprot2models[u].sort(reverse=True)

# Open output for writting
# Why is there no determination of
# interactions that had been previously
# calculated?
output = open(output_file, 'w')
output.write('\t'.join(['ProtA', 'ProtB', 'SubA', 'SubB', 'CovA', 'CovB', 'pCovA', 'pCovB']) + '\n')

# Iterate over each interaction to find the
# best structure / chain to use for each protein
for interaction in interactomes:
	# Obtain UniProt IDs
	p1, p2 = interaction
	
	# Restrict valid ModBase models to use for each
	# protein by removing any models that use excluded
	# PDB structures as their templates
	p1_subunits = [mb for mb in uniprot2models[p1] if mb[-2] not in excluded_pdbs[interaction]]
	p2_subunits = [mb for mb in uniprot2models[p2] if mb[-2] not in excluded_pdbs[interaction]]
	
	# Skip entries that have no valid models left
	if len(p1_subunits) < 1 or len(p2_subunits) < 1:
		continue
	
	# Figure out all of the valid structure combinations
	# that are from difference sources (i.e. cannot look
	# at any combinations that use two ModBase models or
	# two PDB models since that would make this docking set
	# redundant with the previous docking sets)
	subunit_combinations = []
	for s1 in p1_subunits:
		for s2 in p2_subunits:
			
			# Skip structure pairs that are from the same source
			if s1[-1] == s2[-1]:
				continue
			
			# Add structure pair
			subunit_combinations.append((s1[1]+s2[1], s1, s2))
	
	# Skip entries that have no valid structure combinations
	if len(subunit_combinations) == 0:
		continue
	
	# Prioritize Structure option by cummulative
	# coverage
	# Note: This still has the effect of selecting the
	# highest coverage ModBase for P1 matched with the
	# highest coverage PDB for P2 or vice versa
	subunit_combinations.sort(reverse=True)
	
	# Use the best strucutre pair (based on
	# cummulative coverage)
	subA = subunit_combinations[0][1]
	subB = subunit_combinations[0][2]
	
	# For PDB structures, reformat structure names
	# to include the chain identifiers
	if subA[4] != '':
		subA_struct = subA[3] + '_' + subA[4]
	else:
		subA_struct = subA[3]
	
	if subB[4] != '':
		subB_struct = subB[3] + '_' + subB[4]
	else:
		subB_struct = subB[3]
	
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
	output.write('%s\t%s\t%s\t%s\t%s\t%s\t%.4f\t%.4f\n' %(p1, p2, subA_struct, subB_struct, subA[0], subB[0], subA[1], subB[1]))

# Close Output
output.close()
