# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script compiles the dSASA results calculated
# previously in the calc_X_dSASA.py scripts. The
# values are compiled into three final features
# to be used in the Eclair predictions; 1) the
# average dSASA value per residue reported over
# all of the docking results, 2) the maximum dSASA
# value per residue reported over all of the docking
# results, 3) the dSASA value per residue reported
# by the highest ranking docking result.
#
# After reading through parse_priority_dist3d.py
# the first time, I have some additional confusions
# that I need to reconfirm here. Namely, are these
# results calculated base don only one docking
# simulation (compiling all of the 10 top conformations)
# or do they include potentially multiple (max 2?)
# docking simulations? If they include multiple,
# need to double check that the top ranked docking
# result is selected correctly.

# Expected Outcomes:
# The output files should be updated to include
# the final calculated Zdock dSASA based features
# for all of the interactions included in the input 
# file.

# Known Bugs:
# - Confusion why these features are not split
#   to have separate weights for ModBase vs. PDB
#   based docking results rather than joining
#   all of them.
# - It looks like the Mixed docking set is never
#   used to compile any of the final features?

# Imports
import sys, subprocess
from mjm_parsers import parse_dictionary_list
from collections import defaultdict
import numpy as np


# dSASA output files for ModBase and PDB
# docking sets
mb_dsasa_file = 'ires/dSASA_mb_docking.txt'
pdb_dsasa_file = 'ires/dSASA_pdb_docking.txt'

# Prefix name for the output file. Will
# add extension pertaining to the three
# aggregates of the feature (avg, max, top1)
output_prefix = '/home/adr66/eclair/features/per_feature/zdock_combined_0cut'

# Docking score threshold for docks to
# withold. It looks like we accept all
# docks? I don't know if a docking score
# can be negative
docking_score_cutoff = 0

# Generate dictionary to store mapping from
# UniProt Pair --> Docking result --> dSASA values
interaction2sasa = defaultdict(lambda: defaultdict(lambda: ([], [])))

# Store all UniProts seen in each PDB?
pdb2uniprots = defaultdict(set)

# Figure out which interaction pairs have already
# had this feature calculated
already_calculated = set()
with open(output_prefix + '_avg.txt','r') as f:
	for l in f:
		already_calculated.add('_'.join(sorted(l.strip().split('\t')[:2])))

print 'parsing raw dsasa files (pdb + modbase)...'

# I don't think this should happen? And if it does
# happen it should be fixed as by going back to the
# script responsible not with this work around.
#
# Previous Comment:
#
# For some reason, pdb_dsasa_file has empty lines
# which screws up parse_dictionary_list. We need to
# remove empty lines first so the script doesn't crash.
subprocess.call(["sed", "-i", '/^$/d', "/home/adr66/eclair/data/zdock/ires/dSASA_pdb_docking.txt"])

# Iterate over all elements in both the PDB and ModBase
# dSASA files
for e in parse_dictionary_list(pdb_dsasa_file) + parse_dictionary_list(mb_dsasa_file):
	
	# Identify the interaction pair
	interaction = (e['UniProtA'], e['UniProtB'])
	
	# Identify the dock pair
	dock_pair = (e['SubunitA'], e['SubunitB'])
	
	# Identify the docking score
	zdock_score = float(e['ZDOCK_Score'])
	
	# Skip entries that have already had this feature calculated
	if '_'.join(sorted(interaction)) in already_calculated:
		continue
	
	# Skip entries that do not meet the minimum docking
	# score criteria
	if zdock_score < docking_score_cutoff:
		continue
	
	# Parse out the dSASA values for each protein
	dsasas1 = np.array([ float(r) if r != 'nan' else np.nan for r in e['UniProtA_dSASA'].split(';') ])
	dsasas2 = np.array([ float(r) if r != 'nan' else np.nan for r in e['UniProtB_dSASA'].split(';') ])
	
	# Skip entries that either have no change in SASA
	# or had no reported dSASA values
	# These entries reflect poor docking results?
	if all((dsasas1==0.0) | np.isnan(dsasas1)) or all((dsasas2==0.0) | np.isnan(dsasas2)):
		continue
	
	# Special case for Homodimers since they
	# need to accomodate the same values for
	# each side of the interaction
	if e['UniProtA'] == e['UniProtB']:
		# Obtain list of dSASAs from both chains
		both_dsasas = [dsasas1, dsasas2]
		
		# Generate a mask for all the NaN values
		nan_mask = np.ma.masked_array(both_dsasas, np.isnan(both_dsasas))
		
		# Obtain maximum dSASA per residue taking into
		# account information from both versions of
		# the protein
		dsasa_max = np.max(nan_mask, axis=0)
		
		# Add back NaN values for any positions that
		# had no values
		dsasa_max = np.array([dsasa_max.data[r] if dsasa_max.mask[r]==False else np.nan for r in range(len(dsasa_max.data))])
		
		# Use this cummulative dSASA value for both
		# sides of the interaction
		interaction2sasa[interaction][dock_pair][0].append(dsasa_max)
		interaction2sasa[interaction][dock_pair][1].append(dsasa_max)
	# Otherwise, just save the raw dSASA values for each
	# of the two unique proteins
	else:
		interaction2sasa[interaction][dock_pair][0].append(dsasas1)
		interaction2sasa[interaction][dock_pair][1].append(dsasas2)

print 'writing feature files...'

# Open output files for writting
# We report three output features...
#
#  1) avg - The average dSASA per residue as reported
#           over all of the docking results
#
# 2) max - The highest dSASA per residue as reported
#          over all of the docking results
#
# 3) top1 - The dSASA per residue as reported by the
#           best docking result
#
output_avg = open(output_prefix + '_avg.txt', 'a')
output_max = open(output_prefix + '_max.txt', 'a')
output_top1 = open(output_prefix + '_top1.txt', 'a')

# Iterate over each interaction
for i in interaction2sasa:
	# Obtain UniProt IDs
	p1, p2 = i
	
	# Skip entries that have already had this feature calculated
	if '_'.join(sorted([p1,p2])) in already_calculated:
		continue
	
	# Empty lists for saving the dSASA values of each protein
	p1_dsasas, p2_dsasas = [], []
	
	# Iterate over all of the valid docks that had
	# dSASA values for this interaction
	for docki, dock_pair in enumerate(interaction2sasa[i]):
		# If this is the top ranked docking result
		# then we save these outputs for the top1
		# feature
		if docki == 0:
			# Write output files
			# Format is...
			# P1	P2	Feat1	Feat2
			# Where Feat1 and Feat2 are semi-colon delimited
			# lists of feature values for each residue
			output_top1.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in interaction2sasa[i][dock_pair][0][0]]), ';'.join(['%.2f' %(res) for res in interaction2sasa[i][dock_pair][1][0]])))
		
		# Add the dSASAs for this docking result to the lsit
		p1_dsasas += interaction2sasa[i][dock_pair][0]
		p2_dsasas += interaction2sasa[i][dock_pair][1]
	
	# Calculate features for protein 1
	# Obtain masked np array that covers all of the
	# dSASA values accross all docks / ignores NaNs
	#
	# rows = Docking result
	# columns = Residue position
	mdat = np.ma.masked_array(p1_dsasas, np.isnan(p1_dsasas))
	
	# Obtain the mean / max of dSASA value per residue
	# accross all of the docking results / fill with NaN
	# anywhere there was no coverage of a residue
	p1_means = np.mean(mdat, axis=0)
	p1_means = [p1_means.data[r] if p1_means.mask[r]==False else np.nan for r in range(len(p1_means.data))]
	p1_max = np.max(mdat, axis=0)
	p1_max = [p1_max.data[r] if p1_max.mask[r]==False else np.nan for r in range(len(p1_max.data))]
	
	
	# Calculate features for protein 2
	# Obtain masked np array that covers all of the
	# dSASA values accross all docks / ignores NaNs
	#
	# rows = Docking result
	# columns = Residue position
	mdat = np.ma.masked_array(p2_dsasas, np.isnan(p2_dsasas))
	
	# Obtain the mean / max of dSASA value per residue
	# accross all of the docking results / fill with NaN
	# anywhere there was no coverage of a residue
	p2_means = np.mean(mdat, axis=0)
	p2_means = [p2_means.data[r] if p2_means.mask[r]==False else np.nan for r in range(len(p2_means.data))]
	p2_max = np.max(mdat, axis=0)
	p2_max = [p2_max.data[r] if p2_max.mask[r]==False else np.nan for r in range(len(p2_max.data))]
	
	# Write output files
	# Format is...
	# P1	P2	Feat1	Feat2
	# Where Feat1 and Feat2 are semi-colon delimited
	# lists of feature values for each residue
	output_avg.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in p1_means]), ';'.join(['%.2f' %(res) for res in p2_means])))
	output_max.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in p1_max]), ';'.join(['%.2f' %(res) for res in p2_max])))

# Close Outputs
output_avg.close()
output_max.close()
output_top1.close()
