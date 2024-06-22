# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script compiles the dist3D results calculated
# previously in the calc_X_dist3d.py scripts. The
# values are compiled into four final features
# to be used in the Eclair predictions; 1) The
# average dist3D per residue as reported over all
# of the docking results, 2) The highest dist3D per
# residue as reported, over all of the docking results,
# 3) The lowest dist3D per residue as reported, over
# all of the docking results, 4) The dist3D per
# residue as reported by the best docking result.
#
# This script does some very weird stuff dealing with
# the interactions that have co-crystal structures /
# the interactions that do not have co-crystal structures
# apparently differently (but I need to confirm).
# Basically, there is a lot going on her that I did
# not understand on the first pass through, and there
# are also a number of spots where I have identified
# things that look like there are not right.

# Expected Outcomes:
# The output files should be updated to include
# the final calculated Zdock dist3D based features
# for all of the interactions included in the input
# file.

# Known Bugs:
# - Includes static copy of ires_perppi and
#   interactomes from our resources folder
# - It looks like this script explicitly refuses
#   to calculate these features on any interaction
#   pairs that have evidence of being a true
#   interaction based on co-crystal or HINT (?)
#   evidence.
# - In the part that randomly only looks at co-crystal
#   structure information, it aggregates the dist3D
#   values for homodimers by taking the MAXIUMUM
#   rather than the MIMUMUM distance. This differs from
#   the way that it was first done on the PREDICTION
#   set. I don't know if these values actually ever
#   get used or not.
# - Conceptually I have no clue what the purpose of
#   separating out the CC vs PREDICTED interactions
#   in the beginning 66% of this script it. It all
#   seems like it can be cut out to me, but I still
#   need to confirm that there wasn't an informed
#   purpose for doing all of it.
# - Am concerned that all of the outputs are opened
#   for writting, but it looks like the "already
#   computed" results are never re-written. I expect
#   that all old results for this set may be lost?
#   Based on looking at the length of the ouput
#   files, I don't see how that could be happening.

# Imports
from mjm_parsers import parse_dictionary_list, unzip_res_range
from mjm_ml_methods import normalize
from collections import defaultdict
import numpy as np
import random, glob, os, sys, subprocess

'''Parse a single docking experiment per interaction in priority order: (1) PDB docking, (2) Mixed Docking, and (3) ModBase docking.'''

# Helper function to "clean" the dSASA files
# (and also the dist3d files). I don't think
# this is a good way to go about this and
# if a problem really does exist if should be
# handled in the code that generated these files
# not after the fact.
def clean_dsasa_files():
    # For some reason, pdb_dsasa_file has empty lines which screws up parse_dictionary_list.
    # We need to remove empty lines first so the script doesn't crash.
    subprocess.call(["sed", "-i", '/^$/d', "/home/adr66/eclair/data/zdock/ires/dSASA_pdb_docking.txt"])
    subprocess.call(["sed", "-i", '/^$/d', "/home/adr66/eclair/data/zdock/ires/dist3d_mb_docking.txt"])
    subprocess.call(["sed", "-i", '/^$/d', "/home/adr66/eclair/data/zdock/ires/dist3d_pdb_docking.txt"])

# Set some threshold based on the input parameters
# otherwise set to the default.
# For our purposes, it looks like we co not use
# either of these filteres.
if len(sys.argv) > 2:
	docking_score_cutoff = int(sys.argv[2])
	ires_cutoff = 0
else:
	docking_score_cutoff = 0
	ires_cutoff = 0

# Read in Mixed, ModBase, and PDB dist3D results
# I don't know why these are called dSASA in the
# variable names.
mixed_dsasa_file = 'ires/dist3d_mixed_docking.txt'
mb_dsasa_file = 'ires/dist3d_mb_docking.txt'
pdb_dsasa_file = 'ires/dist3d_pdb_docking.txt'

# Read in Mixed, ModBase, and PDB ires results
pdb_ires_file = 'ires/IRES_pdb_docking.txt'
mixed_ires_file = 'ires/IRES_mixed_docking.txt'
mb_ires_file = 'ires/IRES_mb_docking.txt'

# Generate output prefix (we will add
# extension to indicate the type of
# aggregation being performed on each
# of the final features)
output_prefix = '/home/adr66/eclair/features/per_feature/zdock_dist3d_PRIORITY_%icut' %(docking_score_cutoff)

# Generate summary output file
summary_file = 'ires/priority_dist3d_zdock_%icut_%ires.txt' %(docking_score_cutoff, ires_cutoff)

# Generate dictionary for storing dist3D values
# for each interaction.
# Interaction Pair --> Docking Pair --> dist3D values
interaction2sasa = defaultdict(lambda: defaultdict(lambda: ([], [])))

# Store all UniProts seen in Each PDB
pdb2uniprots = defaultdict(set)

# Obtain interactions file from input
interactomes_dir = sys.argv[1]

# Fixed file for pre-calculated ires based on
# cocrystal structure?
# This seems to be a static copy from our resources
# folder
cocrystal_file = '../pdb/ires_perppi_alltax.txt'

# Fixed list of files for HQ interactomes?
# This seems to be a static copy from our resources
# folder
# I have not tracked down the origin of these
# files yet
interactome3d_files = glob.glob('../interactomes/*HQ.txt')

# Obtain list of interactions included in input
# file
interactomes = set()
for f in glob.glob(interactomes_dir):
    interactomes.update(  set([tuple(sorted(l.strip().split('\t')[:2])) for l in open(f)])  )

# Obtain list of reported co-crystal structures
# provided in our ires_perppi resource
cc_interactions = set()
for e in parse_dictionary_list(cocrystal_file):
    cc_interactions.add((e['UniProtA'], e['UniProtB']))

# Obtain list of reported interactions
# provided in some interactome folders
# that I have not tracked down the origin of
int3d_interactions = set()
for f in interactome3d_files:
    int3d_interactions.update(  set([tuple(sorted(l.strip().split('\t')[:2])) for l in open(f)][1:])   )

# Only make predictions on cases that are
# not included in either of these interaction
# files. I think this should be removed from
# the final Eclair-to-go
prediction_interactions = interactomes - cc_interactions - int3d_interactions
print len(prediction_interactions)

# Figure out which interactions have already
# had this feature calculated
# I think this might be broken because it looks
# like it is joining the ENTIRE line rather than
# just the interaction?
# This would explain further confusions that I have
# about the way already_calculated seems to be
# working in this script
already_computed = set()
with open(output_prefix + '_max.txt','r') as f:
	for l in f:
		already_computed.add('_'.join(sorted(l.strip().split('\t'))))

# Dictionary to map docking pair + rank
# to number of interface residues
subunits2ires = {}  #(subA, subB, rank) -> (num_iresA, num_iresB)

# Iterate over all the interface residue
# results for the PDB, Mixed, and ModBase
# docking sets
for e in parse_dictionary_list(pdb_ires_file) + parse_dictionary_list(mixed_ires_file) + parse_dictionary_list(mb_ires_file):
	
	# Report the number of interface residues
	# identified by this dock
	iresA = len(unzip_res_range(e['UniProtA_IRES']))
	iresB = len(unzip_res_range(e['UniProtB_IRES']))
	
	# Save the dock / ires count to the dictionary
	subunits2ires[(e['SubunitA'], e['SubunitB'], e['Rank'])] = (iresA, iresB)
    
# Clean the input files to remove unexpected
# blank lines
clean_dsasa_files()

# Iterate through all of the dist3D results
# for the PDB, Mixed, and ModBase docking sets
# Apparently this is a "first pass" to "collect
# best models for all PREDICTION interactions",
# and "find ratios of docking evidence types"
for e in parse_dictionary_list(pdb_dsasa_file) + parse_dictionary_list(mixed_dsasa_file) + parse_dictionary_list(mb_dsasa_file):
	# Obtain Interaction Pair
	interaction = (e['UniProtA'], e['UniProtB'])
	
	# Obtain Docking Pair
	dock_pair = (e['SubunitA'], e['SubunitB'])
	
	# Obtain Docking Score
	zdock_score = float(e['ZDOCK_Score'])
	
	# Skip Interactions that either have interaction
	# evidence or have already had these features
	# calculated
	if interaction not in prediction_interactions or '_'.join(sorted(interaction)) in already_computed:
		continue
	
	# Skip over Docking results for interactions that
	# we have already saved docking information for.
	# In this way we will only use features calculated
	# from ONE docking simulation, and we will prioritize
	# using results from the PDB, then Mixed, then ModBase
	# docking sets, only going into the lower confidence
	# docking sets when nothing is available in a higher set.
	if interaction in interaction2sasa and dock_pair not in interaction2sasa[interaction]:
		continue
	
	# Skip Docking results that have insufficiently
	# low docking score
	if zdock_score < docking_score_cutoff:
		continue
	
	# Skip Docking results that have insufficiently
	# few interface residues identified in them
	if ires_cutoff > 0:
		try:
			iresA, iresB = subunits2ires[(e['SubunitA'], e['SubunitB'], e['Rank'])]
			if iresA + iresB < ires_cutoff:
				continue
		except KeyError:
			continue
	
	# Parse the dist3D results for each protein
	dsasas1 = np.array([ float(r) if r != 'nan' else np.nan for r in e['UniProtA_dSASA'].split(';') ])
	dsasas2 = np.array([ float(r) if r != 'nan' else np.nan for r in e['UniProtB_dSASA'].split(';') ])
	
	# Skip entries that either have no distance (?)
	# or had no reported distance values
	# These entries reflect poor docking results?
	if all((dsasas1==0.0) | np.isnan(dsasas1)) or all((dsasas2==0.0) | np.isnan(dsasas2)):
		continue
	
	# Special case for Homodimers since they
	# need to accomodate the same values for
	# each side of the interaction
	if e['UniProtA'] == e['UniProtB']:
		# Obtain list of dist3D values from both chains
		both_dsasas = [dsasas1, dsasas2]
		
		# Generate a mask for all the NaN values
		nan_mask = np.ma.masked_array(both_dsasas, np.isnan(both_dsasas))
		
		# Obtain minimum dist3D per residue taking into
		# account information from both versions of
		# the protein
		dsasa_min = np.min(nan_mask, axis=0)
		
		# Add back NaN values for any positions that
		# had no values
		dsasa_min = np.array([dsasa_min.data[r] if dsasa_min.mask[r]==False else np.nan for r in range(len(dsasa_min.data))])
		
		# Use this cummulative dist3D value for both
		# sides of the interaction
		interaction2sasa[interaction][dock_pair][0].append(normalize(dsasa_min))
		interaction2sasa[interaction][dock_pair][1].append(normalize(dsasa_min))
	# Otherwise, just save the raw dist3D values for each
	# of the two unique proteins
	else:	
		interaction2sasa[interaction][dock_pair][0].append(normalize(dsasas1))
		interaction2sasa[interaction][dock_pair][1].append(normalize(dsasas2))


# Calculate fraction of interactions that use
# each docking set source
# I think this is just for our own diagnostic
# purposes?
source_fracs = [0,0,0]

# Iterate over each of the selected Docking results
# to use for each interaction
for interaction in interaction2sasa.keys():
	
	# Make sure there was in fact only one Docking
	# result selected
	assert len(interaction2sasa[interaction]) == 1
	
	# Obtain the names of the two structures involved
	# in the Docking Pair
	subA, subB = interaction2sasa[interaction].keys()[0]
	
	# Increment Counter for appropriate Source
	if len(subA) < 32 and len(subB) < 32:
		# PDB Docking Set
		source_fracs[0] += 1
	elif len(subA) == 32 and len(subB) == 32:
		# ModBase Docking Set
		source_fracs[2] += 1
	else:
		# Mixed Docking Set
		source_fracs[1] += 1

# Display statistics
print 'PREDICT_INTERACTIONS:', sum(source_fracs)
if sum(source_fracs) != 0:
	source_fracs = [float(t)/sum(source_fracs) for t in source_fracs]
	print 'PREDICT RATIOS:', source_fracs
print

# Dictionary to map co-crystal structure
# results --> docking pair + rank --> number of
# interface residues
cc_interaction2sasa = defaultdict(lambda: defaultdict(lambda: ([], [])))

# Clean the input files to remove unexpected
# blank lines (Again?)
clean_dsasa_files()

# Iterate through all of the dist3D results
# for the PDB, Mixed, and ModBase docking sets
# Appearently this is a "second pass" to "collect
# ALL models for all CC interactions", and "try
# to match ratios seen in predicted"
for e in parse_dictionary_list(pdb_dsasa_file) + parse_dictionary_list(mixed_dsasa_file) + parse_dictionary_list(mb_dsasa_file):

	# Obtain Interaction Pair
	interaction = (e['UniProtA'], e['UniProtB'])
	
	# Obtain Docking Pair
	dock_pair = (e['SubunitA'], e['SubunitB'])
	
	# Obtain Docking Pair
	zdock_score = float(e['ZDOCK_Score'])
	
	# Skip Interaction that either fo not have interaction
	# evidence or have already had these features
	# calculated
	if interaction in prediction_interactions or '_'.join(sorted(interaction)) in already_computed:
		continue
	
	# Skip Docking results that have insufficiently
	# low docking score
	if zdock_score < docking_score_cutoff:
		continue
	
	# Skip Docking results that have insufficiently
	# few interface residues identified in them
	if ires_cutoff > 0:
		try:
			iresA, iresB = subunits2ires[(e['SubunitA'], e['SubunitB'], e['Rank'])]
			if iresA + iresB < ires_cutoff:
				continue
		except KeyError:
			continue
	
	# Parse the dist3D results for each protein
	dsasas1 = np.array([ float(r) if r != 'nan' else np.nan for r in e['UniProtA_dSASA'].split(';') ])
	dsasas2 = np.array([ float(r) if r != 'nan' else np.nan for r in e['UniProtB_dSASA'].split(';') ])
	
	# Skip entries that either have no distance (?)
	# or had no reported distance values
	# These entries reflect poor docking results?
	if all((dsasas1==0.0) | np.isnan(dsasas1)) or all((dsasas2==0.0) | np.isnan(dsasas2)):
		continue  #poor docking
	
	# Special case for Homodimers since they
	# need to accomodate the same values for
	# each side of the interaction
	if e['UniProtA'] == e['UniProtB']:
		# Obtain list of dist3D values from both chains
		both_dsasas = [dsasas1, dsasas2]
		
		# Generate a mast for all the NaN calues
		nan_mask = np.ma.masked_array(both_dsasas, np.isnan(both_dsasas))
		
		# Obtain maximum dist3D per residue taking into
		# account information from both versions of
		# the protein
		# NOTE: This is different than how it was calculated
		# above. I am guessing there is an accidental failure
		# to swap from MAX to MIN (since this code was recycled
		# from the dSASA step that used maximum originally)
		# NOTE: I am not sure if the values that are being saved
		# here even matter.
		dsasa_max = np.max(nan_mask, axis=0)
		
		# Add back NaN values for any positions that
		# had no values
		dsasa_max = np.array([dsasa_max.data[r] if dsasa_max.mask[r]==False else np.nan for r in range(len(dsasa_max.data))])
		
		# Use this cummulative dist3D value for both
		# sides of the interaction
		cc_interaction2sasa[interaction][dock_pair][0].append(dsasa_max)
		cc_interaction2sasa[interaction][dock_pair][1].append(dsasa_max)
	# Otherwise, just save teh raw dist3D values for each
	# of the two unique proteins
	else:	
		cc_interaction2sasa[interaction][dock_pair][0].append(dsasas1)
		cc_interaction2sasa[interaction][dock_pair][1].append(dsasas2)

# Adjust the ratios?
# It looks like this is generating a normalized
# count for how many of each type we should
# expect to see in the CC dataset based on
# what we saw in the PREDICTION dataset?
total_interactions = len(cc_interaction2sasa)
target_pdb = total_interactions * source_fracs[0]
target_mixed = total_interactions * source_fracs[1]
target_mb = total_interactions * source_fracs[2]

# Mapping of Interaction --> Docking Pair used
cc_subunits = {}

# Another dictionary?
# That literally never shows up
# in this code again...
# WTF?
final_cc_interaction2sasa = defaultdict(lambda: defaultdict(lambda: ([], []))) #(p1,p2) -> (subA, subB) -> ([dsasas1...], [dsasas2...])

# Iterate over all of the CC interactions sorted
# by the number of sources that they had available
# Has an added element of randomness for sorting among
# the tied interaction?
# My current guess is that this is somehow feeding the
# results from the cc_interaction2sasa subset back into
# the original interaction2sasa dicitonary (since we
# presumably do still want to calculate features on
# these interactions as well?)
for interaction in sorted(cc_interaction2sasa.keys(), key=lambda i: len(cc_interaction2sasa[i].keys())+random.random()):
	# If there is only one Docking Result we have to
	# use it
	if len(cc_interaction2sasa[interaction].keys()) == 1:
		cc_subunits[interaction] = cc_interaction2sasa[interaction].keys()[0]
		# Should this be final_cc_interaction2sasa?
		interaction2sasa[interaction] = cc_interaction2sasa[interaction]
	else:
		# Otherwise break the docking options down into
		# counts for evidence in from each of the three
		# docking sets based on the length of their components
		pdb_subunits = len([s for s in cc_subunits.values() if len(s[0]) + len(s[1]) < 16])
		mixed_subunits = len([s for s in cc_subunits.values() if len(s[0]) + len(s[1]) > 32 and len(s[0]) + len(s[1]) < 64])
		mb_subunits = len([s for s in cc_subunits.values() if len(s[0]) + len(s[1]) == 64])
		
		# Figure out the difference between the expected
		# number of results form that source and the
		# actual number of results from that source?
		# Then this list is sorted in reverse in some
		# meaningful way that I don't understand?
		diffs = sorted([(max(target_pdb-pdb_subunits, 0), 'PDB'), (max(target_mixed-mixed_subunits, 0), 'MIXED'), (max(target_mb-mb_subunits, 0), 'MB')], reverse=True)
		
		# Break the docking pairs out into three lists
		# based on which docking set the pair comes from
		avail_pdb_subunits = [s for s in cc_interaction2sasa[interaction] if len(s[0]) + len(s[1]) < 16]
		avail_mixed_subunits = [s for s in cc_interaction2sasa[interaction] if len(s[0]) + len(s[1]) > 32 and len(s[0]) + len(s[1]) < 64]
		avail_mb_subunits = [s for s in cc_interaction2sasa[interaction] if len(s[0]) + len(s[1]) == 64]
		
		# Now we iterate through each category and set the
		# cc_subunits and interaction2sasa value?
		# I think this just has the effect of actually selecting
		# whichever category get sorted last in the magical
		# diffs sort, and using that as the selected Docking
		# Pair?
		for d in diffs:
			
			if d[1] == 'PDB' and len(avail_pdb_subunits) != 0:
				cc_subunits[interaction] = avail_pdb_subunits[0]
				interaction2sasa[interaction][avail_pdb_subunits[0]] = cc_interaction2sasa[interaction][avail_pdb_subunits[0]]
				#~ print avail_pdb_subunits[0]
				break
			elif d[1] == 'MIXED' and len(avail_mixed_subunits) != 0:
				cc_subunits[interaction] = avail_mixed_subunits[0]
				interaction2sasa[interaction][avail_mixed_subunits[0]] = cc_interaction2sasa[interaction][avail_mixed_subunits[0]]
				#~ print avail_mixed_subunits[0]
				break
			elif d[1] == 'MB' and len(avail_mb_subunits) != 0:
				cc_subunits[interaction] = avail_mb_subunits[0]
				interaction2sasa[interaction][avail_mb_subunits[0]] = cc_interaction2sasa[interaction][avail_mb_subunits[0]]
				#~ print avail_mb_subunits[0]
				break

# Calculate fraction of CC interactions that use
# each docking set source
# I think this is just for our own diagnostic
# purposes?
cc_fracs = [0,0,0]

# Iterate over each of the rselected CC Docking
# results to use for each interaction
for interaction in cc_subunits.keys():
	
	# Make sure there was in fact only one Docking
	# result selected
	assert len(interaction2sasa[interaction]) == 1
	
	# Obtain the names of the two structures involved
	# in the Docking Pair
	subA, subB = interaction2sasa[interaction].keys()[0]
	
	# Increment Counter for appropriate Source
	if len(subA) < 32 and len(subB) < 32:
		# PDB Docking Set
		cc_fracs[0] += 1
	elif len(subA) == 32 and len(subB) == 32:
		# ModBase Docking Set
		cc_fracs[2] += 1
	else:
		# Mixed Docking Set
		cc_fracs[1] += 1

# Display statistics
print 'CC_INTERACTIONS:', sum(cc_fracs)
cc_fracs = [float(t)/sum(cc_fracs) for t in cc_fracs]
print 'CC RATIOS:', cc_fracs
print


print 'writing feature files...'

# Open output files for writting
# This seems like it would overwrite
# all of the "already computed" values
# Note: already_computed is not generated
# correctly when it reads the output file
# since it uses the full line rather than
# just the columns corresponding to the
# interaction pair
#
# We report four output features...
#
# 1) avg - The average dist3D per residue as reported
#           over all of the docking results
#
# 2) max - The highest dist3D per residue as reported
#          over all of the docking results
#
# 3) min - The lowest dist3D per residue as reported
#          over all of the docking results
#
# 4) top1 - The dist3D per residue as reported by the
#           best docking result
#
output_avg = open(output_prefix + '_avg.txt', 'w')
output_max = open(output_prefix + '_max.txt', 'w')
output_min = open(output_prefix + '_min.txt', 'w')
output_top1 = open(output_prefix + '_top1.txt', 'w')

# Open the output summary file for writting
output_summary = open(summary_file, 'w')

# Iterate through each Interaction
for i in interaction2sasa:
	# Obtain UniProt IDs
	p1, p2 = i
	
	# Empty lists for saving the dist3D values of each protein
	p1_dsasas, p2_dsasas = [], []
	
	# Write Interaction Pair / Docking Pair selected to
	# the summery file
	output_summary.write('%s\t%s\t%s\t%s\n' %(p1, p2, interaction2sasa[(p1,p2)].keys()[0][0], interaction2sasa[(p1,p2)].keys()[0][1]))
	
	# Iterate over all of the valid docks that had dSASA values for
	# this interaction
	# Note: This is to say the top 10 docking conformations for
	# the single Docking Pair being used?
	for docki, dock_pair in enumerate(interaction2sasa[i]):
		
		# If this is the top ranked docking result
		# then we save these output for the top1
		# feature
		if docki == 0:   #first pair of subunits docked together (if there are more than one pair of subunits selected for this interaction)
			# Write output files
			# Format is...
			# P1	P2	Feat1	Feat2
			# Where Feat1 and Feat2 are semi-colon delimited
			# lists of feature values for each residue
			output_top1.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in interaction2sasa[i][dock_pair][0][0]]), ';'.join(['%.2f' %(res) for res in interaction2sasa[i][dock_pair][1][0]])))
		
		# Add the dist3D for this docking result to the lsit
		p1_dsasas += interaction2sasa[i][dock_pair][0]
		p2_dsasas += interaction2sasa[i][dock_pair][1]
	
	# Calculate features for protein 1
	# Obtain masked np array that covers all of the
	# dist3D values accross all docks / ignores NaNs
	#
	# rows = DOcking result
	# columns = Residue position		
	mdat = np.ma.masked_array(p1_dsasas, np.isnan(p1_dsasas))
	
	# Obtain the mean / max / min of dist3D value per
	# residue accross all of the docking results / fill
	# with NaN anywhere there was no coverage of a residue
	p1_means = np.mean(mdat, axis=0)
	p1_means = [p1_means.data[r] if p1_means.mask[r]==False else np.nan for r in range(len(p1_means.data))]
	p1_max = np.max(mdat, axis=0)
	p1_max = [p1_max.data[r] if p1_max.mask[r]==False else np.nan for r in range(len(p1_max.data))]
	p1_min = np.min(mdat, axis=0)
	p1_min = [p1_min.data[r] if p1_min.mask[r]==False else np.nan for r in range(len(p1_min.data))]
	
	# Calculate features for protein 2
	# Obtain masked np array that covers all of the
	# dist3D values accross all docks / ignores NaNs
	#
	# rows = DOcking result
	# columns = Residue position		
	mdat = np.ma.masked_array(p2_dsasas, np.isnan(p2_dsasas))

	# Obtain the mean / max / min of dist3D value per
	# residue accross all of the docking results / fill
	# with NaN anywhere there was no coverage of a residue
	p2_means = np.mean(mdat, axis=0)
	p2_means = [p2_means.data[r] if p2_means.mask[r]==False else np.nan for r in range(len(p2_means.data))]
	p2_max = np.max(mdat, axis=0)
	p2_max = [p2_max.data[r] if p2_max.mask[r]==False else np.nan for r in range(len(p2_max.data))]
	p2_min = np.min(mdat, axis=0)
	p2_min = [p2_min.data[r] if p2_min.mask[r]==False else np.nan for r in range(len(p2_min.data))]
	
	# Write output files
	# Format is...
	# P1	P2	Feat1	Feat2
	# Where Feat1 and Feat2 are semi-colon delimited
	# lists of feature values for each residue
	output_avg.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in p1_means]), ';'.join(['%.2f' %(res) for res in p2_means])))
	output_max.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in p1_max]), ';'.join(['%.2f' %(res) for res in p2_max])))
	output_min.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in p1_min]), ';'.join(['%.2f' %(res) for res in p2_min])))

# Close Outputs
output_avg.close()
output_max.close()
output_min.close()
output_top1.close()
output_summary.close()
