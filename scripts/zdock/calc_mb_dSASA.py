# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script calculates the delta SASA for
# all of the docking results in the ModBase docking
# set. Delta SASA are calculated using mjm_tools
# irescalc.py. Wherever possible the script avoids
# recalculating delta SASA values for docking results
# that have already had their delta SASA values
# calculated.

# Expected Outcomes:
# The output file should be updated to include entries
# for delta SASA values of all residues for all of the docking
# results that were included in the ModBase docking
# set.

# Known Bugs:
# - irescalc.py is known to fail to calculate interface
#   residues on particularly large protein structures. This
#   failure was not originally explicily reported. I am
#   guessing that the way I corrected it to handle reporting
#   this error may cause this script to break.
# - Bit of confusion on the bit that converts from ModBase
#   positions to UniProt positions. I think this is just
#   checking to make sure that all the positions not covered
#   by the ModBase models are filled in as NaN but I am not
#   completely sure.

# Imports
import time, sys, os, glob, re, sys
from mjm_parsers import unzip_res_range, parse_dictionary_list
from mjm_tools import fetch_uniprot, xtool, time_elapsed

# Keep track of runtime
start_time = time.time()

# Path to docking results Directory
docking_results_dir = 'modbase_docked_models'

# List of non-redundant ModBase models
# selected for each UniProt
modbase_info_file = '../modbase/models/parsed/select_modbase_models.txt'

# Fixed path to output file
output_file = 'ires/dSASA_mb_docking.txt'

# Obtain list of all interactions in input file
interlist = [sorted(l.strip().split('\t')[:2]) for l in open(sys.argv[1])]

# Obtain list of all docking results
docked_model_files = sorted(glob.glob(os.path.join(docking_results_dir, '*.pdb')))

# Read dictionary of ModBase models
modbase2info = dict((e['modbase_modelID'].upper(), e) for e in parse_dictionary_list(modbase_info_file))

# Obtain list of all docking results that have
# already had their delta SASAs calculated
# Here we save the actual line
already_calculated = {}
if os.path.exists(output_file):
	for l in open(output_file).readlines()[1:]:
		already_calculated[l.split('\t')[5]] = l

print 'Found %i docked models in %s...' %(len(docked_model_files), docking_results_dir)

# Open output for writting
output = open(output_file, 'w')
output.write('\t'.join(['UniProtA', 'UniProtB', 'Rank', 'SubunitA', 'SubunitB', 'File', 'ZDOCK_Score', 'UniProtA_dSASA', 'UniProtB_dSASA']) + '\n')

# Iterate over each docking result
for fi, f in list(enumerate(docked_model_files)):
	# Status Update
	if fi%100==0:
		print '%s/%s complete' %(fi, len(docked_model_files))
	
	# For entries that have already had their interface
	# residues calculated, just re-write the previous
	# results
	if os.path.basename(f) in already_calculated:
		output.write(already_calculated[os.path.basename(f)])
		continue
	
	# Obtain the two Models that made up the dock and
	# the rank of the docking result
	subA, subB, _, num, _ = re.split('--|-|\.', os.path.basename(f))
	
	# Skip entries where either model is not included
	# in the ModBase info
	if subA not in modbase2info or subB not in modbase2info:
		continue
	
	# Obtain the UniProt IDs corresponding to these
	# models
	uniprotA = modbase2info[subA]['uniprot']
	uniprotB = modbase2info[subB]['uniprot']
	
	# Skip entries where the docking results do not
	# represent one of the interactions included in
	# the current input file
	if [uniprotA,uniprotB] not in interlist and [uniprotB,uniprotA] not in interlist:
		continue
	
	# Calculate delta SASA by calling mjm_tools
	# irescalc.py. Here we modify the thresholds for
	# defining and interface so that it reports the
	# dSASA values of all residues. The output looks like...
	#
	# HEAD:
	#
	# Chain	Pos	Resn	dSASA
	# A	233	PHE	53.1
	# B	123	ARG	72.31
	#
	out = [x.strip().split('\t') for x in xtool('irescalc.py', f, c1='A', c2='B', uSASA=0, dSASA=0, o='residue_stats').strip().split('\n')]
	
	# If there is no output, then there the calculation
	# is considered to have failed
	if out == [['']]:
		print 'irescalc failed: %s' %(f)
		continue
	
	# Split the dSASA values by chain / generate
	# dictionary (pos --> dSASA)
	A_dSASA = dict([(int(q[1]), q[3]) for q in out if q[0]=='A'])  #map UniProt residues to SASA values
	B_dSASA = dict([(int(q[1]), q[3]) for q in out if q[0]=='B'])  #map UniProt residues to SASA values
	
	# Filter the ModBase positions so that we only
	# included the positions in the UniProt that
	# were covered by the template. Add in NaN at any
	# positions in teh UniProt that were not covered.
	# I have to think about how this is done more
	# to confirm that I trust the way this is done.
	# I remember having questions about the way the 
	# ModBase positions were mapped originally.
	# Why did this mapping not need to be done previously
	# on the raw interface residues?
	A_dSASA = [str(A_dSASA[r]) if r in A_dSASA else 'nan' for r in range(1, int(modbase2info[subA]['target_length'])+1)]
	B_dSASA = [str(B_dSASA[r]) if r in B_dSASA else 'nan' for r in range(1, int(modbase2info[subB]['target_length'])+1)]
	
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
		zdock_score = open(f[:-7]+'.out').readlines()[int(num)+3].strip().split()[-1]
	# If we can't read the output file or the ouput file is
	# missing, just set the score to n/a
	except:
		zdock_score = 'n/a'
	
	# Write the output file. If the UniProt IDs are not
	# already sorted, flip all of the output columns so
	# that the output corresponds to the sorted interaction
	# pair. The output format is as follows...
	#
	# HEAD:
	# UniProtA UniProtB Rank  SubunitA SubunitB File  ZDOCK_Score UniProtA_dSASA UniProtB_dSASA
	# O00757   P01111   01 0006F41AE20CAED003C84CF1841CE2B7 35BFCA6537FC02352701FDE04AB07347 0006F41AE20CAED003C84CF1841CE2B7--35BFCA6537FC02352701FDE04AB07347--ZDOCK-01.pdb 4.808 nan;nan;nan;nan;nan;nan;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0....
	# O00757   P01111   02 0006F41AE20CAED003C84CF1841CE2B7 35BFCA6537FC02352701FDE04AB07347 0006F41AE20CAED003C84CF1841CE2B7--35BFCA6537FC02352701FDE04AB07347--ZDOCK-02.pdb 4.647 nan;nan;nan;nan;nan;nan;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0....
	# O00757   P01111   03 0006F41AE20CAED003C84CF1841CE2B7 35BFCA6537FC02352701FDE04AB07347 0006F41AE20CAED003C84CF1841CE2B7--35BFCA6537FC02352701FDE04AB07347--ZDOCK-03.pdb 4.090 nan;nan;nan;nan;nan;nan;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0....
	# O00757   P01111   04 0006F41AE20CAED003C84CF1841CE2B7 35BFCA6537FC02352701FDE04AB07347 0006F41AE20CAED003C84CF1841CE2B7--35BFCA6537FC02352701FDE04AB07347--ZDOCK-04.pdb 3.738 nan;nan;nan;nan;nan;nan;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0....
	#
	if uniprotA > uniprotB:
		output.write('\t'.join([uniprotB, uniprotA, num, subB, subA, os.path.basename(f), zdock_score, ';'.join(B_dSASA), ';'.join(A_dSASA)])+'\n')
	else:
		output.write('\t'.join([uniprotA, uniprotB, num, subA, subB, os.path.basename(f), zdock_score, ';'.join(A_dSASA), ';'.join(B_dSASA)])+'\n')

# Close the Output
output.close()

# Final Status Update / Total runtime
print 'Finished calculating dSASA values for %s docked models: %s' %(len(docked_model_files), time_elapsed(start_time))