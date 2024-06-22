# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script calculates the JS Conservation for all positions
# of all UniProt IDs included in a given input file. JS
# Conservation is calculated using the script provided
# by Capra et all 2007, on a slimmed version of the MSAs
# generated in the previous psiblast_update.py script.
# The slimmed MSAs to my knowledge should not be any
# different than the original MSAs except that proteins
# are referred to by their taxon name rather than their
# UniProt ID. The script skips calculating this feature
# wherever a previous result already exists.

# Expected Outcomes:
# All UniProt IDs in the input file should have their
# JS Conservation added to the output file. Additionally
# a slimmed MSA should be created for all original MSA
# files. The exception is any UniProt ID with fewer than
# 50 alignments remaining in the slimmed MSA.

# Known Bugs:
# - No known bugs, but concerns including...
# - The way it skips previously calculated features
#   (no way of knowing if MSA has changed since the
#   last time it was calculated?)
# - Suggest throwing warnings for edge cases
# - Suggest a spot check to ensure JS Score is mapped
#   correctly
# - Possible efficiency improvements

# Imports
import os, glob, subprocess, sys
from mjm_tools import TimeoutCmd
from mjm_parsers import parse_fasta, write_fasta, parse_dictionary_list


#
# Step 1 - Generate Slimmed MSAs
#

# Find all of the MSAs
uniprot_MSAs = glob.glob('psiblast/*aligned.msa')

# Output directory for saving slimmed_MSAs
msa_output_dir = 'slimmed_msas/'

# Set up input file base on input parameter.
# If not input is provided, default to fixed
# input file
if len(sys.argv) > 1:
	uniprot_info_file = sys.argv[1]
else:
	uniprot_info_file = '../uniprot/uniprot_info.txt'

# Fixed output file
output_file = '/home/adr66/eclair/features/per_feature/JSconserve.txt'

# Dictionary mapping UniProt ID to length
uniprot2len = dict((e['id'], int(e['length'])) for e in parse_dictionary_list(uniprot_info_file))

### delete all previous msas
#os.system('rm %s' %(os.path.join(msa_output_dir, '*.msa')))
###

# Figure out which UniProts have already had their JS
# Conservation calculated
already_calculated = set()
for l in open(output_file,'r'):
	already_calculated.add(l.strip().split('\t')[0])

# Iterate over all MSAs (This is potentially redundant because
# we should KNOW which MSAs might be new based on the interactions)
for f in uniprot_MSAs:
	# Read MSA as dictionary / save keys
	ordered_uniprots, msa = parse_fasta(f)
	
	# Here we modify all of the keys in the MSA to be organism
	# names rather than UniProt_Organism as follows...
	#
	# tr|A0A024R087|A0A024R087_HUMAN --> HUMAN
	#
	# This potentially has the side effect (intened?) of removing
	# duplicate entries for the same organism. But I believe that
	# for the use case of Eclair, that should have already been
	# prohibited in the previous psiblast_update.py script.
	#
	# A previous comment also claimed that the last residue of
	# each alignment is removed because UCSC uses a 'Z" to
	# indicate the stop codon. But that line is clearly
	# commented out.
	#
	# We also appear to replace all Selenocysteine (U) with
	# gaps. I don't know why?
	for k in msa.keys():
		#~ msa[k.split('_')[-1]] = msa[k][:-1]   #
		msa[k.split('_')[-1]] = msa[k].replace('U', '-')
		del msa[k]
	
	#remove sequences with low residue coverage of human MSA
	#~ for k in msa.keys():
		#~ perc_missing = msa[k].count('-') / float(len(msa[k]))
		#~ 
		#~ if perc_missing > 0.5:
			#~ del msa[k]
	
	# Skip all MSA with insufficient aligned sequences
	# for good performance in DCA (DCA does not perform
	# well with few sequences in alignment)
	if len(msa) < 50:  #DCA does not perform well if there are few sequences in alignment
		continue
	
	# Write slimmed MSA file from the remaining alignment
	# entries with their modified keys
	write_fasta(msa, os.path.join(msa_output_dir, os.path.basename(f).replace('_aligned', '')))



#
# Step 2 - Calculate Jensen Shannon divergence
#

# Find all slimmed MSA files
slimmed_MSAs = sorted(glob.glob(os.path.join(msa_output_dir, '*.msa')))

print "calculating JS conservation for %s uniprots" %(len(slimmed_MSAs))

# Prepare appending to output_file
output_file = open(output_file, 'a')

# Move into directory with score_conservation script
os.chdir('conservation_code')

# Iterate over each slimed MSA
for i, f in enumerate(slimmed_MSAs):
	
	# Find the UniProt ID
	uniprot = os.path.basename(f)[:-4]
	
	# Skip UniProt IDs that have already had their
	# JS Conservation calculated
	if uniprot in already_calculated:
		continue
	
	print '%s (%s/%s)' %(uniprot, i+1, len(slimmed_MSAs))
	
	# Calculate JS Conservation using what I take to be
	# the script provided in the original publication
	output, err = subprocess.Popen('python score_conservation.py ../%s' %(f), shell=True, stdout=subprocess.PIPE).communicate()
	
	# Parse results by taking the score value for each position
	# Since we have already removed all of the gaps in the original
	# query, these positions should match up, directly with the UniProt
	# sequence positions, but I would encourage a spot-check at some point. 
	#
	# HEAD:
	#
	# # /home/adr66/eclair/data/conservation/slimmed_msas/P00791.msa -- js_divergence - window_size: 3 - background: blosum62 - seq. weighting: True - gap penalty: 1 - normalized: False
	# # align_column_number   score   column
	# 
	# 0       0.63365 MMMMM-M--M-MM-MMMMMM--V--MMM---MMMMMMMM-MMM--MMM-MMMMMMMM-MMMMMMM-MM-MMMMM---M-M--MMMM-MM--MMMM-MMMMMMMMMMMK--M-MMM-
	# 1       0.50970 KKQKK-K--R-KR-KKKRKR--M--NKK---RKKKKKKK-KKK--MKK-KKKKKKKK-KKKKKKK-RK-KRKRK---K-N--KKKKMKK--KKKK-KKKKRKKKKKKK--K-QKK-
	# 2       0.41300 WWSWW-L--C-WC-WVWWWC--W--R-W---CWWWWWWW-W-W--WWW-WWWWLGWW-WIWTWW--VW-WSWST---W-W--WWWWWWW--WWWW-WWWYCWWWWWWL--W-RWLL
	# 3       0.41837 LLLVL-I-LL-LL-LFMLLL--P--LLL---FLALLAFL-LLM--ALL-MLLLILLA-ALLLLLL-LL-LLALL---L-A--LLLLELL--LLLL-LLLLLLLLALLL--L-ALLL
	# 4       0.36371 LLLVL-A-VV-LV-LTVLLV--A--LLL---L-VVVVGV-LLV--VVL-VVLLAGVV-VLLVLLL-FW-LVVVV---V-L--VLVLRVL--LLVL-LWLIVLVVVLVL--V-VLLL
	#
	conservation_scores = [l.strip().split('\t')[1] for l in output.split('\n') if len(l) > 0 and l[0] != '#']
	
	# If the MSA UniProt ID does not match one of the UniProt
	# IDs in the original UniProt Info File, then something
	# has gone wrong. This should probaby be a warning
	if uniprot not in uniprot2len:
		print uniprot, 'not a valid UniProt ID'
		continue
	
	# If the number of JS Conservation score does not match
	# up with the known length of the UniProt ID, then
	# something has gone wronge. This should probably be a 
	# warning. Is this sufficient check that positions map
	# accurately?
	if len(conservation_scores) != uniprot2len[uniprot]:
		print uniprot, 'invalid number of residues'
		continue
	
	# Write the conservation scores in a semi-colon delimited
	# list
	output_file.write('%s\t%s\n' %(uniprot, ';'.join([str(s) for s in conservation_scores])))

# Close the output file
output_file.close()

# Return to starting directory
os.chdir('..')
