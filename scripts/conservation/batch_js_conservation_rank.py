# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script takes the previously calculated JS
# Conservation Scores from batch_js_conservation.py
# and ranks them from most conserved to least conserved.
# The script skips ranking the JS Scores for all
# UniProt IDs that have already had this feature
# calculated.

# Expected Outcomes:
# All UniProt IDs in the input file should have their
# ranked JS Scores written to the output file.

# Known Bugs:
# - No known bugs or concerns

# Imports
#from scipy.stats import rankdata
import numpy as np

# Use fixed input and output files
input_file = '/home/adr66/eclair/features/per_feature/JSconserve.txt'
output_file = '/home/adr66/eclair/features/per_feature/JSconserve_rank.txt'

# Determine all of the previously calculated
# features
output = open(output_file, 'a')
already_calculated = set()
for l in open(output_file, 'r'):
	already_calculated.add(l.strip().split('\t')[0])

# Takes a list of values, and returns their ranking
# relative to one another.
#
# e.g. [7, 2, 9, 4, 1] --> [4.0, 2.0, 5.0, 3.0, 1.0]
#
def rankdata(arr):
	# Make sure array is a np array
	arr = np.array(arr)
	
	# Sort array
	sorted_arr = np.sort(arr)
	
	# Sort indices
	sorted_idx = np.argsort(arr)
	
	# Track current rank
	curr = 1
	
	# Generate result space
	result = np.zeros(arr.shape)
	
	# Rank every position in the array
	for i in range(len(result)-1):
		# Assign the nth rank at the apropriate position
		result[sorted_idx[i]] = curr
		
		# Update the current rank unless there is a tie
		if sorted_arr[i] < sorted_arr[i+1]:
			curr = i+2
	
	# Set last rank
	result[sorted_idx[-1]] = curr
	
	# Return result
	return result

# Iterate over every UniProt with JS Conservation Scores
for l in open(input_file):
	# Obtain UniProt ID and semi-colon delimited list of JS Scores
	prot, js_scores = l.strip().split('\t')
	
	# Skip proteins that have already been calculated
	if prot in already_calculated:
		continue
	
	# Parse scores / multiply by -1 to sort from most conserved to least conserved
	js_scores = [-1*float(f) for f in js_scores.split(';')]
	
	# Rank the scores
	ranked_scores = rankdata(js_scores)
	
	# Append the ranked scores to the output file
	output.write(prot + '\t' + ';'.join([str(int(f)) for f in ranked_scores]) + '\n')

# Close the output
output.close()
