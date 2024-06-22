# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script fetches UniProt information for all
# UniProt IDs in a given input file of interactions
# and saved the information in a specified output
# file. Wherever possible, the script identifies
# UniProt ID information that has previously been
# fetched, so as to avoid fetching the same information
# multiple times.

# Expected Outcomes:
# All UniProt IDs in the input file should have their
# information fetched, and written to the output file.

# Known Bugs:
# - This script depends on mjm_tools fetch_uniprot to
#   obtain UniProt information. This function has been
#   known to break in the past as the uniprot API changes.
# - I am concerned about what may happen when UniProt entries
#   change since this code is set up to ignore re-fetching
#   information for previously handled IDs.

# Imports
import os, time, sys, math, warnings
from mjm_tools import fetch_uniprot

# Set up input / output file based on input parameters
# If no inputs are provided, default to fixed input
# and / or output file
if len(sys.argv) > 1:
	inters_file = sys.argv[1]
	output_file = sys.argv[2]
else:
	inters_file = 'all_uniprots.txt'
	output_file = 'uniprot_info.txt'

# List of attributes for UniProt into to be retrieved
attributes = ['id', 'entry name', 'reviewed', 'organism-id', 'lineage-id(SUPERKINGDOM)', 'length', 'genes', 'feature(ALTERNATIVE SEQUENCE)', 'protein names', 'sequence']

# Read in the UniProt IDs that already have information fetched
# in the output file.
already_fetched = {}
if os.path.exists(output_file):
	with open(output_file) as f:
		f.readline()
		for l in f:
			already_fetched[l.split('\t')[0]] = l

# Obtain a list of all of the Uniprot IDs in the input file
# that we wish to fetch information for
all_uniprots = set()
with open(inters_file, 'r') as ifile:
	for l in ifile:
		all_uniprots.update(l.strip().split('\t')[:2])

# Open up the output file for writting
# Would it be better to open this as an append rather than a re-write?
output = open(output_file, 'w')
output.write('\t'.join(attributes)+'\n')

# Write the previously fetched information for all UniProt IDs
# already in the output file
# This would be uneccessary if output file was opened as append?
for u in already_fetched:
	output.write(already_fetched[u])

# Download the remaining UniProt IDs in batches of 25 at a time
batch_size = 25
to_fetch = list(set(all_uniprots) - set(already_fetched.keys()))
num_batches = int(math.ceil(len(to_fetch)/float(batch_size)))
for batch_num in range(num_batches):
	# Generate batch
	uniprot_batch = to_fetch[batch_num*batch_size:(batch_num+1)*batch_size]
	
	# Status update
	print 'fetching %s-%s of %s: %s' %(batch_num*batch_size, (batch_num+1)*batch_size, len(to_fetch), str(uniprot_batch))
	
	# Fetch results
	result_dict = fetch_uniprot(uniprot_batch, attributes)
	
	# Check for failed IDs / Display a warning
	failed = set(uniprot_batch).difference(set(result_dict))
	if(len(failed) != 0):
		warnings.warn('Unable to fetch information for UniProts: %s' % ", ".join(failed))
	
	# Write results for all fetched UniProts
	for u in result_dict:
		output.write('\t'.join([result_dict[u][a] for a in attributes])+'\n')
	
	# Sleep? Presumably to avoid connection problems?
	time.sleep(1)

# Close output file
output.close()
