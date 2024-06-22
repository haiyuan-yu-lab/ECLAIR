# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script splits the DCA output files generated
# in the batch_dca.py script into a distinct output
# file for each the Direct Information and Mutual
# Information values of the original DCA.

# Expected Outcomes:
# All interactions in the interaction input file
# that had successfully generated DCA outputs should
# have these outputs split into new output files.

# Known Bugs:
# - Interacting protein tuples are not sorted here
#   in constructing the basename, but they were
#   previously sorted for saving the raw DCA
#   results. As a result, any interaction pairs
#   not correctly sorted in the raw input file
#   will never have their parsed DCA results
#   saved correctly.
#   + THIS IS A CONFIRMED BUG with a presumably
#     easy fix, but I haven't fixed it because
#     I don't want to mess with it, and always
#     just sort the interactions correctly for
#     the input list.
# - No bugs, but some interesting recycling of
#   output files from one of Juan's direcoties
#   + This actually was a huge bug based on an
#     inconsistency with batch_dca.py (whether
#     Michael's directory should be used) and
#     led to any interactions that were in mjm
#     but not jfb skipping the coevolution
#     features.

# Imports
import glob, os, sys

# Identify Input Interactions File
inter_file = sys.argv[1]

# Find all of the DCA output files
# Originally Michael's directory was commented out --> some DCA features were
# not being calculated (the DCA files were already in mjm so were not calculated
# fresh, but the parser did not look there so they ended up being skipped)
dca_files = glob.glob('/mnt/data-rs-wihy299nas/mjm659/SCA/dca_july2016_all/*.dca')
dca_files += glob.glob('/mnt/data-rs-wihy299nas/jfb294/dca_aws/*.dca')
dca_files += glob.glob('/home/adr66/eclair/data/sca/dca/*.dca')

# Output directories for DCA MI and DI Scores
mi_dir = '/home/adr66/eclair/data/sca/mi/'
di_dir = '/home/adr66/eclair/data/sca/di/'

# Iterate over each interaction
for l in open(inter_file, 'r'):
	
	# Generate base filename
	# SHOULD BE SORTED OTHERWISE RESULTS IN UNEXPECTED BEHAVIOR?
	basename = '%s_%s' % (l.split('\t')[0], l.split('\t')[1].strip())
	
	# Find original files
	# F0 and F1 are the original DCA files calculated by either Michael or Juan
	# The previous batch_dca.py script first checks these locations to decide if
	# dca needs to be run for this interaction (this is messy caching) so we need
	# to check these locations when we parse as well. We will arbitrarilly prioritize
	# Aaron > Michael > Juan
	f0 = '/mnt/data-rs-wihy299nas/mjm659/SCA/dca_july2016_all/%s.dca' % (basename)
	f1 = '/mnt/data-rs-wihy299nas/jfb294/dca_aws/%s.dca' % (basename)
	f2 = '/home/adr66/eclair/data/sca/dca/%s.dca' % (basename)
	
	# Originally the Michael directory was disabled, to avoid adding a third try / except
	# below I decide which of f0 or f1 should be used here
	if(os.path.exists(f0)):
		f1 = f0
	
	# Create output files
	output_mi = open(os.path.join(mi_dir, basename+'.mi'), 'w')
	output_di = open(os.path.join(di_dir, basename+'.di'), 'w')
	
	# Read in either f1 or f2 and parse it to the split output
	# I'm guessing f1 was just the original DCA files generated
	# before Eclair was better consolidated?
	try:
		# Iterate over each line in the input
		for line in open(f2):
			# Separate line by field
			i, j, mi, di = line.strip().split()
			
			# Write MI to output 1
			output_mi.write('%s\t%s\t%s\n' %(i, j, mi))
			
			# Write DI to output 2
			output_di.write('%s\t%s\t%s\n' %(i, j, di))
		
		# Close outputs
		output_mi.close()
		output_di.close()
	except:
		try:
			# Same as above
			for line in open(f1):
				i, j, mi, di = line.strip().split()
				output_mi.write('%s\t%s\t%s\n' % (i, j, mi))
				output_di.write('%s\t%s\t%s\n' % (i, j, di))
	
			output_mi.close()
			output_di.close()
		except:
			pass
