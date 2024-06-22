# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script continues parsing the previously downloaded
# raw ModBase XML data obtained through modbase_DL.py.
# This script continues by iterating through all of the
# valid models obtained through  modbase_02.py to ensure
# that all of the models have correct indexing that matches
# up 1:1 with the indexing of the UniProt Sequence.

# Expected Outcomes:
# All of the previously generated model PDB files from
# modbase_02.py should be updated to include the proper
# indexing.

# Known Bugs:
# - Strongly suspect that the correctly index may not
#   be calculated correctly, since it determines its
#   offset based on the current index alone and does not
#   take a difference between the current index and the
#   original first index. I have not been able to confirm
#   that this is the case, and must be misinterpretting
#   something because this could not have gone uncaught.
# - Don't know how we can be certain that all of the indices
#   will match up just from having looked at the first ATOM.
# - What happens to the models with "bad_indices"?
# - Should be throwing warnings.

# File Headers (I think this is because this
# script was copied out of the nightly_scripts
# to create a stand alone Eclair ModBase library
#
#!/usr/bin/env python
#DEPENDENCIES: modbase_02.py
#
#EMAIL_LOG:
#EMAIL_ERR: mjm659@cornell.edu

'''Fixes index problems in all hash models'''

# Imports
import glob, os, sys

# Fixed input dir for models
hashDir = 'models/hash'

# Keep track of some counters
num_fixed = 0
bad_indices = 0
model_files = sorted(glob.glob(os.path.join(hashDir, '*.pdb')))

# Iterate through every model to correct
# any index problems
for model in model_files:
	# Open File
	file_handle = open(model, 'r')
	
	# Keep track of where the target UniProt is supposed
	# to start
	model_start = 1
	
	# Flag for special case of first AA
	first_aa = True
	
	# Save corrected lines
	fixed_lines = []
	
	# Iterate through every line in file
	for line in file_handle:
		
		# Save the starting position in the UniProt target
		if 'REMARK 220 TARGET BEGIN:' in line:
			model_start = int(line.replace('REMARK 220 TARGET BEGIN:', '').strip())
		
		# For all ATOM lines
		if line[:4] == 'ATOM':	
			# Split the line in half, first half - atom
			# identifier information, second half - atom
			# position information
			#
			# Note: ATOM entries in PDB are are formatted
			#       as follows...	
			#       
			#       https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
			#
			#			COLUMNS        DATA  TYPE    FIELD        DEFINITION
			#			-------------------------------------------------------------------------------------
			#			1 -  6        Record name   "ATOM  "
			#			7 - 11        Integer       serial       Atom  serial number.
			#			13 - 16        Atom          name         Atom name.
			#			17             Character     altLoc       Alternate location indicator.
			#			18 - 20        Residue name  resName      Residue name.
			#			22             Character     chainID      Chain identifier.
			#			23 - 26        Integer       resSeq       Residue sequence number.
			#			27             AChar         iCode        Code for insertion of residues.
			#			31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
			#			39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
			#			47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
			#			55 - 60        Real(6.2)     occupancy    Occupancy.
			#			61 - 66        Real(6.2)     tempFactor   Temperature  factor.
			#			77 - 78        LString(2)    element      Element symbol, right-justified.
			#			79 - 80        LString(2)    charge       Charge  on the atom.
			#
			first_line_half = line[:26]
			second_line_half = line[26:]
			
			# Get residue sequence number for current atom
			current_index = first_line_half.split()[-1]
			
			# For first AA, save the original index
			if first_aa:
				try:
					original_first_index = int(current_index)
				# If index is poorly formatted, skip this model
				# Should also throw warning?
				# How do we know which ones this happened to?
				except:
					#print model, current_index
					bad_indices += 1
					break
			
			# If it looks like the indices are going to match up
			# then we can skip fixing this model.
			# How do we know there won't be any gaps or jumps
			# in the indexing later on that will throw it off?
			if first_aa and int(current_index) == model_start:
				break
			
			# If line needs to be fixed, trim off the original index
			first_line_half = first_line_half[:-1*len(current_index)] + ' '*len(current_index)
			
			# Calculate the new index
			# Based on the way this is being calculated, it looks like
			# the assumption is that the current_index must always
			# start at 1? But I don't believe that that is correct
			# because the script saves the original_first_index, as
			# if we should have to do some offset calculation using it.
			# I suspect there may be an error here?
			new_index = str(int(current_index) + model_start - 1)
			
			# Add the corrected index
			first_line_half = first_line_half[:-1*len(new_index)] + new_index
			
			# Save this corrected line
			fixed_lines.append(first_line_half + second_line_half)
			
			first_aa = False
		# For other lines, no corrections necessary
		else:
			fixed_lines.append(line)
	
	# If no corrections were necessary, close
	# the files and continue
	if first_aa:
		file_handle.close()
		continue	
	
	# Close the file
	file_handle.close()
	
	# Write corrected file
	print 'fixing:', model, original_first_index, 'to', model_start
	num_fixed += 1
	#print model
	output_file = open(model, 'w')
	output_file.writelines(fixed_lines)
	output_file.close()

# Final Status Update
print '-----------------------------'
print '%s models had non-int indices' %(bad_indices)
print 'fixed %s out of %s total files' %(num_fixed, len(model_files))
