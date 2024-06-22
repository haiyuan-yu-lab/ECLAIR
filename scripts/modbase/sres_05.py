# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script calculates the Solvent Accessible Surface
# Area (SASA) for all of the models remaining in the
# select models summary file. Despite being called a
# surface residue calculating script / calling sres_calc.py
# I believe this is just a fancy wrapper for the Naccess
# program that we use to calculate SASA.

# Expected Outcomes:
# All models included in the select models summary file
# should have their SASA calculated / added to the output
# file.

# Known Bugs:
# - Add warnings

'''Batch calculate all native solvent-accessible surface area (SASA) for residues in all ModBase Models'''

# Imports
import sys, os, time
from mjm_parsers import unzip_res_range, parse_dictionary_list
from mjm_tools import fetch_uniprot, xtool, time_elapsed
start_time = time.time()

# Obtian UniProt info file from input
# parameters otherwise, use default
if len(sys.argv) > 1:
	uniprot_info_file = sys.argv[1]
else:
	uniprot_info_file = '../uniprot/uniprot_info.txt'

# Fixed select models input
model_index_file = 'models/parsed/select_modbase_models.txt'

# Fixed SASA output
output_file = 'SASA_modbase.txt'

# Fixed model PDB input dir
model_file_dir = 'models/hash'

# Read UniProt Info as dictionary
# UniProt ID --> dictionary of attributes
uniprot2info = dict((e['id'], e) for e in parse_dictionary_list(uniprot_info_file))

# Read select models summary as dictionary
# List of dictionaries
model_index = parse_dictionary_list(model_index_file)

# Obtain set of all Models that have previously
# had SASA calculation performed
# This is a dictionary of model : output lines
previously_calculated = {}
if os.path.exists(output_file):
	for l in open(output_file):
		previously_calculated[l.split('\t')[11]] = l


# Open output file for writting
# Why isn't this append?
output = open(output_file, 'w')

# List of columns
headers = ['uniprot', 'template_length', 'target_length', 'template_pdb', 'target_begin', 'target_end', 
			'sequence_identity', 'model_score', 'modpipe_quality_score', 'zDOPE', 'eVALUE', 'modbase_modelID']

# Write header (add SASA as an extra final output column)
output.write('\t'.join(headers+['SASA'])+'\n')

# Iterate over ever select model
for e in model_index:
	
	# Skip models that were previously calculated
	if e['modbase_modelID'] in previously_calculated:
		# Just re-write the old output
		output.write(previously_calculated[e['modbase_modelID']])
		continue
	
	# If the model's UniProt is not in the info file
	# Skip it. (How is this possible?)
	# Throw warnings?
	if e['uniprot'] not in uniprot2info:
		continue
	
	# Try to calculate Surface residues
	try:
		
		# Get the appropriate PDB file for this model
		model_file = os.path.join(model_file_dir, e['modbase_modelID']+'.pdb')
		
		# Call srescalc.py on the model the net effect
		# is that this returns a list that looks like...
		#
		# [ [Unknown, Resi, Resn, SASA] ... ]
		#
		# The raw output from srescalc.py is as follows...
		#
		# HEAD:
		#
		#	Unknown	Resi	Resn	SASA
		# 				116	THR	97.3
		# 				117	ALA	78.3
		# 				118	SER	70.6
		# 				119	VAL	45.6
		# 				120	ALA	49.4
		#
		out = [x.replace('\n', '').split('\t') for x in xtool('srescalc.py', model_file, uSASA=-1, o='residue_stats').split('\n')][:-1]
		
		# Generate dictionary mapping UniProt residue to SASA value
		uniprotSASA = dict([(int(q[1]), q[3]) for q in out])
		
		# Generate SASA string for the model
		# This is a semi-colon delimited list of SASA values
		# for ALL positions in the UniProt. Wherever this
		# particular model does not cover a given residue
		# its SASA is listed as NaN
		SASAs = [str(uniprotSASA[r]) if r in uniprotSASA else 'NaN' for r in range(1, int(uniprot2info[e['uniprot']]['length'])+1)]
		
		# Write the output
		output.write('\t'.join([e[header] for header in headers]+ [';'.join(SASAs)])+'\n')
	
	# If something goes wrong just skip it
	# Warnings?
	except:
		continue

# Close the output
output.close()

print 'Finished calculating SASA values for %s ModBase models. Time elapsed: %s' %(len(model_index), time_elapsed(start_time))
