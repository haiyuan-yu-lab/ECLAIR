# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script begins to process the raw ModBase
# XML downloads obtained in modbase_DL.py. This
# script performs a first pass through the data
# with the main purpose of identifying all of the
# models that are proposed to represent each
# UniProt ID, filtering out models that either
# are improperly formatted, or have information
# that is not consistent with our own query from
# the UniProt API (e.g. length). All invalid
# models are removed while the valid models are
# written to individual PDB files and have their
# header information saved to text files.

# Expected Outcomes:
# All UniProt IDs that were included in the previous
# raw ModBase download should have their models
# parsed, and PDB files and HEADER text should be
# saved for all valid models.

# Known Bugs:
# - There should be warnings / optional prints
#   for all of the cases where a model may be
#   deemed invalid
# - There is not currently any check to delete
#   previously downloaded models that are no
#   longer valid (maybe this is yet to come in
#   a further modbase step?)

# File Headers (I think this is because this
# script was copied out of the nightly_scripts
# to create a stand alone Eclair ModBase library
#
#!/usr/bin/env python
#DEPENDENCIES: modbase_DL.py
#
#EMAIL_LOG:
#EMAIL_ERR: mjm659@cornell.edu

'''Keeps data/hash and data/headers directories updated with latest models from the data/uniprot models.'''

# Imports
import os, glob, sys, re
from os import path, popen, system
from subprocess import call
from mjm_parsers import parse_dictionary_list

# Set up direcoties paths / UniProt Information
# File based on input parameters, otherwise
# use defaults.
if len(sys.argv) > 1:
	uniprotDir = sys.argv[1]
	hashDir = sys.argv[2]
	headersDir = sys.argv[3]
	uniprotInfoFile = sys.argv[4]
else:
	uniprotDir = 'models/uniprot'
	hashDir = 'models/hash'
	headersDir = 'models/headers'
	uniprotInfoFile = '../uniprot/uniprot_info.txt'


# Read UniProt Info File as dictionary each UniProt ID
# maps to a dictionary of its attributes.
uniprot2info = dict((e['id'], e) for e in parse_dictionary_list(uniprotInfoFile))

# Save a dictionary that will map the unique model
# identifiers for each of the ModBase models to
# the UniProt ID that they actually represent.
hash_id2uniprot = {}

# Maintain a set of all of the unique ModBase
# models that have already been downloaded / parsed
# out. When we compare the existing set of downloaded
# hash ids to the updated raw ModBase download files
# and we encounter duplicates we will compare them
# to decide which to delete and which to parse.
downloaded_hash_ids = set([path.splitext(os.path.basename(p))[0] for p in glob.glob(path.join(hashDir, '*.pdb'))])

# All of the ModBase raw files downloaded for each UniProt
uniprot_files = sorted(glob.glob(path.join(uniprotDir, '*.pdb')))

# Keep track of some counters
total_models = 0
valid_models = set()
valid_uniprots = set()

# Iterate over each UniProt ModBad file.
# This is an original base to determin which model / hash ids
# are contained in each UniProt file.
for i, f in enumerate(uniprot_files):
	# Progress Message
	if i%1000==0:
		print i, f
	
	# Obtain UniProt ID
	uniprot = path.splitext(os.path.basename(f))[0]
	
	# Skip if UniProt ID not in UniProt info file
	# Not sure how this would happen?
	if uniprot not in uniprot2info:
		continue
	
	# The downloaded ModBase files are originally
	# in XML formating. This line just parses
	# the original file to pull out all of the
	# actual PDB content. These original files can
	# contain several ModBad models within the same
	# file under different <pdbfile> tags.
	#
	# HEAD:
	# <?xml version="1.0" encoding="ISO-8859-1"?>
	# <files>
	# <pdbfile>
	# 		<model_id>08928893801e4b64765dc2680b54c9a1</model_id>
	# 		<content>
	# HEADER    ModPipe Model of 57162264 (X6RE90)      2005-06-0
	# TITLE     Model of
	# SOURCE    Homo sapiens
	# AUTHOR     URSULA PIEPER, BENJAMIN WEBB, EASHWAR NARAYANAN, ANDREJ SALI
	# REMARK 220 EXPERIMENTAL DETAILS
	# ...
	# ATOM      1  N   LYS   120      25.875  -7.545  26.901  1.00 90.23       1SG   2
	# ATOM      2  CA  LYS   120      26.165  -8.891  26.359  1.00 90.23       1SG   3
	# ...
	# ATOM    578  O   ASN   194      44.885  -4.279   1.687  1.00 73.28           O
	# ATOM    579  OXT ASN   194      43.158  -2.990   1.133  1.00 73.28           O
	# TER     580      ASN   194
	# END
	# 		</content>
	# </pdbfile>
	# <pdbfile>
	# ...
	# </pdbfile>
	# </files>
	#
	models = [m.replace('<content>', '').replace('</content>', '').strip() for m in re.findall(r'<content>.*?</content>', open(f).read(), flags=re.DOTALL)]
	
	# Iterate over each model
	for model in models:
		
		# Update model counter
		total_models += 1
		
		#
		# Sanity check to make sure ModBase target length matches UniProt length
		#
		
		# Note: Target length is being used here in the sense
		# of the length for a BLAST result (PDB template was the
		# query and the UniProt ID was the reference / target).
		# This is NOT the length of the actual model.
		target_length = [d.replace('TARGET LENGTH:', '').strip() for d in re.findall(r'TARGET LENGTH:.*', model)]
		
		# Skip entries that have more than one length remark
		# Throw warning?
		if len(target_length) != 1:
			continue
		
		# Skip entries with invalid target length formating
		# Throw warning?
		try:
			target_length = int(target_length[0])
		except:
			continue
		
		# Skip entries whose target length does not match
		# the expected UniProt length. If these values ever
		# do not match up, it likely indicates that the two
		# sources are relying on different versions of the
		# UniProt entry and hence the results cannot be matched
		# up together.
		# Throw warning
		if target_length != int(uniprot2info[uniprot]['length']):
			continue
		
		#
		# Parse Other Fields in Model
		#
		
		# Obtain Model ID
		# Skip Entries with multiple remarks
		hash_id = [d.replace('MODPIPE MODEL ID:', '').strip() for d in re.findall(r'MODPIPE MODEL ID:.*', model)]
		if len(hash_id) > 1:
			continue
		hash_id = hash_id[0]
		
		# Obtain Sequence Identity (i.e. BLAST pident)
		# Skip Entries with multiple remarks
		sequence_identity = [d.replace('SEQUENCE IDENTITY:', '').strip() for d in re.findall(r'SEQUENCE IDENTITY:.*', model)]
		if len(sequence_identity) != 1:
			continue
		
		# Obtain Model Score
		# Skip Entries with multiple remarks
		model_score = [d.replace('MODEL SCORE:', '').strip() for d in re.findall(r'MODEL SCORE:.*', model)]
		#print model_score
		if len(model_score) != 1:
			continue
		
		# Obtain Quality Score
		# Skip Entries with multiple remarks
		MPQS = [d.replace('ModPipe Quality Score:', '').strip() for d in re.findall(r'ModPipe Quality Score:.*', model)]
		if len(MPQS) != 1:
			continue
		
		# Obtain zDOPE Score
		# Skip Entries with multiple remarks
		zDOPE = [d.replace('zDOPE:', '').strip() for d in re.findall(r'zDOPE:.*', model)]
		if len(zDOPE) != 1:
			continue
		
		# Obtain Evalue
		# Skip Entries with multiple remarks
		evalue = [d.replace('EVALUE:', '').strip() for d in re.findall(r'EVALUE:.*', model)]
		if len(evalue) != 1:
			continue
		
		# Obtain Template
		# Skip Entries with multiple remarks
		pdb = [d.replace('TEMPLATE PDB:', '').strip() for d in re.findall(r'TEMPLATE PDB:.*', model)]
		if len(pdb) != 1:
			continue
		
		# Obtain Target Start (i.e. region covered
		# by model)
		# Skip Entries with multiple remarks
		target_begin = [d.replace('TARGET BEGIN:', '').strip() for d in re.findall(r'TARGET BEGIN:.*', model)]
		if len(target_begin) != 1:
			continue
		
		# Obtain Target End (i.e. region covered
		# by model)
		# Skip Entries with multiple remarks
		target_end = [d.replace('TARGET END:', '').strip() for d in re.findall(r'TARGET END:.*', model)]
		if len(target_end) != 1:
			continue
		
		# Obtain Template Start (i.e. part of template
		# that is homologous)
		# Skip Entries with multiple remarks
		template_begin = [d.replace('TEMPLATE BEGIN:', '').strip() for d in re.findall(r'TEMPLATE BEGIN:.*', model)]
		if len(template_begin) != 1:
			continue
		
		# Obtain Template End (i.e. part of template
		# that is homologous)
		# Skip Entries with multiple remarks
		template_end = [d.replace('TEMPLATE END:', '').strip() for d in re.findall(r'TEMPLATE END:.*', model)]
		if len(template_end) != 1:
			continue
		#-----------------------------------------
		
		# Calculate template length
		# Skip Entries with poorly formated
		# target start / end
		try:
			template_length = int(target_end[0]) -  int(target_begin[0])
		except:
			continue
		
		# Update that this is a model that can be
		# parsed accurately / that this UniProt
		# ID has valid models
		valid_models.add(hash_id)
		valid_uniprots.add(uniprot)
		
		# Write all of this header information
		# to a header file for the model
		header = [uniprot, template_length, target_length, pdb[0], target_begin[0], target_end[0], sequence_identity[0], model_score[0], MPQS[0], zDOPE[0], evalue[0], hash_id]
		header_handle = open(path.join(headersDir, hash_id+'.txt'), 'w')
		header_handle.write('\t'.join([str(i) for i in header]))
		header_handle.close()
		
		
		# Write the PDB contents of each model
		# to a model file
		model_handle = open(path.join(hashDir, hash_id+'.pdb'), 'w')
		model_handle.write(model)
		model_handle.close()

print '%s hash_ids already parsed in ../hash dir\n' %(len(downloaded_hash_ids))
print '%s total hash_ids in <uniprot>.pdb modbase dump files' %(total_models)
print '   %s valid models (%s uniprots)\n' %(len(valid_models), len(valid_uniprots))


# Remove any previously downloaded models
# that are no longer valid on this most
# recent pass through the data.
#
# Why is this commented out?
#
# I guess it must have something to do with not
# wanting to delete all of the models for every
# UniProt ID not included in the current set of
# interactions, but this still seems like a bad
# idea to have no check.
#~ to_remove = downloaded_hash_ids - valid_models
#~ 
#~ print 'Removing %s hash models...' %(len(to_remove))
#~ 
#~ for hash_id in to_remove:
	#~ #print 'removing %s' %(hash_id)
	#~ call(['rm', path.join(hashDir, hash_id+'.pdb')])
	#~ call(['rm', path.join(headersDir, hash_id+'.txt')])
