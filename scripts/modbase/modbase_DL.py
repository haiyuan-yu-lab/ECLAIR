# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script fetches ModBase information for all
# UniProt IDs in a given input file of interactions
# and saves the information in a specified output
# file. Wherever possible, the script identifies
# UniProt ID information that has previously been
# fetched, so as to avoid fetching the same information
# multiple times.

# Expected Outcomes:
# All UniProt IDs in the input file should have their
# information fetched, and saved to a UniProt specific
# output file wherever it is available.

# Known Bugs:
# - None

# Imports
from subprocess import call
import glob, os, sys
from time import sleep

# Set up input file /output directory
# based on input parameters
# If no input is provided, default to fixed
# input file / output directory
if len(sys.argv) > 1:
	input_file = sys.argv[1]
	outputDir = sys.argv[2]
else:
	input_file = '../uniprot/all_uniprots.txt'
	outputDir = '../models/'

# Obtain set of all UniProt IDs
all_uniprots = set([line.strip().split('\t')[0] for line in open(input_file)])

# Obtain set of all UniProt IDs that have previously
# had their ModBase data downloaded
already_downloaded = set([os.path.basename(f).split('.')[0] for f in glob.glob(os.path.join(outputDir, '*.pdb'))])

# Obtain set of UniProt IDs that still need
# to be downloaded
to_download = all_uniprots - already_downloaded

print 'downloading models for %s new uniprots.' %(len(to_download))

# Iterate over each UniProt ID
for uniprot in to_download:
	
	# Fetch ModBase data for the UniProt
	# These files are NOT actually PDB files.
	# The files are downloaded in an XML tag
	# format and each file can contain the contents
	# of several PDB files under different
	# <pdbfile> tags. These files still need to be
	# parsed to separate out each of the ModBase
	# models for each UniProt into individual files.
	#
	# HEAD:
	# <?xml version="1.0" encoding="ISO-8859-1"?>
	# <files>
	# <pdbfile>
	#     <model_id>08928893801e4b64765dc2680b54c9a1</model_id>
	#     <content>
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
	#     </content>
	# </pdbfile>
	# <pdbfile>
	# ...
	# </pdbfile>
	# </files>
	#
	command = ["wget", "--quiet", "http://salilab.org/modbase/retrieve/modbase/?databaseID="+uniprot, "-O", os.path.join(outputDir, uniprot+".pdb")]
	print ' '.join(command)
	call(command)
	
	# Sleep? Presummably to avoid timeout errors
	sleep(4)

print 'done.'
