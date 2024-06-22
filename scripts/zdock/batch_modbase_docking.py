# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script performs docking simulations on all of
# the docking pairs included in the modbase_docking_set.txt.
# Docking is carried out using mjm_tools wrapped_zdock.py
# script to generate the top 10 docked conformations
# for each docking pair. Wherever possible, the script
# avoid re-doing docking simulations that have been
# completed previously.

# Expected Outcomes:
# Docking results should be generated for each of the
# docking pairs included in the docking set.

# Known Bugs:
# - None

# Imports
import time, subprocess, os, glob, sys
from multiprocessing.pool import ThreadPool
from collections import defaultdict
from mjm_parsers import parse_dictionary_list, unzip_res_range

# Obtain list of interactions
interlist = [sorted(l.strip().split('\t')[:2]) for l in open(sys.argv[1])]

# Obtain number of cores from user provided
# input, otherwise default to 1
if len(sys.argv) > 2:
	num_cores = int(sys.argv[2])
else:
	num_cores = 1

# Log the docking set information and the total run_time for a command
def log_time(interval, e):
	with open('modbase_timelog.txt','a') as outf:
		outf.write('%s\t%0.2f\n' % ('\t'.join([e['ProtA'],e['ProtB'],e['SubA'],e['SubB'],e['CovA'],e['CovB'],e['pCovA'],e['pCovB']]), interval))

# Runs a command / Logs its runtime
def do_command(comm):
	(c, e) = comm
	#~ out, err = subprocess.Popen(c, stdout=subprocess.PIPE).communicate()
	starttime = time.time()
	os.system(c)
	interval = time.time() - starttime
	log_time(interval, e)
	return 'Command complete: %s' %(c)

# Input docking set
index_file = 'modbase_docking_set.txt'

# Fixed Output Directory
output_folder = 'modbase_docked_models'

# Input ModBase Models Directory
modbase_model_dir = '../modbase/models/hash'

# Figure out which pairs of Models have already had
# docking simulations run based on presence of the
# corresponding Zdock output.
already_docked = set([tuple(os.path.basename(f).split('.')[0].split('--')[:2]) for f in glob.glob(os.path.join(output_folder, '*--ZDOCK.out'))])

# Generate Queue of commands
command_queue = []

# Iterate over each row of the docking set
for e in parse_dictionary_list(index_file):
	
	# Set up the dock so that the Model with
	# the higher coverage (presumably the larger
	# model?) is treated as the static receptor
	# and the smaller Model is treated as the
	# ligand to dock.
	if int(e['CovA']) >= int(e['CovB']):
		receptor, ligand = e['SubA'], e['SubB']
	else:
		receptor, ligand = e['SubB'], e['SubA']
	
	# Skip entries where the model pair does not
	# correspond to one of the interactions in the
	# current input file
	if sorted([e['ProtA'],e['ProtB']]) not in interlist:
		continue
	
	# Skip entries where the docking simulations have
	# already been carried out
	if (receptor.upper(), ligand.upper()) in already_docked:
		print 'ALREADY DOCKED:', receptor, ligand
		continue
	
	# Obtain path to the input receptor and ligand
	receptor_file = os.path.join(modbase_model_dir, receptor+'.pdb')
	ligand_file = os.path.join(modbase_model_dir, ligand+'.pdb')
	# The chains are null since these are ModBase models
	# that only have one chain
	receptor_chain, ligand_chain = '_', '_'
	
	# Add Zdock command to Queue
	# Use mjm_tools wrapped_zdock.py to dock the
	# ligand to the receptor, selecting the specified
	# chain in each. Keep the receptor structure fixed.
	# Fix mod_base PDB files (fields 73-80?). Output
	# results to the output_folder.
	#
	# The output for the wrapped_zdock.py includes 10
	# PDB files corresponding to the 10 best scoring
	# docks, and a .out file, which I have no idea
	# how to interpret currently.
	command_queue.append(('wrapped_zdock.py %s %s -RC %s -LC %s -F --fix_modbase --outdir %s' %(receptor_file, ligand_file, receptor_chain, ligand_chain, output_folder), e))

# Submit all of the docking jobs to a threadpool
my_pool = ThreadPool(num_cores)
command_outputs = my_pool.map(do_command, command_queue)
