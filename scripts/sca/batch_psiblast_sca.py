# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script uses Statistical Coupling Analysis to
# calculate the positional correlation matrix for all
# interaction pairs specified in the input file. This
# positional correlation matrix quantitatively indicates
# the correlated evolution for all pairs of positions
# in the joined MSA for the two UniProt IDs in each
# interaction pair. The script skip recalculating this
# feature wherever it has already been calculated.

# Expected Outcomes:
# A joiend MSA and SCA output file should be generated
# corresponding to each of the interaction pairs in
# the input file. The exception to this is any interaciton
# pairs that have a joined MSA alignment length greater
# than 10,000 or any pairs that have fewer than 50
# sequences in their final joined MSA.

# Known Bugs:
# - No bugs, some concerns including the regeneration
#   of Slimmed MSAs generated previously, and some
#   implementation / threshold setting questions.

# Imports
from mjm_tools import TimeoutCmd
from mjm_parsers import parse_fasta, write_fasta
import os, glob, shutil, sys, time
from multiprocessing.pool import ThreadPool

# Log the input interaction (?) and the total run_time for a command
def log_time(interval, e):
	with open('sca_timelog.txt','a') as outf:
		outf.write('%s\t%0.2f\n' % (e, interval))

# Runs a command / Logs its runtime
def do_command(comm):
	(c,e) = comm
	starttime = time.time()
	os.system(c)
	interval = time.time() - starttime
	log_time(interval, e)

# Fixed path to sca
sca_path = '/home/mjm659/ires_ml/data/sca/SCA5_forDist/sca5/'

# Fixed path for input MSAs
uniprot_msa_dir = '../conservation/psiblast/'

# Identifies the interaction file(s)
interactome_files = glob.glob(sys.argv[1])

# Fixed path to pre-computed Interface Residues
# defined per UniProt pair, using a combination
# of PDB sources.
# HEAD:
#
# UniProtA        UniProtB        TaxIDA  TaxIDB  NumIresA        NumIresB        UniProtIresA    UniProtIresB    NumAllUniProtResA     NumAllUniProtResB       UniProtAllResA  UniProtAllResB  PDBSources
# P03502  P05161  518987  9606    30      28      [8,14-16,18-19,29,32-34,36-39,52,75,84,87-88,91-95,97-102]      [8-13,36-37,39,44,46,48-53,72,74-81,99,101]   97      153     [7-103] [3-155] 3R66-A:C;3R66-A:D;3R66-B:C;3R66-B:D;3RT3-C:B;3SDL-A:C;3SDL-A:D;3SDL-B:C;3SDL-B:D
# P55407  P55407  394     394     25      25      [8,12-13,15-16,82-84,121-124,147-149,151-152,155-156,159-160,162-163,171,173] [8,12-13,15-16,82-84,121-124,147-149,151-152,155-156,159-160,162-163,171,173]   235     235     [2-236] [2-236]       2Q0O-A:B
# A0QS40  A0QS41  246196  246196  29      39      [24,27,30-34,36,38,48,50,57-58,61-62,65,81-82,84-93,96] [2-3,20-21,26,29-30,32-40,46,50,59-60,63-65,67-68,71-72,75-78,80-84,86-87,90] 135     141     [7-115,119-144] [2-142] 4RV2-A:B
# E2IPT4  E2IPT4  493803  493803  8       8       [269,273,276,309,311-313,315]   [269,273,276,309,311-313,315]   120  120      [231-350]       [231-350]       3QFQ-B:E
# P45448  Q61066  10090   10090   27      15      [373,380,389-391,394-395,398,440-442,445,483-484,534,537-538,541-544,547-548,550,553,556-557] [270,273-276,278-279,281,283,397-401,408]       242     183     [318-524,526-560]       [251-313,353-472]     3F5C-A:B;3F5C-A:C
#
priority_file = '../pdb/ires_perppi_alltax.txt'


# Output dirs?
msa_output_dir = 'joined_msas'
read_dca_dir = '/mnt/data-rs-wihy299nas/mjm659/SCA/sca_corr_psiblast'
dca_output_dir = 'sca_corr_psiblast'

# Hardcoded number of threads to use?
# BAD
nthreads = 16

# Other fixed features
min_sequence_coverage = 0 #0.25
min_num_sequences = 50 #25

### delete all previous msas ###
#shutil.rmtree(msa_output_dir)
#os.mkdir(msa_output_dir)
###

# Obtain list of all interactions
interactions = set()
for f in interactome_files:
	print 'parsing', f, '...'
	interactions.update(  set([tuple(sorted(l.strip().split('\t')[:2])) for l in open(f)])  )
print len(interactions), 'interactions'

# Mark all interactions included in PDB Ires file
# as "priority interactions"
priority_interactions = set([tuple(sorted(l.split('\t')[:2])) for l in open(priority_file)][1:])

#
# Step 1 - Generate joined MSAs
#

# Find all the aligned MSAs from the psiblast_update.py script
uniprot_MSAs = [os.path.basename(f).replace('_aligned.msa', '') for f in glob.glob(os.path.join(uniprot_msa_dir, '*_aligned.msa'))]

# Keep track of some counters
to_process = [] 
num = 0
heterodimers = 0
het_w_msas = 0
num_priority = 0

# Iterate over each interaction
for i in interactions:
	# Obtain the two UniProt IDs
	p1, p2 = i
	
	# Skip all Homodimers
	# A meaningful SCA cannot be performed with homodimers
	if p1 == p2 : continue
	
	# Update heterodimer counter
	heterodimers += 1
	
	# Skip all interactions that do not both have MSAs
	if p1 not in uniprot_MSAs or p2 not in uniprot_MSAs:
		continue
	
	# Update heterodimer + MSA counter
	het_w_msas += 1
	
	# Read both MSAs as dictionaries
	_, p1_msa = parse_fasta(os.path.join(uniprot_msa_dir, p1+'_aligned.msa'))
	_, p2_msa = parse_fasta(os.path.join(uniprot_msa_dir, p2+'_aligned.msa'))
	
	
	# It looks like we are just recreating the Slimmed MSAs
	# that were generated in the batch_js_conservation.py
	# script. Why don't we just use the Slimmed MSAs here?
	for k in p1_msa.keys():
		p1_msa[k.split('_')[-1]] = p1_msa[k].replace('U', '-')
		del p1_msa[k]
	
	for k in p2_msa.keys():
		p2_msa[k.split('_')[-1]] = p2_msa[k].replace('U', '-')
		del p2_msa[k]
	
	# Skip all pairs with an empty MSA at this step?
	# How is this even possible?
	# Suggest adding a print statement and / or warning
	if len(p1_msa) == 0 or len(p2_msa) == 0:
		continue
	
	# Skip MSAs that tare considered too long
	# and will "take forever to run"
	# The definition of too long here is based on the
	# cummulative length of the alignments and this
	# check is justified since the run_sca should
	# run in O(n**2) time or worse
	p1_len = len(p1_msa[p1_msa.keys()[0]])
	p2_len = len(p2_msa[p2_msa.keys()[0]])
	if p1_len + p2_len > 10000:
		continue
	
	# Begin constructing joined MSA
	joined_msa = {}
	
	# Iterate over all species that are shared between the two
	# MSA. THIS is why we cared about dropping the UniProt IDs
	# and just keeping the species for the key
	for k in set(p1_msa.keys()).intersection(set(p2_msa.keys())):
		
		# Calculate the percentage of each alignment that is gap
		p1_perc_missing = p1_msa[k].count('-') / float(len(p1_msa[k]))
		p2_perc_missing = p2_msa[k].count('-') / float(len(p2_msa[k]))
		
		# Skip entries that have insufficient coverage
		# There is not currently any filtering done
		# based on coverage
		if p1_perc_missing > (1.0-min_sequence_coverage) or p2_perc_missing > (1.0-min_sequence_coverage):
			print '(%s, %s) Not enough sequence coverage for DCA' % (p1, p2)
			continue
		
		# Generate joined MSA by just adding the two
		# individual MSA together
		joined_msa[k] = p1_msa[k] + p2_msa[k]
	
	# Skip entries that have an insufficient number of
	# alignments remaining in the joined MSA
	# DCA does not perform well if there are few
	# sequences in the alignment	
	if len(joined_msa) < min_num_sequences:
		print '(%s, %s) Not enough sequences for DCA' % (p1, p2)
		continue
	
	# Write the joined MSA file
	joined_msa_file = os.path.join(msa_output_dir, '%s_%s.msa' % (p1, p2))
	write_fasta(joined_msa, joined_msa_file)
	
	# Prioritize processing the interactions that have
	# PDB evidence?
	# I think this was probably just something that they
	# wanted to do when they were originally training the
	# model, but not something that needs to be taken
	# into consideration for making predictions
	if (p1, p2) in priority_interactions:
		to_process.insert(0, joined_msa_file)
		num_priority += 1
	else:
		to_process.append(joined_msa_file)
	
	# Update total number of joined MSAs to process
	num += 1

# Status update
print 'Out of %i interactions:\n  %i heterodimers\n  %i with both msas\n  %i passing all filters\n  %i high priority' %(len(interactions), heterodimers, het_w_msas, num, num_priority)


# Figure out which DCA .sca files have been generated previously
finished_dca = set([os.path.basename(f).split('.')[0] for f in glob.glob(os.path.join(dca_output_dir, '*_*.sca'))])

# Figure out which DCA .sca files have been generated previously (on the NAS?)
# I'm guessing maybe this is where they were generated originally, and at some
# point they wanted to consolidate them all in the eclair folder without having
# to regenerate everything? So they just have two caches?
# This actually does not get used at all in the current version of the code?
old_dca = set([os.path.basename(f).split('.')[0] for f in glob.glob(os.path.join(read_dca_dir, '*_*.sca'))])

# Figure out all of the joined MSAs that can be processed
to_process = glob.glob(os.path.join(msa_output_dir, '*_*.msa'))

# Figure out all the joiend MSAs that we actually want to process based on
# the interactions files
want_to_process = set(map(lambda inter: '%s_%s' % (inter[0],inter[1]),interactions))

#
# Step 2 - Generate SCA
#

# Generate a list of commands to run
command_queue = []

# Iterate over all of the joined MSAs
for i, msa in enumerate(to_process):
	print '%s/%s %s' %(i, num, msa)
	
	# Figure out the interaction for this MSA
	interaction = os.path.basename(msa).split('.')[0]
	
	# Skip interactions that we have already processed
	if interaction in finished_dca: # or interaction in old_dca:
		print 'already computed', interaction
		continue
	
	# Skip interactions that are not relevant to the specific
	# interaction file
	if interaction not in want_to_process:
		continue
	
	# Generate output filename
	dca_output = os.path.join(dca_output_dir, '%s.sca' %(interaction))
	
	# Add command to queue
	# The output for the run_sca command should give a positional
	# correlation matrix (Cp) that which quantitatively indicates the
	# correlated evolution of all pairs of positions in the alignment.
	# We can use this on our joined MSAs to figure out if there are
	# possitions accross the interaction pairs that co-evolve.
	#
	# HEAD:
	#
	# i       j       correlation
	# 1       2       2.419753
	# 1       3       2.568013
	# 1       4       2.812190
	# 1       5       3.094885
	# 1       6       2.778559
	#
	command_queue.append(('''matlab -singleCompThread -nodisplay -nosplash -r "addpath('%s');run_sca('%s', '%s');exit"'''  %(sca_path, msa, dca_output), interaction))

# Submit all commands using ThreadPool
my_pool = ThreadPool(nthreads)
command_outputs = my_pool.map(do_command, command_queue)
