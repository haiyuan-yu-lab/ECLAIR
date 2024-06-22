# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script performs Direct Couplind Analysis (DCA) in
# order to determin the Mutual and Direct information
# between any two positions in a MSA. For our purposes
# we can use this to estimate the direct evolutionary
# pressure acting jointly on two residues between a
# pair of interacting proteins. This script skips
# recalculating these features wherever it has already
# done so.

# Expected Outcomes:
# The DCA output should be calculated for every pair
# of interactions in the input file. I believe the
# restriction that the alignment must not exceed
# 10000 residues in length and that there must be at
# least 50 sequences in the alignment indirectly apply
# but they are not hardcoded here.

# Known Bugs:
# - The timeout on the DCA threads sometimes fails and
#   the process must be killed manually
# - Problem with the way the number of threads to be used
#   is hardcoded + general difference in opinions of the
#   way the multithreading is handled
# - Confusion on whether or not the filters from the
#   psiblas_sca.py script indirectly carry over or
#   if they 

# Imports
from mjm_tools import TimeoutCmd
from mjm_parsers import parse_fasta, write_fasta
import os, glob, shutil, time, sys, threading, subprocess
from multiprocessing.pool import ThreadPool
from Queue import Queue

# Limit maximum running time for a given DCA process
# to 30 minutes. My understanding is that DCA can
# run indefinitely and needs to be manually stopped
# sometimes. Further, sometimes enforcing the timeout
# does not work.
MAX_TIME_DCA = 1800

# Logs the amount of time that a given interaction
# runs for
def log_time(interval, inter):
	with open('dca_timelog.txt','a') as outf:
		outf.write('%s\t%0.2f\n' % (inter, interval))

# Runs dca on arguments in the cmd_queue, stops when it gets a sentinel
def dca(cmd_queue, sentinel):
	while True:
		c = cmd_queue.get()
		if c is sentinel:
			break
		
		# Get the command and interaction
		cmd = c[0]
		interaction = c[1]
		
		# Get start time
		starttime = time.time()
		print('Doing command for interaction: %s' % (interaction))
		
		# Start executing command in a Thread
		thread = SingleCommandThread(cmd)
		thread.start()
		
		# Timeout the thread after max time
		thread.join(MAX_TIME_DCA)
		
		# If the thread has not finished executing
		# print an error / kill the thread
		if (thread.is_alive()):
			print('Terminating for interaction: %s' % (interaction))
			sys.stderr.write("Terminating dca for interaction: %s - timeout for max_time %s \n" % (interaction, MAX_TIME_DCA))
			thread.p.terminate()
		# Otherwise, DCA calculation is successfull
		else:
			print('Succeeded for interaction: %s' % (interaction))
		cmd_queue.task_done()
		
		# Log the command
		log_time(time.time() - starttime, interaction)

# Thread class that runs and waits on a single command when started
class SingleCommandThread(threading.Thread):
	def __init__(self, cmd):
		super(SingleCommandThread, self).__init__()
		self.cmd = cmd
        
	def run(self):
		self.p = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		self.p.wait()


#priority_interactions = set(['_'.join(sorted(l.strip().split('\t')[:2])) for l in open('../../features/per_feature/ires50.txt')])

# Identify interaction input files
interactome_files = glob.glob(sys.argv[1])

# Find all of the joined MSAs
msas = glob.glob('/home/adr66/eclair/data/sca/joined_msas/*.msa')

# Fixed directory paths
# EDIT - SHAYNE 2019_12_04 Disable usage of old caches
#read_dir = '/mnt/data-rs-wihy299nas/mjm659/SCA/dca_july2016_all/'
#read_dir2= '/mnt/data-rs-wihy299nas/jfb294/dca_aws/'
output_dir = '/mnt/data-rs-wihy299nas/adr66/data/sca/dca/'
read_dir3 = output_dir
dca_path = '/home/adr66/eclair/data/sca/'

# Hardcoded number of threads to use
# BAD
num_threads = 10
sentinel = None

# Figure out all of the interactions in the input files
interactions = set()
for f in interactome_files:
	interactions.update(set(['_'.join(sorted(l.strip().split('\t')[:2])) for l in open(f)]))

# Sort the MSAs by length of alignment so that
# we start with the fastest running DCA calculations
msas.sort(key=lambda f: len(open(f).readlines()[1]))

# Obtain a set of all of the completed DCA results
# I am guessing that all of the read_dirx directories
# are direcoties where some of these results were saved
# before everything was consolidated in one place
finished = set(glob.glob(os.path.join(output_dir, '*.dca')))
# EDIT - SHAYNE 2019_12_04 Disable usage of old caches
#finished.update(set(glob.glob(os.path.join(read_dir, '*.dca'))))
#finished.update(set(glob.glob(os.path.join(read_dir2,'*.dca'))))
#finished.update(set(glob.glob(os.path.join(read_dir3,'*.dca'))))

# Figure out which interactions have already had this feature calculated
already_calculated = set(map(lambda f: os.path.basename(f)[:-4],finished))

# Generate a Queue for running commands in
command_queue = Queue()

print("Starting to look through msa files")
# Iterate through all MSA files
for i, msa in enumerate(msas):
	# Obtain interaction (P1_P2)
	interaction = os.path.basename(msa)[:-4]
	
	# Generate the outpue filename
	output_file = os.path.join(output_dir, '%s.dca' %(interaction))
    
	#~ if interaction not in priority_interactions:
		#~ continue
	
	# Skip interactions that have already had their DCA calculated
	if interaction in already_calculated:
		continue
	
	# Skip interactions that are not in the interaction lists
	if interaction not in interactions:
		continue
	
	
	print("Running dca for interaction: %s" % (interaction))
	
	# Add the command to the Queue
	# The output for the lowram_dca command should give a matrix
	# containing analogous to that provided through the SCA
	# outputs. These outputs provide two values for the relationship
	# between all pairs of positions in the MSA: 1) MI, the mutual
	# information between i and j, 2) DI, the direct information
	# between i and j. These results can be parsed to obtain
	# estimates of the direct evolutionary pressure between
	# any two positions between the interacting proteins. Note,
	# the format is actually space delimited rather than tab
	# delimited as shown below.
	#
	# HEAD:
	#
	# i	j	MI			DI
	# 1	2	1.41257	0.0350167
	# 1	3	0.937296	0.00696671
	# 1	4	1.23918	0.0073415
	# 1	5	0.907056	0.0105828
	# 1	6	1.1267	0.00986768
	#
	cmd = ["matlab", "-singleCompThread", "-nodisplay", "-nosplash", "-r", "addpath('%s');lowram_dca('%s', '%s');exit" % (dca_path, msa, output_file)]
	command_queue.put((cmd, interaction))

# Generate / start 10 threads to handle executing DCA commands
# Each thread will fetch the next interaction to run from
# a common pool, so this should be a safe way to split up
# the inputs
# Note: Since this is hard coded at 10 instead of num_threads
#       this script will break if num_threads is ever set
#       to anything other than 10
threads = [threading.Thread(target=dca, args=(command_queue, sentinel)) for n in range(10)]
for t in threads:
	t.start()

# Wait for all the inputs to be removed from the queue, then add
# the sentinel to kill all of the threads
# I don't like this approach at all
command_queue.join()
for i in range(num_threads):
	command_queue.put(sentinel)

# Join all of the threads to "double check that they
# were killed"
# I'm not sure if this actually works, because if for
# some reason they were not killed this would just
# hang forever?
for t in threads:
	t.join()
