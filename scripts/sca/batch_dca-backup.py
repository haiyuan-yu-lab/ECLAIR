from mjm_tools import TimeoutCmd
from mjm_parsers import parse_fasta, write_fasta
import os, glob, shutil, time, sys
from multiprocessing.pool import ThreadPool

def log_time(interval, inter):
	with open('dca_timelog.txt','a') as outf:
		outf.write('%s\t%0.2f\n' % (inter, interval))

def do_command((c, interaction)):
	starttime = time.time()
	#os.system(c)
	c.Run()
	interval = time.time() - starttime
	log_time(interval, interaction)


#priority_interactions = set(['_'.join(sorted(l.strip().split('\t')[:2])) for l in open('../../features/per_feature/ires50.txt')])
interactome_files = glob.glob(sys.argv[1])

msas = glob.glob('/home/adr66/eclair/data/sca/joined_msas/*.msa')
read_dir = '/mnt/data-rs-wihy299nas/mjm659/SCA/dca_july2016_all/'
read_dir2= '/mnt/data-rs-wihy299nas/jfb294/dca_aws/'
output_dir = '/mnt/data-rs-wihy299nas/adr66/data/sca/dca/'
read_dir3 = output_dir
dca_path = '/home/adr66/eclair/data/sca/'

nthreads = 10

interactions = set()
for f in interactome_files:
	interactions.update(set(['_'.join(sorted(l.strip().split('\t')[:2])) for l in open(f)]))

msas.sort(key=lambda f: len(open(f).readlines()[1]))   #sort by total length of alignment, low to high

finished = set(glob.glob(os.path.join(output_dir, '*.dca')))
finished.update(set(glob.glob(os.path.join(read_dir, '*.dca'))))
finished.update(set(glob.glob(os.path.join(read_dir2,'*.dca'))))
finished.update(set(glob.glob(os.path.join(read_dir3,'*.dca'))))
already_calculated = set(map(lambda f: os.path.basename(f)[:-4],finished))
command_queue = []

# Testing
msas = ["/home/adr66/eclair/data/sca/joined_msas/O43236_Q8IYM1.msa"]

for i, msa in enumerate(msas):
	interaction = os.path.basename(msa)[:-4]  #P1_P2
	output_file = os.path.join(output_dir, '%s.dca' %(interaction))
    
	#~ if interaction not in priority_interactions:
		#~ continue
	
	if interaction in already_calculated:
		continue
	if interaction not in interactions:
		continue
	timeout = 1800
	myCmd = TimeoutCmd(['matlab','-singleCompThread','-nodisplay','-nosplash','-r','''"addpath('%s');lowram_dca('%s', '%s');exit"''' % (dca_path, msa, output_file)], timeout)
	command_queue.append((myCmd, interaction))
	#command_queue.append(('''matlab -singleCompThread -nodisplay -nosplash -r "addpath('%s');lowram_dca('%s', '%s');exit"'''  %(dca_path, msa, output_file), interaction))

my_pool = ThreadPool(nthreads)
command_outputs = my_pool.map(do_command, command_queue)
