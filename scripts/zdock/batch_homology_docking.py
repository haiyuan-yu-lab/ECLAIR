
import time, subprocess, os, glob, sys
from multiprocessing.pool import ThreadPool
from collections import defaultdict
from mjm_parsers import parse_dictionary_list, unzip_res_range

interlist = [sorted(l.strip().split('\t')[:2]) for l in open(sys.argv[1])]

if len(sys.argv) > 2:
	num_cores = int(sys.argv[2])
else:
	num_cores = 1

def do_command(c):
	os.system(c)

index_file = 'homology_docking_set.txt'
output_folder = 'pdb_docked_models'

already_docked = set([tuple(os.path.basename(f).split('.')[0].split('--')[:2]) for f in glob.glob(os.path.join(output_folder, '*--ZDOCK.out'))])

#------------------------------------

command_queue = []
for e in parse_dictionary_list(index_file):
	
	if int(e['CovA']) >= int(e['CovB']):
		receptor, ligand = e['SubA'], e['SubB']
	else:
		receptor, ligand = e['SubB'], e['SubA']
		
	#if sorted([e['ProtA'],e['ProtB']]) not in interlist:
	#	continue

	if (receptor, ligand) in already_docked:
		print 'ALREADY DOCKED:', receptor, ligand
		#continue
	
	receptor_pdb, receptor_chain = receptor.split('_')
	ligand_pdb, ligand_chain = ligand.split('_')
	
	#~ print 'DOCKING:', receptor, ligand
	command_queue.append('wrapped_zdock.py %s %s -RC %s -LC %s -N 2000 -F --outdir %s' %(receptor_pdb, ligand_pdb, receptor_chain, ligand_chain, output_folder))
	

my_pool = ThreadPool(num_cores)
command_outputs = my_pool.map(do_command, command_queue)
