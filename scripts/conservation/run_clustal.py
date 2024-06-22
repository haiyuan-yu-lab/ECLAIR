import subprocess, argparse, sys, os, glob
from helper import *
from multiprocessing.pool import ThreadPool
from mjm_parsers import parse_fasta, parse_dictionary_list

def restricted_float(x):
	x = float(x)
	if x < 0.0 or x > 1.0:
		raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % x)
	return x

def num_cores(x):
	x = int(x)
	if x < 1 or x > 144:
		raise argparse.ArgumentTypeError('%d not in range [1, 144]' % x)
	return x

parser = argparse.ArgumentParser(description='Find and align similar sequences to a protein using PSIBLAST and Clustal Omega.')
parser.add_argument('-c', '--cutoff', help='Cutoff e-value for all PSIBLAST matches used for alignment', type=float, default=0.05, required=False)
parser.add_argument('-r', '--range', help='Fraction of the query sequence required to consider PSIBLAST match for alignment', type=restricted_float, default=0.5)
parser.add_argument('--mult', help='Flag to include multiple proteins for single species', action='store_true')
parser.add_argument('--prok', help='Flag to consider just prokaryotic proteins in Uniprot DB', action='store_true')
parser.add_argument('--euk', help='Flag to consider just eukaryotic proteins in Uniprot DB', action='store_true')
parser.add_argument('--sprot', help='Flag to exclude unreviewed Uniprot entries', action='store_true')
parser.add_argument('-nc', help='Number of cores to use for PSIBLAST and Clustal Omega calculations', type=num_cores, default=1, required=False)

args = vars(parser.parse_args())
nthreads = args['nc']


dbstr = '/home/resources/uniprot/databases/uniprot_all.fasta'
if args['euk']:
	dbstr = '/home/resources/uniprot/databases/uniprot_euk.fasta'
if args['sprot']:
	dbstr = '/home/resources/uniprot/databases/uniprot_sprot.fasta'

outstr = 'psiblast/rawoutb.txt'

def get_uniprot(line):
	bar1 = line.find('|')
	bar2 = line.find('|', bar1+1)
	return line[bar1+1:bar2]

def get_taxa(line):
	und = line.find('_')
	return line[und+1:]

eutaxa = set()
with open('/home/resources/uniprot/databases/eukaryota.dat', 'r') as efile:
	efile.readline()
	for line in efile:
		eutaxa.add(line.split('\t')[1])
efile.close()

protlist = list(set([l.split('\t')[0] for l in open(outstr,'r')]))
prokset = set()
eukset = set(protlist)

#
# STEP 3 - Parse PSIBLAST output
#

keys, unidb = parse_fasta(dbstr)

match_data = []
with open(outstr, 'r') as pfile:
	curr = None
	curr_list = []
	length = 0.
	curr_set = set()
	added_query = False

	cutoff = args['cutoff']
	# Use a match for Clustal Omega if all of the following apply:
	# 1) The e-value of the PSIBLAST match < cutoff
	# 2) The protein belongs to a eukaryote species, unless --prok is set
	# 3) No protein has already been added from this species, unless --mult is set
	# 4) The range of the match spans at least args['range'] of the query
	# 5) The protein is not a perfect match (unless it is the query)

	for line in pfile:
		info = line.split('\t')
		if info[0] == curr:
			tx = get_taxa(info[1])
			good1 = float(info[-2]) < cutoff
			good2 = (tx in eutaxa and info[0] in eukset) or (tx not in eutaxa and info[0] in prokset)
			good3 = tx not in curr_set or args['mult']
			good4 = float(info[3]) > args['range'] * length
			good5 = float(info[2]) < 100

			if good1 and good2 and good3 and good4 and good5:
				curr_list.append(info[1])
				curr_set.add(tx)
			elif get_uniprot(info[1]) == curr and not added_query:
				curr_list.insert(0, info[1])
				curr_set.add(tx)
				added_query = True
		else:
			if curr is not None:
				match_data.append((curr, curr_list))
			
			# This must be a perfect match
			length = float(info[3])
			curr_list = []
			curr_set = set()
			curr = info[0]
			added_query = False
			if get_uniprot(info[1]) == curr:
				curr_list.insert(0, info[1])
				curr_set.add(get_taxa(info[1]))
				added_query = True

	if curr is not None:
		match_data.append((curr, curr_list))

pfile.close()

for prot, prot_list in match_data:
	seqsfile = 'psiblast/clustal_input/%s.fasta' % prot
	with open(seqsfile, 'w') as inpfile:
		for match in prot_list:
			inpfile.write('> %s\n' % match)
			inpfile.write('%s\n' % unidb[match])
	inpfile.close()


#
# STEP 4 - Run Clustal Omega and Format Alignments
#

def do_command(protid):
	seqsfile = 'psiblast/clustal_input/%s.fasta' % protid
	clustalout = 'psiblast/%s_clustal.msa' % protid
	alignedout = 'psiblast/%s_aligned.msa' % protid

	if os.path.exists(alignedout):
		return
	
	if not os.path.exists(clustalout):
		numseqs = 0
		try:
			with open(seqsfile, 'r') as sf:
				for i, l in enumerate(sf):
					numseqs = i
			sf.close()
		except IOError:
			return
		numseqs /= 2

		if numseqs > 1:
			# Run Clustal Omega
			os.system('clustalo -i %s -o %s --wrap 100000 --force' % (seqsfile, clustalout))
			#subprocess.call(['clustalo', '-i', seqsfile, '-o', clustalout, '--wrap', '100000', '--force'])
		else:
			# Clustal Omega will fail, copy single sequence to output file
			os.system('cp %s %s' % (seqsfile, clustalout))
			#subprocess.call(['cp', seqsfile, clustalout])
	
	# Format Alignment

	try:
		with open(clustalout, 'r') as clout:
			outtxt = ''
			gaps = []
			for idx, line in enumerate(clout):
				line = line.strip()
				if idx % 2 == 0: # Header
					outtxt += line
					outtxt += '\n'
				elif idx == 1: # Query
					for i in range(len(line)):
						gaps.append(line[i] == '-')
				if idx % 2 == 1: # Query or Match
					newseq = ''
					for i in range(len(gaps)):
						if not gaps[i]:
							if i < len(line):
								newseq += line[i]
							else:
								newseq += '-'
					outtxt += newseq
					outtxt += '\n'
		clout.close()
		
		with open(alignedout, 'w') as alout:
			alout.write(outtxt)
		alout.close()
	except IOError:
		pass

my_pool = ThreadPool(nthreads)
command_outputs = my_pool.map(do_command, protlist)
