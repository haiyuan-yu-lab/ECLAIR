# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script generates multiple sequence Clustal Omega
# multiple sequence alignment files for all UniProt IDs
# in a given input file of interactions. The sequences
# to be included in each MSA are selected by filtering
# the results from PSIBLAST based on five conditions...
#
#    1) The e-value of the PSIBLAST match < cutoff
#    2) The protein belongs to a eukaryote species, unless --prok is set
#    3) No protein has already been added from this species, unless --mult is set
#    4) The range of the match spans at least args['range'] of the query
#    5) The protein is not a perfect match (unless it is the query)
#
# All Clustal Omega outputs are formatted to only
# include the positions that actually match the UniProt
# sequences of interest (i.e. remove all gaps with respect
# to the querry). There are currently no checks to
# avoid re-calculating MSA on previously handled UniProts
# and there are no checks for the consistency of the MSAs
# between runs.

# Expected Outcomes:
# All UniProt IDs in the input file should have their
# MSA generated and saved for later use. I do not fully
# understand the conditions under which a protein
# could fail to generate a MSA. My understanding is
# that even if there are no appropriate matches in
# the PSIBLAST output, every protein should at least
# generate a MSA that contains only itself.

# Known Bugs:
# - There are no known bugs, but there are some concerns
# - Unchecked edge cases based on our local uniprot
#   resources not matching up with the fetched uniprot
#   information?
# - Uncaught IOError in the final alignment parsing
# - Concerns about stocasticity between different
#   runs of PSIBLAST / how ties are handled?
# - Suspected potential to reduce PSI-BLAST runtime
#   by doing a PSI-BLAST on the whole of UniProt / saving
#   results?
# - Suspected potential to reduce need to regenerate
#   MSA by first checking if the sequences that would
#   be contained in the MSA have changed or not?

# Imports
import subprocess, argparse, sys, os, glob
from helper import *
from multiprocessing import Pool
from mjm_parsers import parse_fasta, parse_dictionary_list

# Convert x to float, assert x in range [0.0, 1.0]
# Only used for input parameter constraints
def restricted_float(x):
	x = float(x)
	if x < 0.0 or x > 1.0:
		raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % x)
	return x

# Process number of cores, assert x within multivac's capacity
# Only used for input parameter constraints
def num_cores(x):
	x = int(x)
	if x < 1 or x > 144:
		raise argparse.ArgumentTypeError('%d not in range [1, 144]' % x)
	return x

# Handle input parameters
parser = argparse.ArgumentParser(description='Find and align similar sequences to a protein using PSIBLAST and Clustal Omega.')
parser.add_argument('-f', metavar='/path/to/protfile', help='File containing list of input Uniprot ids')
parser.add_argument('-ints', metavar='/path/to/interfile', help='File containing list of interactions')
parser.add_argument('-p', metavar='UP', nargs='*', help='List of input Uniprot ids')
parser.add_argument('-c', '--cutoff', help='Cutoff e-value for all PSIBLAST matches used for alignment', type=float, default=0.05, required=False)
parser.add_argument('-r', '--range', help='Fraction of the query sequence required to consider PSIBLAST match for alignment', type=restricted_float, default=0.5)
parser.add_argument('--mult', help='Flag to include multiple proteins for single species', action='store_true')
parser.add_argument('--prok', help='Flag to consider just prokaryotic proteins in Uniprot DB', action='store_true')
parser.add_argument('--euk', help='Flag to consider just eukaryotic proteins in Uniprot DB', action='store_true')
parser.add_argument('--sprot', help='Flag to exclude unreviewed Uniprot entries', action='store_true')
parser.add_argument('-nc', help='Number of cores to use for PSIBLAST and Clustal Omega calculations', type=num_cores, default=1, required=False)

# Parse parameters
args = vars(parser.parse_args())
protstr = args['f']
protids = args['p']
interstr = args['ints']

# Ensure that either a input file or input list is provided
if protstr is None and protids is None:
	print 'Nothing inputted, nothing to do.\nUse python psiblast_update.py -h for help.'
	sys.exit()


# Select the UniProt database to use
dbstr = '/home/resources/uniprot/databases/uniprot_all.fasta'
if args['euk']:
	# Use only Eukaryotic proteins
	dbstr = '/home/resources/uniprot/databases/uniprot_euk.fasta'
if args['sprot']:
	# Use unreview UniProt entries
	dbstr = '/home/resources/uniprot/databases/uniprot_sprot.fasta'
# Why isn't there a check for the --prok flag?

# Presumably returns the UniProt identifier from an enty in the uniprot fasta?
def get_uniprot(line):
	# WHY?
	bar1 = line.find('|')
	bar2 = line.find('|', bar1+1)
	return line[bar1+1:bar2]

# Presumably return the taxon identifier from an entry in the uniprot fasta?
def get_taxa(line):
	# WHY?
	und = line.find('_')
	return line[und+1:]

#
# STEP 1 - Write all protein sequences to a file for PSIBLAST
#

# I HAVE NO CLUE WHAT THIS IS OR WHY IT IS COMMENTED OUT
# I WILL PRETEND IT DOES NOT EXIST
'''
protid_to_key = {}
for key in keys:
	bar1 = key.find('|')
	bar2 = key.find('|', bar1+1)
	protid = key[bar1+1:bar2]
	protid_to_key[protid] = key
'''

# The filename that we will use for the PSIBLAST input Fasta
seqsstr = 'seqfile.txt'

# Obtain a list of all Eukaryota taxa in our uniprot resources
# HEAD:
#Taxon   Mnemonic        Scientific name Common name     Synonym Other Names     Reviewed        Rank    Lineage ParentVirus hosts
#651506          Aa calceata                     Aa calceata (Rchb.f.) Schltr.   annotated       Species Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; Liliopsida; Asparagales; Orchidaceae; Orchidoideae; Cranichideae; Cranichidinae; Aa        152839
#415388          Aa colombiana                   Aa colombiana Schltr.   annotated       Species Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; Liliopsida; Asparagales; Orchidaceae; Orchidoideae; Cranichideae; Cranichidinae; Aa        152839
#415389          Aa hartwegii                    Aa hartwegii Garay      annotated       Species Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; Liliopsida; Asparagales; Orchidaceae; Orchidoideae; Cranichideae; Cranichidinae; Aa        152839
eutaxa = set()
with open('/home/resources/uniprot/databases/eukaryota.dat', 'r') as efile:
	efile.readline()
	for line in efile:
		# This isn't actually the Taxon ID?
		eutaxa.add(line.split('\t')[1])
efile.close()


# List of Proteins
protlist = []
# Set of Prokaryotic Proteins
prokset = set()
# Set of Eukaryotic Proteins
eukset = set()

# Parse List of Interactions from File
# Flatten the list to a set of unique UniProt IDs included
# in the user's interaction file
inters = [l.strip().split('\t')[:2] for l in open(interstr)]
inters = set([item for sublist in inters for item in sublist])

# This is always True (in the current use case at least)
if protstr is not None:
	
	# Parse UniProt info file as a list of dictionaries (each
	# row is one dictionary)
	uniprot_dict = parse_dictionary_list(protstr)
	
	# Write sequence of each relevant UniProt ID to PSIBLAST input
	# Fasta
	with open(seqsstr, 'w') as seqfile:
		for e in uniprot_dict:
			# Skip all UniProts not included in user's interactions
			if e['id'] not in inters:
				continue
			
			# Write UniProt to fasta
			seqfile.write('> %s\n' % e['id'])
			seqfile.write('%s\n' % e['sequence'])
			
			# Add UniProt to appropriate kingdom set (proks, vs. euks)
			# It seems like this is a haphazard way of defining prokaryotics
			# since it says anything that is not explicity eukaryotic is prokaryotic?
			# Are there other options? Suggest looking at UniProt documentation to
			# find more specific ID mapping here.
			if e['lineage-id(SUPERKINGDOM)'] == '2759':
				eukset.add(e['id'])
			else:
				prokset.add(e['id'])
			
			# Add UniProt to full protein list
			protlist.append(e['id'])
	
	# Close PSIBLAST input Fasta
	seqfile.close()

# I think this is just an alternate way of writting the
# PSIBLAST input Fasta if the Protein IDs are provided
# as a list instead of through a file
'''
if protids is not None:
	method = 'a' if len(protlist) > 0 else 'w'
	with open(seqsstr, method) as seqfile:
		for protid in protids:
			seqfile.write('> %s\n' % protid)
			seqfile.write('%s\n' % unidb[protid_to_key[protid]])
			protlist.append(protid)
	seqfile.close()
'''

#
# STEP 2 - Run PSIBLAST
#

# PSIBLAST Output File
outstr = 'psiblast/rawout.txt'

# Run PSIBLAST
os.system('psiblast -query %s -db %s -dbsize 20905025825 -outfmt 6 -out %s -num_threads %d' % (seqsstr, dbstr, outstr, args['nc']))
#subprocess.call(['psiblast', '-query', seqsstr, '-db', dbstr, '-outfmt', '6', '-out', outstr, '-num_threads', str(args['nc'])])

#
# STEP 3 - Parse PSIBLAST output
#


# Parse PSIBLAST database Fasta as a dictionary
keys, unidb = parse_fasta(dbstr)

# List of all UniProt IDs and their PSIBLAST matches
match_data = []

# Read PSIBLAST output
with open(outstr, 'r') as pfile:
	# Current UniProt being Processed
	curr = None
	# Current List of matches?
	curr_list = []
	# Length of Currnet UniProt
	length = 0.
	# Current Set of species that have a match for the Current UniProt
	curr_set = set()
	# Flag to indicate if the Current UniProt has already added a match to itself?
	added_query = False
	
	# User specified cutoff
	cutoff = args['cutoff']
	# Use a match for Clustal Omega if all of the following apply:
	# 1) The e-value of the PSIBLAST match < cutoff
	# 2) The protein belongs to a eukaryote species, unless --prok is set
	# 3) No protein has already been added from this species, unless --mult is set
	# 4) The range of the match spans at least args['range'] of the query
	# 5) The protein is not a perfect match (unless it is the query)
	
	# Iterate through lines in the PSIBLAST output (I think this is the same as regular BLAST format)
	# HEAD:
	#qseqid  sseqid                  pident  length  mis     gap     qstart  qend    sstart  send    evalue  bitscore
	#P31749  tr|G1S2V0|G1S2V0_NOMLE  100.00  480     0       0       1       480     1       480     0.0     1004
	#P31749  tr|G3RB32|G3RB32_GORGO  100.00  480     0       0       1       480     1       480     0.0     1004
	#P31749  tr|K7AGW5|K7AGW5_PANTR  100.00  480     0       0       1       480     1       480     0.0     1004
	#P31749  tr|B0LPE5|B0LPE5_HUMAN  100.00  480     0       0       1       480     1       480     0.0     1004
	#P31749  sp|P31749|AKT1_HUMAN    100.00  480     0       0       1       480     1       480     0.0     1004
	for line in pfile:
		# Split data out into liens
		info = line.split('\t')
		
		if("covid" in info[0].lower()):
			continue
		
		# Continue working on current UniProt ID if next line matches it
		if info[0] == curr:
			# Get the taxa of the reference match
			tx = get_taxa(info[1])
			
			# Check 5 condition above to see if this match can be retained
			# for Clustal Omega
			
			# 1) The e-value of the PSIBLAST match < cutoff
			good1 = float(info[-2]) < cutoff
			# 2) The protein belongs to a eukaryote species, unless --prok is set
			good2 = (tx in eutaxa and info[0] in eukset) or (tx not in eutaxa and info[0] in prokset)
			# 3) No protein has already been added from this species, unless --mult is set
			good3 = tx not in curr_set or args['mult']
			# 4) The range of the match spans at least args['range'] of the query
			good4 = float(info[3]) > args['range'] * length
			# 5) The protein is not a perfect match (unless it is the query)
			good5 = float(info[2]) < 100
			
			# Add the sseqid and taxa to the current list / set if all checks pass
			if good1 and good2 and good3 and good4 and good5:
				curr_list.append(info[1])
				curr_set.add(tx)
			# Otherwise, if the sseqid matches the current Protein, add this match?
			elif get_uniprot(info[1]) == curr and not added_query:
				curr_list.insert(0, info[1])
				curr_set.add(tx)
				added_query = True
		# Set up for a new UniProt ID if the next lines does not match the current
		else:
			# If there is a current UniProt, save its ID and list of matches
			if curr is not None:
				match_data.append((curr, curr_list))
			
			# I Don't know what this means? Maybe it is with regard to being
			# able to assume that the length of this match reflects the whole
			# length of the UniProt
			# # This must be a perfect match
			
			# Update all of the current UniProt variables
			length = float(info[3])
			curr_list = []
			curr_set = set()
			curr = info[0]
			added_query = False
			
			# If the sseqid matches the current Protein, add this match?
			if get_uniprot(info[1]) == curr:
				curr_list.insert(0, info[1])
				curr_set.add(get_taxa(info[1]))
				added_query = True
	
	# If there is a current UniProt, save its ID and list of matches
	if curr is not None:
		match_data.append((curr, curr_list))

# Close the PSIBLAST Output
pfile.close()

# Write individual Fasta file for each UniProt in the PSIBLAST
# matches list (should be a subset of the user's interactions?)
# This Fasta file will include all of the matches that will
# be used for clustal omega multiple sequence alignment.
# NOTE: May be able to add a check here to figure out if the
#       proteins involved in the multiple sequence alignment
#       have changed. If they have not, should not need to
#       regenerate the MSA for this protein?
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


# Possibly add a check that set(protlist) == inters
# Catches the case where a user UniProt ID somehow
# is not represented in our local resources, and should
# throw a warning. Might be better suited to include this
# earlier in the code.


# Iterate through every protein in the protlist (this should basically
# be all of the user's interaction UniProt IDs?)
for protid in protlist:
	# Identify the fasta file for this protein
	seqsfile = 'psiblast/clustal_input/%s.fasta' % protid
	
	# Identify the clustal output file for this protein
	clustalout = 'psiblast/%s_clustal.msa' % protid
	
	# Figure out how many matches were included for this protein
	# This could be done in other, better ways?
	numseqs = 0
	try:
		with open(seqsfile, 'r') as sf:
			for i, l in enumerate(sf):
				numseqs = i
		sf.close()
	except IOError:
		pass
	numseqs /= 2
	
	# General Clustal Omega MSA using the specified fastas using matches
	if numseqs > 1:
		# Run Clustal Omega
		os.system('clustalo -i %s -o %s --wrap 100000 --force' % (seqsfile, clustalout))
		#subprocess.call(['clustalo', '-i', seqsfile, '-o', clustalout, '--wrap', '100000', '--force'])
	else:
		# Clustal Omega will fail, copy single sequence to output file
		os.system('cp %s %s' % (seqsfile, clustalout))
		#subprocess.call(['cp', seqsfile, clustalout])
	
	# Format Alignment
	alignedout = 'psiblast/%s_aligned.msa' % protid
	
	try:
		# Read clustal MSA
		with open(clustalout, 'r') as clout:
			# Text for writting to formatted alignment
			outtxt = ''
			
			# Keeps track of gaps?
			gaps = []
			
			# Iterate over each line
			for idx, line in enumerate(clout):
				line = line.strip()
				
				# Add Header lines as they are
				if idx % 2 == 0:
					outtxt += line
					outtxt += '\n'
				# Special case for the first entry in the alignment (I think this should be the current UniProt?)
				# Find all of the gaps in the alignment since we only care
				# about using the MSA with regard to the current UniProt
				# querry, we don't care about any of the positions where
				# the querry has a gap
				elif idx == 1: # Query
					# Find all the Gaps
					for i in range(len(line)):
						gaps.append(line[i] == '-')
				
				# For all matches
				if idx % 2 == 1:
					# Update the sequence by removing all of the positions
					# that were a gap in the current UniProt alignment
					newseq = ''
					for i in range(len(gaps)):
						if not gaps[i]:
							if i < len(line):
								newseq += line[i]
							else:
								newseq += '-'
					
					# Write the formatted alignment sequence
					outtxt += newseq
					outtxt += '\n'
		# Close the clustal output
		clout.close()
		
		# Write all of the formatted alignment lines
		# to the final alignment output
		with open(alignedout, 'w') as alout:
			alout.write(outtxt)
		alout.close()
	# Should print some warning here?
	except IOError:
		pass
