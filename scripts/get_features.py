#! /home/jc2375/.conda/envs/adr66-copy/bin/python

import os, subprocess, glob, argparse, sys, time, logging

# Set up parser to accept input parameters
parser = argparse.ArgumentParser(description='Gather features for a list of Uniprot ids.')
parser.add_argument('-f', metavar='/path/to/protfile', help='Absolute path of file containing list of input interactions')
parser.add_argument('--prok', help='Run PSIBLAST on only prokaryotic organisms', action='store_true')
parser.add_argument('-nc', help='Number of cores to use for generating predictions', type=int, default=1, required=False)
parser.add_argument('-predpath', metavar='/path/to/prediction/folder', help='Absolute path to folder for prediction')
parser.add_argument('--step', help='Step to start from (default 1)', type=int, default=1, required=False)

# Parse input
args = vars(parser.parse_args())

# Helper functions for logging
def log_time(name, interval):
	with open('timelog.txt','a') as tf:
		tf.write('%s: %0.2f seconds\n' % (name, interval))

def run_with_logging(commands, logger, verbose=False):
	p = subprocess.Popen(commands,
			     stdout=subprocess.PIPE,
			     stderr=subprocess.PIPE)
	
	# Added option to print all outputs from all subprocesses rather
	# than only printing stderr
	if(verbose):
		stdout, stderr = p.communicate()
		logger.info(stdout)
		if stderr:
			logger.error(stderr)
	else:
		_, stderr = p.communicate()
		if stderr:
			logger.error(stderr)

# GLOBALS
inter_file = args['f']
uniprot_file = '/home/adr66/eclair/data/uniprot_info.txt'
LOGGER_NAME = 'get_features_logger'
LOG_DEST = '/home/jc2375/outputs/get_features.log'

#########################################
#     STEP 0 - Code Setup               #
#########################################

# DO NOT uncomment this unless you want to reset all files!
#run_with_logging(['python', 'reset_files.py'])

# Move to correct directory for execution
os.chdir('/home/adr66/eclair/')

# Initialize logging
logger = logging.getLogger(LOGGER_NAME)
handler = logging.FileHandler(LOG_DEST)
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

logger.info('------------- Running get_features.py ----------------')
logger.info('Step Selected: {0}'.format(args['step']))
logger.info('STEP 0 - Code Setup: Complete')

#########################################
#     STEP 1 - List new Uniprot Ids     #
#########################################

if args['step']<=1:
	os.chdir('data/uniprot/')
	logger.info('Running: update_uniprot_info.py')
	
	# Purpose:
	# Fetches UniProt information for all UniProt IDs included in the
	# inter_file that have not previously been fetched.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/uniprot_info.txt
	#
	# Expected Outcome:
	# The uniprot_file should contain entries for all UniProt IDs
	# included in the inter_file.
	#
	# Known Bugs:
	# - Can break with changes to the UniProt API.
	# - Concerned about caching / updates to existing UniProt entries.
	run_with_logging(['python', 'update_uniprot_info.py', inter_file, uniprot_file], logger)
	os.chdir('../..')
	
	logger.info('STEP 1 - List new Uniprot IDs: Complete')
#exit()
#########################################
#   STEP 2 - ExPasy and Pfam features   #
#########################################

if args['step']<=2:
	os.chdir('data/isoforms/')
	logger.info('Running: parse_isoforms.py')
	# Purpose:
	# Generates a semi-colon delimited isoform mask (0s indicate positions
	# with no isoform, 1s indicate positions with isforms) for all UniProt
	# IDs in the inter_file that have not previously been processed.
	#
	# Files Affected:
	# - /home/adr66/eclair/features/per_feature/isoforms.txt
	#
	# Expected Outcome:
	# The /home/adr66/eclair/features/per_feature/isoforms.txt should contain
	# entries for all of the UniProt IDs included in the inter_file.
	#
	# Known Bugs:
	# - Concerned about caching / updates to existing UniProt entries.
	# - Potential efficiency improvements have been identified in the code.
	run_with_logging(['python', 'parse_isoforms.py', uniprot_file], logger)
	
	os.chdir('../expasy/')
	logger.info('Running: populate_expasy_per_uniprot.py')
	# Purpose:
	# Calculate seven pre-selected expasy features for all UniProt
	# IDs in the inter_file that have not previously been processed.
	# Expasy features are simply a semi-colon delimited list of the
	# Amino acid scale values for each residue in the UniProt sequence.
	#
	# Files Affected:
	# - /home/adr66/eclair/features/per_feature/expasy_ACCE.txt
	# - /home/adr66/eclair/features/per_feature/expasy_AREA.txt
	# - /home/adr66/eclair/features/per_feature/expasy_BULK.txt
	# - /home/adr66/eclair/features/per_feature/expasy_COMP.txt
	# - /home/adr66/eclair/features/per_feature/expasy_HPHO.txt
	# - /home/adr66/eclair/features/per_feature/expasy_POLA.txt
	# - /home/adr66/eclair/features/per_feature/expasy_TRAN.txt
	#
	# Expected Outcome:
	# The corresponding expasy feature files should each contain entries for
	# all of the UniProt IDs included in the inter_file. The exception is for
	# any UniProt IDs shorter than 10 residues.
	#
	# Known Bugs:
	# - No bugs, confusing feature setup, potential for a larger
	#   lab resource to be developed based on this script.
	run_with_logging(['python', 'populate_expasy_per_uniprot.py', uniprot_file], logger)
	
	os.chdir('../domains/')
	logger.info('Running: parse_domains.py')
	# Purpose:
	# Generates a domain annotation mask for all UniProt
	# IDs in the inter_file that have not previously been processed.
	# The mask is a semi-colon delimited list of 0s, no domain, and
	# 1s, domain, for every position in the UniProt sequence.
	#
	# Files Affected:
	# - /home/adr66/eclair/features/per_feature/pfam_domains.txt
	#
	# Expected Outcome:
	# The output file should contain entries for all of the UniProt IDs
	# included in the inter_file. The exceptions are any UniProt ID that
	# either 1) is longer than 10000 residues, 2) has no domain annotations
	# or 3) has inconsistent information between UniProt and pfam.
	#
	# Known Bugs:
	# - No bugs, possible efficiency improvements have been identified
	#   in the code.
	run_with_logging(['python', 'parse_domains.py', uniprot_file], logger)
	
	os.chdir('../../')
	
	logger.info('STEP 2 - ExPasy and Pfam features: Complete')

#########################################
#    STEP 3 - Conservation Features     #
#########################################

if args['step']<=3:

	os.chdir('data/conservation/')
	psiblast_call = ['python', 'psiblast_update.py', '-f', uniprot_file, '-ints', inter_file, '-nc', str(args['nc'])]
	if args['prok']:
		psiblast_call.append('--prok')
	logger.info('Running: psiblast_update.py')
	# Purpose:
	# Generate Multiple Sequence Alignments for all UniProt
	# IDs in the inter_file. MSAs are formatted so that all
	# positions that are a gap with respect to the querry
	# UniProt are removed.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/conservation/seqfile.txt
	# - /home/adr66/eclair/data/conservation/psiblast/rawout.txt
	# - /home/adr66/eclair/data/conservation/psiblast/clustal_input/[UNIPROTID].fasta # For each UniProt in interactions
	# - /home/adr66/eclair/data/conservation/psiblast/[UNIPROTID]_clustal.msa         # For each UniProt in interactions
	# - /home/adr66/eclair/data/conservation/psiblast/[UNIPROTID]_aligned.msa         # For each UniProt in interactions
	#
	# Expected Outcome:
	# A corresponding [UNIPROTID]_aligned.msa file shoulud be generated
	# for every UniProt ID included in the inter_file. I do not understand
	# the cases where this could fail to happen. Additionally, every
	# UniProt ID should generate a [UNIPROTID]_clustal.msa file and a
	# [UNIPROTID].fasta containing the raw MSA and clustal input Fasta
	# respectively. The PSIBLAST input file, seqfile.txt, should be
	# updated to include only the UniProt IDs in the inter_file, and
	# a new output file, rawout.txt, should have been generated with
	# this input.
	#
	# Known Bugs:
	# - What if any of the user's interaction UniProt IDs for some
	#   reason don't show up in our local resources folder? Where is
	#   not currently way to catch this and a warning should be added?
	# - Can probably add some checks to avoid re-calculating previous
	#   features. E.g. only calculate new MSA if the list of sequences
	#   included has changed? What about stocasticity between different
	#   runs of the PSIBLAST / tied alignments?
	# - Should print warning on IOError when generating final aligned.msa files?
	run_with_logging(psiblast_call, logger)
	#exit()
	logger.info('Running: batch_js_conservation.py')
	
	# Purpose:
	# Calculates the JS Conservation for all UniProt IDs in the
	# inter_file. The JS Conservation file is formatted as a
	# semi-colon delimited list of the JS Scores for each position
	# in the UniProt.
	#
	# Files Affected:
	# - /home/adr66/eclair/features/per_feature/JSconserve.txt
	# - /home/adr66/eclair/data/conservation/slimmed_msas/[UNIPROTID].msa         # For each UniProt in interactions
	#
	# Expected Outcome:
	# All UniProt IDs in the inter_file should have their JS
	# Conservation scores added to the output file. Slimmed MSA
	# files should be generated for every MSA file. The exception
	# is any UniProt ID that has fewer than 50 alignments remaining
	# in the slimmed MSA.
	#
	# Known Bugs:
	# - Concern about how it skips previously calculated features?
	#   What if the MSA has changed?
	# - Suggest throwing warnings for some edge cases rather than just
	#   print statements (depends on how common they are?)
	# - Suggest a spot check that JS Score positions are being mapped
	#   correctly
	run_with_logging(['python', 'batch_js_conservation.py', uniprot_file], logger)
	
	logger.info('Running: batch_js_conservation_rank.py')
	# Purpose:
	# Ranks the JS Scores calculated previously such that
	# the most conserved position is ranked first and the 
	# least conserved position is ranked last. In the case
	# of ties, each tied position receives the same rank
	# and the ranking skips ahead to account for the number
	# of ties. Outputs are written in semi-colon delimited
	# format, ranking every position in the protein.
	#
	# Files Affected:
	# - /home/adr66/eclair/features/per_feature/JSconserve_rank.txt
	#
	# Expected Outcome:
	# - All UniProt IDs that were able to calculate the regular
	#   JS Conservation feature should have their ranked
	#   conservation scores for each position added to the
	#   output file.
	#
	# Known Bugs:
	# - None
	run_with_logging(['python', 'batch_js_conservation_rank.py'], logger)
	#exit()
	
	os.chdir('../sca/')
	
	logger.info('Running: batch_psiblast_sca.py')
	# Purpose:
	# Uses Statistical Coupling Analysis (SCA) to calculate the
	# positional correlation matrix for all itneraction pairs
	# included in the inter_file.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/sca/joined_msas/[UNIPROTID1]_[UNIPROTID2].msa         # For each UniProt in interactions
	# - /home/adr66/eclair/data/sca/sca_corr_psiblast/[UNIPROTID1]_[UNIPROTID2].sca   # For each UniProt in interactions
	# - /home/adr66/eclair/data/sca/sca_timelog.txt
	#
	# Expected Outcome:
	# Joined MSA and SCA output files should be generated
	# for each interaction included in the inter_file. The
	# exception to this is interactions with fewer than 50
	# sequences in their final joined MSA or with joined
	# MSA alignment length exceeding 10000.
	#
	# Known Bugs:
	# - Efficiency: We regenerate the Slimmed MSAs
	#   instead of just reusing them
	# - Why is the coverage threshold set to 0?
	run_with_logging(['python', 'batch_psiblast_sca.py', inter_file], logger)
	
	logger.info('Running: batch_dca.py')
	# Purpose:
	# Uses Direct Couplin Analysis (DCA) to calculate for all
	# the Mutual and Direct information between any pair
	# of position for all of the interaction pairs included
	# in the inter_file.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/sca/dca/[UNIPROTID1]_[UNIPROTID2].dca   # For each UniProt in interactions
	# - /home/adr66/eclair/data/sca/dca/dca_timelog.txt
	#
	# Expected Outcome:
	# DCA files should be generated for each interaction
	# included in the inter_file. I'm not sure, but I think
	# that the preceding exceptions for betch_spiblast_sca.py
	# still apply (not hardcoded, but I think the proper input
	# will be missing if the interaction did not pass the
	# previous filter)?
	#
	# Known Bugs:
	# - Maximum timeout on DCA processes sometimes
	#   fails and script continues indefinitely / must
	#   be manually killed.
	# - Poor inconsistency between using a hard coded
	#   num_threads = 10 and and actual hard coded 10.
	#   This could cause some problems if these values
	#   were ever changed / is generally not good practice
	# - Why is there no equivalent check that the MSA have
	#   a minimum number of alignments or that the length
	#   of the alignments be below a certain limit? This
	#   might have something to do with the joined MSAs
	#   that were filtered based on these criteria never
	#   having been written in the previous batch_psiblast_sca.py
	#   step, so the filter carries over indirectly?
	run_with_logging(['python', 'batch_dca.py', inter_file], logger)
	
	logger.info('Running: parse_sca_results.py')
	# Purpose:
	# Calculates the final six features to be used in
	# Eclair based on different aggregates of the SCA
	# results. All features are calculated per residue
	# and saved in semi-colon delimited lists. Since
	# these features are specific to the interaction,
	# the the files are formatted as...
	#
	# P1	P2	Features (P1)	Features (P2)
	#
	# Files Affected:
	# - /home/adr66/eclair/features/per_feature/SCA_PB_max.txt
	# - /home/adr66/eclair/features/per_feature/SCA_PB_max_norm.txt
	# - /home/adr66/eclair/features/per_feature/SCA_PB_max_rank.txt
	# - /home/adr66/eclair/features/per_feature/SCA_PB_max_rank_percentile.txt
	# - /home/adr66/eclair/features/per_feature/SCA_PB_top10.txt
	# - /home/adr66/eclair/features/per_feature/SCA_PB_all.txt
	#
	# Expected Outcome:
	# The corresponding output files for each of the six
	# SCA features should be updated to include entries
	# for each interaction that was able to generate
	# SCA results.
	#
	# Known Bugs:
	# - Some warnings should be thrown that are not
	# - There is a bit in here where features are written
	#   based on explicit sorting of the UniProt IDs. I have
	#   not seen this anywhere else in the pipeline. I can't
	#   believe this has not been caught if it introduces an
	#   error so it must be ok?
	run_with_logging(['python', 'parse_sca_results.py', inter_file, uniprot_file], logger)
	
	logger.info('Running: split_dca_mi_di.py')
	# Purpose:
	# Separates the DCA output into two separate files
	# one for the Direct Information and one for the
	# Mutual Information. 
	#
	# Files Affected:
	# - /home/adr66/eclair/data/sca/mi/[UNIPROTID1]_[UNIPROTID2].mi   # For each UniProt in interactions
	# - /home/adr66/eclair/data/sca/di/[UNIPROTID1]_[UNIPROTID2].di   # For each UniProt in interactions
	#
	# Expected Outcome:
	# All interactions in the inter_file that had
	# successfully generated DCA outputs previously
	# should have these outputs split to generate
	# the corresponding .mi and .di files.
	#
	# Known Bugs:
	# - No bugs, some pulling of old results from Juan's
	#   direcoties
	run_with_logging(['python', 'split_dca_mi_di.py', inter_file], logger)
	
	logger.info('Running: parse_sca_results2.py')
	# Purpose:
	# Calculates the final ten features to be used in
	# Eclair based on different aggregates of the DCA
	# results. All features are calculated per residue
	# and saved in semi-colon delimited lists. Since
	# these features are specific to the interaction,
	# the the files are formatted as...
	#
	# P1	P2	Features (P1)	Features (P2)
	#
	# Files Affected:
	# - /home/adr66/eclair/features/per_feature/DDI_max.txt
	# - /home/adr66/eclair/features/per_feature/DDI_mean.txt
	# - /home/adr66/eclair/features/per_feature/DDI_top10.txt
	# - /home/adr66/eclair/features/per_feature/DDI_max_localZ.txt
	# - /home/adr66/eclair/features/per_feature/DDI_max_globalZ.txt
	# - /home/adr66/eclair/features/per_feature/DMI_max.txt
	# - /home/adr66/eclair/features/per_feature/DMI_mean.txt
	# - /home/adr66/eclair/features/per_feature/DMI_top10.txt
	# - /home/adr66/eclair/features/per_feature/DMI_max_localZ.txt
	# - /home/adr66/eclair/features/per_feature/DMI_max_globalZ.txt
	#
	# Expected Outcome:
	# The corresponding output files for each of the ten
	# DCA features should be updated to include entries
	# for each interaction that was able to generate
	# DCA results.
	#
	# Known Bugs:
	run_with_logging(['python', 'parse_sca_results2.py', uniprot_file, inter_file], logger)
	os.chdir('../../')
	
	logger.info('STEP 3 - Conservation Features: Complete')


#exit()

#########################################
#        STEP 4 - MODBASE Models        #
#########################################

if args['step']<=4:
	os.chdir('data/modbase/')
	logger.info('Running: modbase_DL.py')
	# Purpose:
	# Fetches ModBase information for each Uniprot ID
	# included in the inter_file. Note, the information
	# files downloaded here are still in raw XML format
	# with many tags for unique <pdbfiles> that serve
	# as ModBase models. Considerable parsing is sill
	# required before this information may be used.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/modbase/models/uniprot/[UNIPROTID].pdb   # For each UniProt in interactions
	#
	# Expected Outcome:
	# Each UniProt ID that has information available
	# in ModBase should have its data downloaded and
	# saved to a corresponding output file.
	#
	# Known Bugs:
	# - None
	run_with_logging(['python', 'modbase_DL.py', uniprot_file, 'models/uniprot'], logger)
	
	logger.info('Running: modbase_02.py')
	# Purpose:
	# Begins to process previously fetched raw ModBase
	# XML downloads. Performs a first pass through the
	# data to identify all of the valid models for each
	# UniProt ID, and saves their PDB content and header
	# information.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/modbase/models/header/[MODELID].txt   # For each valid model of each UniProt in interactions
	# - /home/adr66/eclair/data/modbase/models/hash/[MODELID].pdb     # For each valid model of each UniProt in interactions
	#
	# Expected Outcome:
	# Header and PDB files should be created for all valid
	# models in the raw ModBase download. Models may
	# be deemed invalid for 1) having UniProt information
	# (i.e. length) that is inconsistent with the local
	# UniProt info file, 2) having multiple remarks for
	# a certain header value, or 3) having improperly
	# formatted header values.
	#
	# Known Bugs:
	# - Should implement warnings
	# - Why is the code to remove previously
	#   downloaded models that are no longer
	#   valid commented out? Probably related
	#   to wanting to keep the cache / not lose
	#   it just because a certain UniProt ID
	#   was not included in a new run?
	run_with_logging(['python', 'modbase_02.py', 'models/uniprot', 'models/hash', 'models/headers', uniprot_file], logger)
	
	logger.info('Running: modbase_03.py')
	# Purpose:
	# Continues parsing the ModBase information for each
	# UniProt. In this step all of the valid models
	# remaining from modbase_02.py are processed to make
	# sure that their indices match up directly with the
	# position indices for the position in the UniProt ID
	# that is being modelled.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/modbase/models/hash/[MODELID].pdb     # For each valid model of each UniProt in interactions
	#
	# Expected Outcome:
	# All of the previously generated PDB files for
	# valid models left after modbase_02.py should
	# be updated to have correct indexing wherever
	# the original downloads did not match up 1:1 with
	# the UniProt indices.
	#
	# Known Bugs:
	# - Should implement warnings
	# - How do we catch bad_indices templates?
	# - How do we know that just because the first
	#   index is correct ALL of the indices will be
	#   correct?
	# - It looks like the way that the corrected index
	#   is calculated may be wrong? Shouldn't it add
	#   (current_index - original_first_index) rather
	#   than just current_index?
	run_with_logging(['python', 'modbase_03.py'], logger)
	
	logger.info('Running: modbase_04.py')
	# Purpose:
	# Continues parsing the ModBase information for each
	# Uniprot. In this step all of the headers for the valid
	# models remaining from modbase_02.py are pooled together
	# into two summary files; a complete summary, and a
	# non-redundant, high quality summary file.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/modbase/models/parsed/all_modbase_models.txt
	# - /home/adr66/eclair/data/modbase/models/parsed/select_modbase_models.txt
	#
	# Expected Outcome:
	# The all models and select models summary files should
	# be updated to include the appropriate header information
	# for all of the valid models corresponding to the UniProt
	# IDs included in the inter_file.
	#
	# Known Bugs:
	# - No bugs, but questions...
	# - For filtering of redundant models in the
	#   select_modbase_models file are there ever
	#   cases where there are two overlapping, but
	#   not completely overlapping models?
	run_with_logging(['python', 'modbase_04.py'], logger)
	
	logger.info('Running: sres_05.py')
	# Purpose:
	# Calculates the SASA for each model in the select
	# model set generated previously in modbas_04. It is
	# currently unclear to me if this actually does anything
	# to look at surface residues (e.g. only returns SASA
	# values for the surface residues) or is just a SASA
	# calculator being called a surface residue calculator.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/modbase/SASA_modbase.txt
	#
	# Expected Outcome:
	# All new models included in the select models summary
	# file generated previously should have their SASA
	# values for all positions calculated / written to the
	# output file. Output is a semi-colon delimited list
	# of SASA values, all positions in the UniProt not covered
	# by a given model are marked NaN.
	#
	# Known Bugs:
	# - Should have warnings
	run_with_logging(['python', 'sres_05.py', uniprot_file, inter_file], logger)
	
	logger.info('Running: batch_modbase_sasa.py')
	# Purpose:
	# Compiles the individual SASA_modbase information
	# generated for each ModBase model previously
	# in sres_05.py to provide an aggregate SASA value
	# per residue per UniProt. The average of the SASA
	# values for each model that covers a particular
	# position is reported at each position in the
	# UniProt. The output is reported in semi-colon
	# delimited lists. For some reason this output
	# is saved per interaction rather than per
	# protein.
	#
	# Files Affected:
	# - /home/adr66/eclair/features/per_feature/SASA_MB.txt
	#
	# Expected Outcome:
	# The output file should be updated to include SASA
	# estimates for all interactions listed in the
	# inter_file that had valid ModBase models.
	#
	# Known Bugs:
	# - Interaction order sorting again (code is under the
	#   assumption that interactions should already be
	#   sorted, but does it again anyway. Nothing I have
	#   seen indicates that there is any restriction on
	#   the interaction pairs being sorted unless the user
	#   imposed that restriction on the input manually)
	# - Throwing warnings
	# - Confusion about why the output here is per interaction
	#   rather than per protein
	run_with_logging(['python', 'batch_modbase_sasa.py', uniprot_file, inter_file], logger)
	os.chdir('../..')
	
	logger.info('STEP 4 - MODBASE Models: Complete')

#########################################
#        STEP 5 - PDB Structures        #
#########################################

if args['step']<=5:
	os.chdir('data/pdb/')
	logger.info('Running: batch_pdb_sasa.py')
	# Purpose:
	# Calculates the average SASA per residue for each
	# UniProt ID in each interaction in the inter_file
	# as is done in batch_modbase_sasa.py, except that
	# this time it uses precomputed SASA information
	# for the whole PDB. The output is reported in
	# semi-colon delimited lists. For some reason this
	# output is saved per interaction rather than per
	# protein.
	#
	# Files Affected:
	# - /home/adr66/eclair/features/per_feature/SASA_PDB.txt
	#
	# Expected Outcome:
	# The output file should be updated to include SASA
	# estimates for all interactions listed in the
	# inter_file that mapped to PDB structures.
	#
	# Known Bugs:
	# - Depends on several static copies of files from out
	#   resources folder. I am now also concerned that I may
	#   have missed some dependence on static resources copies
	#   in my first pass through the earlier files
	# - Sorts interaction pair order
	# - Should add warnings
	run_with_logging(['python', 'batch_pdb_sasa.py', uniprot_file, inter_file], logger)
	
	logger.info('Running: batch_combined_sasa.py')
	# Purpose:
	# Calculates the combined SASA per residue for
	# each UniProt ID in each interaction in the
	# inter_file. This combined SASA is reported
	# as both the maximum and the average SASA
	# per residue as reported in either ModBase
	# or PDB. The output is reported in semi-colon
	# delimited lists. For some reason this output
	# is saved per interaction rather than per
	# protein.
	#
	# Files Affected:
	# - /home/adr66/eclair/features/per_feature/SASA_combined_max.txt
	# - /home/adr66/eclair/features/per_feature/SASA_combined_avg.txt
	#
	# Expected Outcome:
	# The output files should be updated to include SASA
	# estimates for all interactions listed in the
	# inter_file that had PDB or ModBase structures.
	#
	# Known Bugs:
	# - Depends on static copies of files from the
	#   resources folder.
	# - Sorts interaction pair order (I am now feeling
	#   somewhat more confidence that this is not
	#   a problem, because if looks like I may have
	#   been missing that the interactions are
	#   sorted each time the inter_file is read?)
	# - Add warnings
	run_with_logging(['python', 'batch_combined_sasa.py', uniprot_file, inter_file], logger)
	os.chdir('../../')
	
	logger.info('STEP 5 - PDB Structures: Complete')
#exit()
#########################################
#           STEP 6 - DOCKING            #
#########################################

if args['step']<=6:
	
	os.chdir('data/zdock')
	
	logger.info('Running: compile_modbase_docking_set.py')
	# Purpose:
	# Generates the ModBase docking sets to be used for
	# performing docking simulations for each of the
	# interactions included in the inter_file. The
	# docking sets contain a single pair of the highest
	# coverage ModBase models for each UniProt protein
	# in the interaction.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/zdock/modbase_docking_set.txt
	#
	# Expected Outcome:
	# The output file should be updated to include entries
	# for every interaction in the inter_file that has
	# valid ModBase Models.
	#
	# Known Bugs:
	# - Dependence on static resources copies
	# - Add warning
	# - Question whether using a single docking simulation
	#   per interaction is really the best option. Why
	#   not do multiple docks to achieve maximum total
	#   coverage and then combine?
	#run_with_logging(['python', 'compile_modbase_docking_set.py', uniprot_file, inter_file], logger)
	
	logger.info('Running: compile_pdb_docking_set.py')
	# Purpose:
	# Generates the PDB docking sets to be used for
	# performing docking simulations for each of the
	# interactions included in the inter_file. The
	# docking sets contain a single pair of the highest
	# coverage PDB structures for each UniProt protein
	# in the interaction.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/zdock/pdb_docking_set.txt
	#
	# Expected Outcome:
	# The output file should be updated to include entries
	# for every interaction in the inter_file that has
	# valid PDB Structures.
	#
	# Known Bugs:
	# - Dependence on static resource copies
	# - Questions about skipping large complexes as
	#   docking options
	#run_with_logging(['python', 'compile_pdb_docking_set.py', uniprot_file, inter_file], logger)
	
	logger.info('Running: compile_mixed_docking_set.py')
	# Purpose:
	# Generates a mixed ModBase / PDB docking set to be
	# used for performing docking simulations for each of
	# the interactions included in the inter_file. The
	# docking sets contain a single pair of the highest
	# coverage combination of one ModBase and one PDB
	# structure to represent the interaction.
	#
	# Files Affected:
	# /home/adr66/eclair/data/zdock/mixed_docking_set.txt
	#
	# Expected Outcome:
	# The output file should be updated to inlclude entries
	# for every interaction in the inter_file that has
	# valid PDB and ModBase Structures.
	#
	# Known Bugs:
	# - It looks like the select_modbase_models
	#   file that it uses is not actually the file
	#   that is actively being updated by Eclair
	#   but rather an old one in mjm659
	# - Dependence on static resource copies
	# - Questions about skipping large compelxes
	#run_with_logging(['python', 'compile_mixed_docking_set.py', uniprot_file, inter_file], logger)
	
	logger.info('Running: batch_modbase_docking.py')
	# Purpose:
	# Performs docking simulations on all of the docking
	# pairs included in the modbase_docking_set.txt
	# list generated in compile_modbase_docking_set.py.
	# Docking is carried out using a wrapper for Zdock.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/zdock/modbase_docked_models/[MODELID1]--[MODELID2]--ZDOCK.out    # For each pair of models included in the docking set
	# - /home/adr66/eclair/data/zdock/modbase_docked_models/[MODELID1]--[MODELID2]--ZDOCK-XX.pdb # For each pair of models included in the docking set. Where XX ranges from 01-10
	# - /home/adr66/eclair/data/zdock/modbase_timelog.txt
	#
	# Expected Outcome:
	# Docking results should be generated for each of
	# the docking pairs included in the docking set.
	#
	# Known Bugs:
	# - Zdock Errors? Need to look into. I suspect if this is a true error it affects all
	#   of the Zdock steps.
	#		head: cannot open '/tmp/tmpLeVtRq/zdock.out' for reading: No such file or directory
	#		Traceback (most recent call last):
	#		File "/home/resources/mjm_tools/wrapped_zdock.py", line 98, in <module>
	#		copyfile(zdock_output_file, new_zdock_outfile)
	#		File "/home/adr66/.conda/envs/jp/lib/python2.7/shutil.py", line 82, in copyfile
	#		with open(src, 'rb') as fsrc:
	#		IOError: [Errno 2] No such file or directory: '/tmp/tmpLeVtRq/zdock.out'
	#run_with_logging(['python', 'batch_modbase_docking.py', inter_file, str(args['nc'])], logger)
	
	logger.info('Running: batch_pdb_docking.py')
	# Performs docking simulations on all of the docking
	# pairs included in the pdb_docking_set.txt
	# list generated in compile_pdb_docking_set.py.
	# Docking is carried out using a wrapper for Zdock.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/zdock/pdb_docked_models/[MODELID1]--[MODELID2]--ZDOCK.out    # For each pair of models included in the docking set
	# - /home/adr66/eclair/data/zdock/pdb_docked_models/[MODELID1]--[MODELID2]--ZDOCK-XX.pdb # For each pair of models included in the docking set. Where XX ranges from 01-10
	# - /home/adr66/eclair/data/zdock/pdb_timelog.txt
	#
	# Expected Outcome:
	# Docking results should be generated for each of
	# the docking pairs included in the docking set.
	#
	# Known Bugs:
	# - None
	#run_with_logging(['python', 'batch_pdb_docking.py', inter_file, str(args['nc'])], logger)
	
	logger.info('Running: batch_mixing_docking.py')
	# Performs docking simulations on all of the docking
	# pairs included in the mixed_docking_set.txt
	# list generated in compile_mixed_docking_set.py.
	# Docking is carried out using a wrapper for Zdock.
	#
	# Files Affected:
	# - /home/adr66/eclair/data/zdock/mixed_docked_models/[MODELID1]--[MODELID2]--ZDOCK.out    # For each pair of models included in the docking set
	# - /home/adr66/eclair/data/zdock/mixed_docked_models/[MODELID1]--[MODELID2]--ZDOCK-XX.pdb # For each pair of models included in the docking set. Where XX ranges from 01-10
	# - /home/adr66/eclair/data/zdock/mixed_timelog.txt
	#
	# Expected Outcome:
	# Docking results should be generated for each of
	# the docking pairs included in the docking set.
	#
	# Known Bugs:
	# - None
	#run_with_logging(['python', 'batch_mixed_docking.py', inter_file, str(args['nc'])], logger)
	
	docks = ['mb','pdb','mixed']
	types = ['ires','dSASA','dist3d']
	
	for dock in docks:
		for type in types:
			call_array = ['python', 'calc_%s_%s.py' % (dock, type), inter_file]
			if dock != 'mb':
				call_array.append(uniprot_file)
			logger.info('Running: calc_%s_%s.py' % (dock, type))
			# Purpose:
			# Calculates the interface residues, delta
			# SASA, and 3D atomic distances features, on
			# either the ModBase, PDB, or mixed docking
			# set results. 3D atomic distances are reported
			# as the average distance between the CA of the
			# residue in question to all of the CA atoms
			# on the other protein.
			#
			# Files Affected:
			# - /home/adr66/eclair/data/zdock/ires/IRES_mb_docking.txt
			# - /home/adr66/eclair/data/zdock/ires/IRES_pdb_docking.txt
			# - /home/adr66/eclair/data/zdock/ires/IRES_mixed_docking.txt
			# - /home/adr66/eclair/data/zdock/ires/dSASA_mb_docking.txt
			# - /home/adr66/eclair/data/zdock/ires/dSASA_pdb_docking.txt
			# - /home/adr66/eclair/data/zdock/ires/dSASA_mixed_docking.txt
			# - /home/adr66/eclair/data/zdock/ires/dist3d_mb_docking.txt
			# - /home/adr66/eclair/data/zdock/ires/dist3d_pdb_docking.txt
			# - /home/adr66/eclair/data/zdock/ires/dist3d_mixed_docking.txt
			#
			# Expected Outcome:
			# All of the relevant features for all of the
			# interactions included in the inter_file should
			# be written to the relevant output files.
			#
			# Known Bugs:
			# - irescalc.py is known to fail to calculated
			#   interface residues on particularly large
			#   protein structures. I suspect that the fix
			#   I added to report this error may actually cause
			#   the calc_X_ires.py scripts to crash.
			# - Mapping using the SIFTs file in calc_pdb_ires.py
			#   only maps any given PDB chain to the highest
			#   coverage UniProt ID. I have a feeling this could
			#   drop some information on UniProt IDs whose best
			#   structural match is a suboptimal match?
			# - There is a difference between the dSASA value
			#   calculation and the ires calculation in which
			#   we take a step to explicitly map the ModBase
			#   residues back to UniProt residues (make sure
			#   we only include residues that are included in
			#   UniProt / that we set NaN for any values
			#   not covered by the template). I am not sure
			#   why we do not need to do this same mapping
			#   on the interface residues. We can definitely
			#   have UniProt residues that do not map to
			#   the ModBase, but can we have ModBase residues
			#   that do not map to the UniProt? Probably
			#   not? Need to look into what the ModBase
			#   models actually are.
			# - Due to the repetitive nature of this code
			#   some of the varriable names are not reflective
			#   of the actual values being calculated (e.g.
			#   idist3D is still called dSASA in the idist3D
			#   scripts). The same hold for some of my comments
			#   where I did not catch places to change as I
			#   copied them over between scripts.
			# - Looks like the ires features calculated
			#   here are never actually used downstream?
			# - Looks like the Mixed docking set is
			#   never actually used downstream?
			#run_with_logging(call_array, logger)
	
	logger.info('Running: parse_combined_dSASA.py')
	# Purpose:
	# Calculates the final aggregated zDOCK dSASA
	# based Eclair features for all of the interactions
	# included in the inter_file. This aggregate
	# feature set includes three values; 1) the
	# average dSASA value per residue reported
	# over all of the docking results, 2) the maximum
	# dSASA value per residue reported over all of
	# the docking results, 2) the dSASA value per
	# residue reported by the highest ranking
	# docking result.
	#
	# Files Affected:
	# - /home/adr66/eclair/features/per_feature/zdock_combined_0cut_avg.txt
	# - /home/adr66/eclair/features/per_feature/zdock_combined_0cut_max.txt
	# - /home/adr66/eclair/features/per_feature/zdock_combined_0cut_top1.txt
	#
	# Expected Outcome:
	# The output files should be updated to
	# include feature calculations for all of the
	# interactions included in the inter_file
	# that had docking simulations performed on
	# them.
	#
	# Known Bugs:
	# - Confused why these features are not split to
	#   have separate weights for ModBase vs. PDB based
	#   docking results rather than joining all of them
	# - It looks like the Mixed docking set is never
	#   used to compile any of the final features?
	# - It also looks like the ires values are never
	#   going to be used to compile any of the final
	#   features?
	# - I am now confused as to whether or not this
	#   calculation is based on the 10 docking results
	#   from a single docking simulation or if it can
	#   encorporate multiple docking simulations. If
	#   it can encorporate multiple, I need to check to
	#   make sure that the true highest scored docking
	#   result accross all simulations gets used for the
	#   top1 feature.
	run_with_logging(['python', 'parse_combined_dSASA.py'], logger)
	
	logger.info('Running: parse_priority_dist3d.py')
	# Purpose:
	# Calculates the final aggregated zDOCK dist3d
	# based Eclair features for all of the interactions
	# included in the inter_file. This aggregate
	# feature set included four values; 1) average,
	# 2) max, 3) min, 4) top docking result.
	#
	# Files Affected:
	# - /home/adr66/eclair/features/per_feature/zdock_dist3d_PRIORITY_0cut_avg.txt
	# - /home/adr66/eclair/features/per_feature/zdock_dist3d_PRIORITY_0cut_max.txt
	# - /home/adr66/eclair/features/per_feature/zdock_dist3d_PRIORITY_0cut_min.txt
	# - /home/adr66/eclair/features/per_feature/zdock_dist3d_PRIORITY_0cut_top1.txt
	# - /home/adr66/eclair/data/zdock/ires/priority_dist3d_zdock_0cut_0res.txt
	#
	# Expected Outcome:
	# The output files should be updated to
	# include feature calculations for all of the
	# interactions included in the inter_file
	# that had docking simulations performed on
	# them.
	#
	# Known Bugs:
	# - Includes static copy of ires_perppi and
	#   interactomes from our resources folder
	# - It looks like this script explicitly refuses
	#   to calculate these features on any interaction
	#   pairs that have evidence of being a true
	#   interaction based on co-crystal or HINT (?)
	#   evidence.
	# - In the part that randomly only looks at co-crystal
	#   structure information, it aggregates the dist3D
	#   values for homodimers by taking the MAXIUMUM
	#   rather than the MIMUMUM distance. This differs from
	#   the way that it was first done on the PREDICTION
	#   set. I don't know if these values actually ever
	#   get used or not.
	# - Conceptually I have no clue what the purpose of
	#   separating out the CC vs PREDICTED interactions
	#   in the beginning 66% of this script it. It all
	#   seems like it can be cut out to me, but I still
	#   need to confirm that there wasn't an informed
	#   purpose for doing all of it.
	# - Am concerned that all of the outputs are opened
	#   for writting, but it looks like the "already
	#   computed" results are never re-written. I expect
	#   that all old results for this set may be lost?
	#   Based on looking at the length of the ouput
	#   files, I don't see how that could be happening.
	run_with_logging(['python', 'parse_priority_dist3d.py', inter_file], logger)
	os.chdir('../..')
	
	logger.info('STEP 6 - DOCKING: Complete')

#########################################
#    STEP 7 - Make Interaction Files    #
#########################################

if args['step']<=7:
	os.chdir('/home/adr66/eclair/features')
	logger.info('Running: compile_features.py')
	# Purpose:
	# Compiles all of the features generated in the
	# rest of the pipeline into the final input
	# features to be used in classification.
	#
	# Files Affected:
	# - /home/adr66/eclair/features/per_interaction/[UNIPROTID1]_[UNIPROTID2]_0.pkl   # For each UniProt in interactions
	# - /home/adr66/eclair/features/per_interaction/[UNIPROTID1]_[UNIPROTID2]_1.pkl   # For each UniProt in interactions
	#
	# Expected Outcome:
	# All interactions included in the inter_file
	# should have their final compiled feature
	# files generated.
	#
	# Known Bugs:
	# - Confused why this script does not use all
	#   of the features that have been calculated
	# - Concerned about the way the script refuses
	#   to every update the compiled features once
	#   an interaction has been calculated once
	# - The save_ires_column flag seems to be
	#   redundant since the script seems to loop
	#   over all of the features (including IRES)
	#   afterwards anyway
	run_with_logging(['python', 'compile_features.py', inter_file, uniprot_file], logger)
	os.chdir('..')
	
	logger.info('STEP 7 - Make Interaction Files: Complete')

#########################################
#           STEP 8 - Predict            #
#########################################


if args['step']<=8:
	os.chdir('classifiers')
	logger.info('Running: predict_all.py')
	if (args['predpath']):
		# Purpose:
		# Performs the final interface residue
		# prediction using pre-defined classifiers
		# and the compiled feature sets generated
		# previously in compile_features.py.
		#
		# Files Affected:
		# - [PREDPATH]/[UNIPROTID1]_[UNIPROTID2].pkl                                    # For each UniProt in interactions
		#
		# Expected Outcome:
		# Final prediction pickels should be
		# generated for each interaction in the
		# inter_file.
		#
		# Known Bugs:
		# - Updated predictions are never made
		#   even if the feature availability
		#   has changed
		run_with_logging(['python', 'predict_all.py', uniprot_file, inter_file, args['predpath']], logger)
	else:
		# Purpose:
		# Same as above...
		#
		# Files Affected:
		# - /home/adr66/eclair/classifiers/predictions/[UNIPROTID1]_[UNIPROTID2].pkl    # For each UniProt in interactions
		#
		# Expected Outcome:
		# Same as above...
		#
		# Known Bugs:
		# Same as above...
		run_with_logging(['python', 'predict_all.py', uniprot_file, inter_file], logger)
	os.chdir('..')
	
	logger.info('Step 8 - Predict: Complete)')

logger.info('------------- Finished running get_features.py ----------------')
