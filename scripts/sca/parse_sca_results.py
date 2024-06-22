# Authors:
# - Presumably written by some combination of Michael, Juan
#   Jay, and Aaron
# - Comments / Modifications by Shayne

# Purpose:
# This script parses the SCA output generated by
# psiblast_sca.py to calculate the final six SCA
# features used in Eclair. The script skips
# recalculating these features for any interactions
# that have previously had them calculated.

# Expected Outcomes:
# All interactions included in the input file should
# have the six SCA features calculated and should have
# corresponding entries in the relevant per_feature
# output file. The six features are...
#
# 1) SCA_PB_max - The maximum SCA score for each
#                 position in each protein
#
# 2) SCA_PB_max_norm - SCA_PB_max normalized based
#                      on the average SCA score
#                      for the entire interaction
#
# 3) SCA_PB_max_rank - A ranking of all positions
#                      in each protein based on how
#                      their SCA_PB_max compared to
#                      all other positions in the protein
#
# 4) SCA_PB_max_rank_percentile - SCA_PB_max_rank
#                                 converted to a percentile
#
# 5) SCA_PD_top10 - The average of the top 10 SCA
#                   score for each position in each
#                   protein
#
# 6) SCA_PD_all - The average SCA score for each position
#                 in each protein

# Known Bugs:
# - No known bugs, but concerns including...
# - Warnings that should be thrown on edge cases
# - Uncertainty on if it is safe to sort the UniProt
#   IDs that make up each interaction in the final
#   feature outputs

# Imports
import os, glob, sys
from numpy import mean
from mjm_parsers import parse_dictionary_list

# Find all of the SCA files / sort them based on modification date
sca_files = sorted(glob.glob('sca_corr_psiblast/*_*.sca'))
sca_files.sort(key=os.path.getmtime, reverse=True)

# Obtain list of all interactions in input files
interactome_files = glob.glob(sys.argv[1])
inters = []
for f in interactome_files:
	for l in open(f):
		inters.append(sorted(tuple(l.strip().split('\t')[:2])))


# Fixed output directories where SCA files may be found
# My understanding is that sca_dir1 is from an earlier stage
# of development. All new Eclair predictions should be in
# sca_dir2
sca_dir1 = '/mnt/data-rs-wihy299nas/mjm659/SCA/sca_corr_psiblast'
sca_dir2 = 'sca_corr_psiblast'

# Identify all SCA files corresponding to input interactions
sca_files = []
for p1, p2 in inters:
	# Search both output directories for the interaction
	f1 = glob.glob('%s/%s_%s.sca' % (sca_dir1, p1, p2))
	f2 = glob.glob('%s/%s_%s.sca' % (sca_dir2, p1, p2))
	
	# Append any results that are found
	if len(f2) > 0:
		sca_files.append(f2)
	else:
		sca_files.append(f1)

# Flatten the list
sca_files = [item for sublist in sca_files for item in sublist]

# Set up input file based on input parameter
# If no input is provided, default to fixed
# input file
if len(sys.argv) > 1:
	uniprot_info_file = sys.argv[2]
else:
	uniprot_info_file = '../uniprot/uniprot_info.txt'

# Fixed output directory
output_dir = '/home/adr66/eclair/features/per_feature/'

# Names for all of the output files, matched to their
# feature name
out_files = {'max': 'SCA_PB_max.txt',
			 'max_norm': 'SCA_PB_max_norm.txt',
			 'max_rank': 'SCA_PB_max_rank.txt',
			 'max_rank_perc': 'SCA_PB_max_rank_percentile.txt',
			 'top10': 'SCA_PB_top10.txt',
			 'all': 'SCA_PB_all.txt' }

# Generate dictionary mapping UniProt ID to sequence
uniprot2aaseq = dict((e['id'], e['sequence']) for e in parse_dictionary_list(uniprot_info_file))

#-----------------------------------------------------------------

#create and clear files for later appending
#for f in out_files.values():
#	output = open(os.path.join(output_dir, f), 'w')
#	output.close()

# Generate list of all interactions that have already
# had these features calculated (specifically based on
# the SCA_PB_max.txt file)
already_computed = set(['_'.join(l.strip().split('\t')[:2]) for l in open(os.path.join(output_dir,out_files['max']))])

#-----------------------------------------------------------------

# Break SCA files into batches of 200
batch_size = 200
batches = [sca_files[i*batch_size:(i+1)*batch_size] for i in range(0,len(sca_files)/batch_size+1)]

# Iterate over every batch
for batch_num, sca_batch in enumerate(batches):
	
	# List of all interactions?
	all_interactions = []
	
	# Map interaction to mean of DI value in SCA?
	meandi = {}
	
	# Iterate over each of the individual SCA files in the batch
	for num, f in enumerate(sca_batch):
		print num+1+batch_num*batch_size,'/', len(sca_files), f
		
		# Obtain UniProt IDs for the interaction
		u1, u2 = os.path.basename(f)[:-4].split('_')
		
		# Skip homodimers
		if u1==u2:
			continue
		
		# Skip interactions that have already been calculated
		if '_'.join([u1,u2]) in already_computed:
			continue
		
		# Skip interactions that do not show up in the UniProt
		# Info File (How is this possible?)
		# This should throw a warning
		if u1 not in uniprot2aaseq or u2 not in uniprot2aaseq:
			continue
		
		# Obtain the length of each member of the interaction
		len1, len2 = len(uniprot2aaseq[u1]), len(uniprot2aaseq[u2])
		
		# Read the SCA file
		data = open(f).readlines()
		
		# Ensure that there is consistency between the lengths from the
		# SCA alignments and the cummulative length of the individual
		# UniProt sequences
		#
		# Apparently there are 8/923 interactions where this was the case
		# in some test case
		#
		# This should definitely throw a warning
		tot_len = int(data[-1].strip().split()[1])
		if tot_len != len1 + len2:
			print 'Protein lengths unexpected. Skipping ', f
			continue
		
		# Create emtpy dictionaries to store a mapping from (UniProt, pos)
		# to SCA score. We will later use this to generate aggregate features.
		#
		# e.g {P1 : {1 : SCA(1), 2 : SCA(2), ...}, P2 : {...}}
		#
		res2sca = {u1: dict([(i,[]) for i in range(len(uniprot2aaseq[u1]))]),
				 u2: dict([(i,[]) for i in range(len(uniprot2aaseq[u2]))]) }
		
		# Keep track of all of the DIs
		all_dis = []
		
		# Iterate over the DCA files
		for l in data:
			# Obtain the two residue positions in question + the score
			res1, res2, di = [float(t) for t in l.strip().split('\t')]
			
			# Ensure that we only keep information for the inter-protein
			# relationship, and skip all intra-protein position pairs
			if res2 <= len1 or res1 > len1:
				continue
			
			# Convert from 1 index-ed list to 0 indexed positions for python
			res1 = int(res1-1)
			res2 = int(res2-1) - len1
			
			# Add the SCA scores to the proper dictionary
			res2sca[u1][res1].append(di)
			res2sca[u2][res2].append(di)
			all_dis.append(di)
		
		# Reverse the lists of SCA values so that the
		# most correlated values appear first
		for res in res2sca[u1]:
			res2sca[u1][res].sort(reverse=True)
		for res in res2sca[u2]:
			res2sca[u2][res].sort(reverse=True)
		
		# Save the mean SCA score for the accross the
		# entire interaction
		meandi[tuple(sorted([u1,u2]))] = mean(all_dis)
		
		# Save parsed SCA score lists
		all_interactions.append(res2sca)
	
	
	
	#---------------write output--------------
	print '... writing batch output for %i entries' %(len(sca_batch))
	
	# Now we write to all of the output feature files.
	# This bit is fairly straightforward so I have not
	# explicitly commented each of these steps. The format
	# for all of the feature files is as follows...
	#
	# HEAD:
	#
	# P1		P2			Feature Per Res (P1)						Feature Per Res (P2)
	# F4K567 O64852	2.056926;0.538821;0.269124;1.101526	0.492727;0.099668;0.092347;0.091123
	#
	
	# MAX
	# Write output for max SCA score for each residue
	output = open(os.path.join(output_dir, out_files['max']), 'a')
	
	for i in all_interactions:
		u1, u2 = sorted(i.keys())
		output.write('%s\t%s\t%s\t%s\n' %(u1, u2, ';'.join(['%.6f' %(i[u1][res][0]) for res in i[u1]]), ';'.join(['%.6f' %(i[u2][res][0]) for res in i[u2]])))
	
	output.close()
	
	# MAX_NORM
	# Write output for max SCA score for each residue normalized
	# by average SCA score for the interaction
	output = open(os.path.join(output_dir, out_files['max_norm']), 'a')
	
	for i in all_interactions:
		u1, u2 = sorted(i.keys())
		output.write('%s\t%s\t%s\t%s\n' %(u1, u2, ';'.join(['%.6f' %(i[u1][res][0]/meandi[(u1,u2)]) for res in i[u1]]), ';'.join(['%.6f' %(i[u2][res][0]/meandi[(u1,u2)]) for res in i[u2]])))
	
	output.close()
	
	# MAX_RANK
	# Write output for ranked max SCA score for each residue
	# i.e. we rank the residues based on their max SCA score
	# relative to the rest of the protein
	output = open(os.path.join(output_dir, out_files['max_rank']), 'a')
	
	for i in all_interactions:
		u1, u2 = sorted(i.keys())
		
		rank1 = sorted([i[u1][r][0] for r in i[u1]], reverse=True)
		rank2 = sorted([i[u2][r][0] for r in i[u2]], reverse=True)
		
		output.write('%s\t%s\t%s\t%s\n' %(u1, u2, ';'.join(['%i' %(rank1.index(i[u1][res][0])) for res in i[u1]]), ';'.join(['%i' %(rank2.index(i[u2][res][0])) for res in i[u2]])))
	
	output.close()
	
	# MAX_RANK_PERC
	# Write output for ranked percentile max SCA score for each residue
	# i.e. we just take the previously calculated ranks and turn
	# them into percentiles
	output = open(os.path.join(output_dir, out_files['max_rank_perc']), 'a')
	
	for i in all_interactions:
		u1, u2 = sorted(i.keys())
		
		rank1 = sorted([i[u1][r][0] for r in i[u1]], reverse=True)
		rank2 = sorted([i[u2][r][0] for r in i[u2]], reverse=True)
		
		output.write('%s\t%s\t%s\t%s\n' %(u1, u2, ';'.join(['%.6f' %(rank1.index(i[u1][res][0])/float(len(i[u1]))) for res in i[u1]]), ';'.join(['%.6f' %(rank2.index(i[u2][res][0])/float(len(i[u2]))) for res in i[u2]])))
	
	output.close()
	
	# TOP10
	# Write ouput for the mean of the top 10 highest SCA scores
	# for each residue
	output = open(os.path.join(output_dir, out_files['top10']), 'a')
	
	for i in all_interactions:
		u1, u2 = sorted(i.keys())
		output.write('%s\t%s\t%s\t%s\n' %(u1, u2, ';'.join(['%.6f' %(mean(i[u1][res][0:10])) for res in i[u1]]), ';'.join(['%.6f' %(mean(i[u2][res][0:10])) for res in i[u2]])))
	
	output.close()
	
	# ALL
	# Write output for the mean SCA score per residue
	output = open(os.path.join(output_dir, out_files['all']), 'a')
	
	for i in all_interactions:
		u1, u2 = sorted(i.keys())
		output.write('%s\t%s\t%s\t%s\n' %(u1, u2, ';'.join(['%.6f' %(mean(i[u1][res])) for res in i[u1]]), ';'.join(['%.6f' %(mean(i[u2][res])) for res in i[u2]])))
	
	output.close()