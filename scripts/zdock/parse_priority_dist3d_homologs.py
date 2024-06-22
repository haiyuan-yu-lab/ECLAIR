from mjm_parsers import parse_dictionary_list, unzip_res_range
from mjm_ml_methods import normalize
from collections import defaultdict
import numpy as np
import random, glob, os, sys

'''Parse a single docking experiment per interaction in priority order: (1) PDB docking, (2) Mixed Docking, and (3) ModBase docking.'''

if len(sys.argv) > 2:
	docking_score_cutoff = int(sys.argv[2])
	ires_cutoff = 0
else:
	docking_score_cutoff = 0
	ires_cutoff = 0

pdb_dsasa_file = 'ires/dist3d_homology_docking.txt'

pdb_ires_file = 'ires/IRES_homology_docking.txt'

output_prefix = '/home/adr66/eclair/features/per_feature/zdock_dist3d_hom_PRIORITY_%icut' %(docking_score_cutoff)
summary_file = 'ires/priority_dist3d_hom_zdock_%icut_%ires.txt' %(docking_score_cutoff, ires_cutoff)

interaction2sasa = defaultdict(lambda: defaultdict(lambda: ([], []))) #(p1,p2) -> (subA, subB) -> ([dsasas1...], [dsasas2...])
pdb2uniprots = defaultdict(set)  #store all uniprots seen in each PDB

#----------------------------------------------
interactomes_dir = sys.argv[1]
cocrystal_file = '../pdb/ires_perppi_alltax.txt'
interactome3d_files = glob.glob('../interactomes/*HQ.txt')

interactomes = set()
for f in glob.glob(interactomes_dir):
    interactomes.update(  set(map(lambda e: (e['UniProt_1'],e['UniProt_2']), parse_dictionary_list(f)))  )

cc_interactions = set()
for e in parse_dictionary_list(cocrystal_file):
    cc_interactions.add((e['UniProtA'], e['UniProtB']))
    
int3d_interactions = set()
for f in interactome3d_files:
    int3d_interactions.update(  set([tuple(sorted(l.strip().split('\t')[:2])) for l in open(f)][1:])   )
   
prediction_interactions = interactomes - cc_interactions - int3d_interactions
print len(prediction_interactions)

already_computed = set()
if os.path.exists(output_prefix + '_max.txt'):
	with open(output_prefix + '_max.txt','r') as f:
		for l in f:
			already_computed.add('_'.join(sorted(l.strip().split('\t'))))
#----------------------------------------------


#----------------------------------------------
subunits2ires = {}  #(subA, subB, rank) -> (num_iresA, num_iresB)
for e in parse_dictionary_list(pdb_ires_file):

	iresA = len(unzip_res_range(e['UniProtA_IRES']))
	iresB = len(unzip_res_range(e['UniProtB_IRES']))
	
	subunits2ires[(e['SubunitA'], e['SubunitB'], e['Rank'])] = (iresA, iresB)
#----------------------------------------------
    

#First pass, collect best models for all PREDICTION interactions, find ratios of docking evidence types
for e in parse_dictionary_list(pdb_dsasa_file):
	
	interaction = (e['UniProtA'], e['UniProtB'])
	dock_pair = (e['SubunitA'], e['SubunitB'])
	zdock_score = float(e['ZDOCK_Score'])
	
	#first round, only get features for interactions to be predicted
	if interaction not in prediction_interactions or '_'.join(sorted(interaction)) in already_computed:
		continue
	
	#if we've already seen this interaction docked using a different set of subunit models (i.e. from a different source) continue
	if interaction in interaction2sasa and dock_pair not in interaction2sasa[interaction]:
		continue
	
	if zdock_score < docking_score_cutoff:
		continue
	
	if ires_cutoff > 0:
		try:
			iresA, iresB = subunits2ires[(e['SubunitA'], e['SubunitB'], e['Rank'])]
			if iresA + iresB < ires_cutoff:
				continue
		except KeyError:
			continue
	
	dsasas1 = np.array([ float(r) if r != 'nan' else np.nan for r in e['UniProtA_dSASA'].split(';') ])
	dsasas2 = np.array([ float(r) if r != 'nan' else np.nan for r in e['UniProtB_dSASA'].split(';') ])
	
	if all((dsasas1==0.0) | np.isnan(dsasas1)) or all((dsasas2==0.0) | np.isnan(dsasas2)):
		continue  #poor docking
	
	if e['UniProtA'] == e['UniProtB']:  #homodimers have the same features for both proteins, save minimum dist3d for either subunit
		both_dsasas = [dsasas1, dsasas2]
		nan_mask = np.ma.masked_array(both_dsasas, np.isnan(both_dsasas))
		dsasa_min = np.min(nan_mask, axis=0)
		dsasa_min = np.array([dsasa_min.data[r] if dsasa_min.mask[r]==False else np.nan for r in range(len(dsasa_min.data))])
		interaction2sasa[interaction][dock_pair][0].append(normalize(dsasa_min))
		interaction2sasa[interaction][dock_pair][1].append(normalize(dsasa_min))
	
	else:	
		interaction2sasa[interaction][dock_pair][0].append(normalize(dsasas1))
		interaction2sasa[interaction][dock_pair][1].append(normalize(dsasas2))


#get source fractions:
         #pdb, mixed, modbase
source_fracs = [0,0,0]
for interaction in interaction2sasa.keys():

	assert len(interaction2sasa[interaction]) == 1

	subA, subB = interaction2sasa[interaction].keys()[0]

	if len(subA) < 32 and len(subB) < 32:
		source_fracs[0] += 1  #both pdb subunits
	elif len(subA) == 32 and len(subB) == 32:
		source_fracs[2] += 1  #both modbase subunits
	else:
		source_fracs[1] += 1  #one pdb and one modbase subunit

print 'PREDICT_INTERACTIONS:', sum(source_fracs)
if sum(source_fracs) != 0:
	source_fracs = [float(t)/sum(source_fracs) for t in source_fracs]
	print 'PREDICT RATIOS:', source_fracs
print

#SECOND pass, collect ALL models for all CC interactions, then try to match ratios seen in predicted

cc_interaction2sasa = defaultdict(lambda: defaultdict(lambda: ([], []))) #(p1,p2) -> (subA, subB) -> ([dsasas1...], [dsasas2...])

for e in parse_dictionary_list(pdb_dsasa_file):
	
	interaction = (e['UniProtA'], e['UniProtB'])
	dock_pair = (e['SubunitA'], e['SubunitB'])
	zdock_score = float(e['ZDOCK_Score'])
	
	#first round, only get features for interactions to be predicted
	if interaction in prediction_interactions or '_'.join(sorted(interaction)) in already_computed:
		continue
	
	if zdock_score < docking_score_cutoff:
		continue
	
	if ires_cutoff > 0:
		try:
			iresA, iresB = subunits2ires[(e['SubunitA'], e['SubunitB'], e['Rank'])]
			if iresA + iresB < ires_cutoff:
				continue
		except KeyError:
			continue
	
	dsasas1 = np.array([ float(r) if r != 'nan' else np.nan for r in e['UniProtA_dSASA'].split(';') ])
	dsasas2 = np.array([ float(r) if r != 'nan' else np.nan for r in e['UniProtB_dSASA'].split(';') ])
	
	if all((dsasas1==0.0) | np.isnan(dsasas1)) or all((dsasas2==0.0) | np.isnan(dsasas2)):
		continue  #poor docking
	
	if e['UniProtA'] == e['UniProtB']:  #homodimers have the same features for both proteins
		both_dsasas = [dsasas1, dsasas2]
		nan_mask = np.ma.masked_array(both_dsasas, np.isnan(both_dsasas))
		dsasa_max = np.max(nan_mask, axis=0)
		dsasa_max = np.array([dsasa_max.data[r] if dsasa_max.mask[r]==False else np.nan for r in range(len(dsasa_max.data))])
		cc_interaction2sasa[interaction][dock_pair][0].append(dsasa_max)
		cc_interaction2sasa[interaction][dock_pair][1].append(dsasa_max)
	
	else:	
		cc_interaction2sasa[interaction][dock_pair][0].append(dsasas1)
		cc_interaction2sasa[interaction][dock_pair][1].append(dsasas2)

#adjust ratios
total_interactions = len(cc_interaction2sasa)
target_pdb = total_interactions * source_fracs[0]
target_mixed = total_interactions * source_fracs[1]
target_mb = total_interactions * source_fracs[2]

cc_subunits = {}

final_cc_interaction2sasa = defaultdict(lambda: defaultdict(lambda: ([], []))) #(p1,p2) -> (subA, subB) -> ([dsasas1...], [dsasas2...])
#interactions sorted by number of sources available, low to high, with an element of randomness
for interaction in sorted(cc_interaction2sasa.keys(), key=lambda i: len(cc_interaction2sasa[i].keys())+random.random()):
	#~ print interaction
	#if only one source of zdock (one pair of subunits docked), just use it
	if len(cc_interaction2sasa[interaction].keys()) == 1:
		cc_subunits[interaction] = cc_interaction2sasa[interaction].keys()[0]
		interaction2sasa[interaction] = cc_interaction2sasa[interaction]
		#~ print 'one choice:', cc_interaction2sasa[interaction].keys()[0]
		#~ print
	else:
		pdb_subunits = len([s for s in cc_subunits.values() if len(s[0]) + len(s[1]) < 16])
		mixed_subunits = len([s for s in cc_subunits.values() if len(s[0]) + len(s[1]) > 32 and len(s[0]) + len(s[1]) < 64])
		mb_subunits = len([s for s in cc_subunits.values() if len(s[0]) + len(s[1]) == 64])
		
		diffs = sorted([(max(target_pdb-pdb_subunits, 0), 'PDB'), (max(target_mixed-mixed_subunits, 0), 'MIXED'), (max(target_mb-mb_subunits, 0), 'MB')], reverse=True)
		
		#~ print diffs
		
		avail_pdb_subunits = [s for s in cc_interaction2sasa[interaction] if len(s[0]) + len(s[1]) < 16]
		avail_mixed_subunits = [s for s in cc_interaction2sasa[interaction] if len(s[0]) + len(s[1]) > 32 and len(s[0]) + len(s[1]) < 64]
		avail_mb_subunits = [s for s in cc_interaction2sasa[interaction] if len(s[0]) + len(s[1]) == 64]
		
		#~ print avail_pdb_subunits, avail_mixed_subunits, avail_mb_subunits
		
		for d in diffs:
			
			if d[1] == 'PDB' and len(avail_pdb_subunits) != 0:
				cc_subunits[interaction] = avail_pdb_subunits[0]
				interaction2sasa[interaction][avail_pdb_subunits[0]] = cc_interaction2sasa[interaction][avail_pdb_subunits[0]]
				#~ print avail_pdb_subunits[0]
				break
			elif d[1] == 'MIXED' and len(avail_mixed_subunits) != 0:
				cc_subunits[interaction] = avail_mixed_subunits[0]
				interaction2sasa[interaction][avail_mixed_subunits[0]] = cc_interaction2sasa[interaction][avail_mixed_subunits[0]]
				#~ print avail_mixed_subunits[0]
				break
			elif d[1] == 'MB' and len(avail_mb_subunits) != 0:
				cc_subunits[interaction] = avail_mb_subunits[0]
				interaction2sasa[interaction][avail_mb_subunits[0]] = cc_interaction2sasa[interaction][avail_mb_subunits[0]]
				#~ print avail_mb_subunits[0]
				break
		#~ print
					
			
#get final fractions of cc interactions:
         #pdb, mixed, modbase
cc_fracs = [0,0,0]
for interaction in cc_subunits.keys():

	assert len(interaction2sasa[interaction]) == 1

	subA, subB = interaction2sasa[interaction].keys()[0]

	if len(subA) < 32 and len(subB) < 32:
		cc_fracs[0] += 1  #both pdb subunits
	elif len(subA) == 32 and len(subB) == 32:
		cc_fracs[2] += 1  #both modbase subunits
	else:
		cc_fracs[1] += 1  #one pdb and one modbase subunit
	

print 'CC_INTERACTIONS:', sum(cc_fracs)
cc_fracs = [float(t)/sum(cc_fracs) for t in cc_fracs]
print 'CC RATIOS:', cc_fracs
print

#-----------------------------------------------------------------------------

print 'writing feature files...'

output_avg = open(output_prefix + '_avg.txt', 'w')
output_max = open(output_prefix + '_max.txt', 'w')
output_min = open(output_prefix + '_min.txt', 'w')
output_top1 = open(output_prefix + '_top1.txt', 'w')

output_summary = open(summary_file, 'w')

for i in interaction2sasa:
	p1, p2 = i
	
	p1_dsasas, p2_dsasas = [], []
	
	output_summary.write('%s\t%s\t%s\t%s\n' %(p1, p2, interaction2sasa[(p1,p2)].keys()[0][0], interaction2sasa[(p1,p2)].keys()[0][1]))
	
	for docki, dock_pair in enumerate(interaction2sasa[i]):
		
		if docki == 0:   #first pair of subunits docked together (if there are more than one pair of subunits selected for this interaction)
			output_top1.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in interaction2sasa[i][dock_pair][0][0]]), ';'.join(['%.2f' %(res) for res in interaction2sasa[i][dock_pair][1][0]])))
		
		p1_dsasas += interaction2sasa[i][dock_pair][0]
		p2_dsasas += interaction2sasa[i][dock_pair][1]
	
	#~ print i, len(p1_dsasas), len(p2_dsasas), len(p1_dsasas[0]), len(p2_dsasas[0])
	
	mdat = np.ma.masked_array(p1_dsasas, np.isnan(p1_dsasas))
	p1_means = np.mean(mdat, axis=0)
	p1_means = [p1_means.data[r] if p1_means.mask[r]==False else np.nan for r in range(len(p1_means.data))]
	p1_max = np.max(mdat, axis=0)
	p1_max = [p1_max.data[r] if p1_max.mask[r]==False else np.nan for r in range(len(p1_max.data))]
	p1_min = np.min(mdat, axis=0)
	p1_min = [p1_min.data[r] if p1_min.mask[r]==False else np.nan for r in range(len(p1_min.data))]
	
	mdat = np.ma.masked_array(p2_dsasas, np.isnan(p2_dsasas))
	p2_means = np.mean(mdat, axis=0)
	p2_means = [p2_means.data[r] if p2_means.mask[r]==False else np.nan for r in range(len(p2_means.data))]
	p2_max = np.max(mdat, axis=0)
	p2_max = [p2_max.data[r] if p2_max.mask[r]==False else np.nan for r in range(len(p2_max.data))]
	p2_min = np.min(mdat, axis=0)
	p2_min = [p2_min.data[r] if p2_min.mask[r]==False else np.nan for r in range(len(p2_min.data))]
	
	output_avg.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in p1_means]), ';'.join(['%.2f' %(res) for res in p2_means])))
	output_max.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in p1_max]), ';'.join(['%.2f' %(res) for res in p2_max])))
	output_min.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in p1_min]), ';'.join(['%.2f' %(res) for res in p2_min])))
	
output_avg.close()
output_max.close()
output_min.close()
output_top1.close()
output_summary.close()
