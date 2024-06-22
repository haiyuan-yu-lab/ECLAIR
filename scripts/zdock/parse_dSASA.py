from mjm_parsers import parse_dictionary_list
from collections import defaultdict
import numpy as np


docking_dsasa_file = 'ires/dSASA_mb_docking.txt'
output_prefix = '../../features/per_feature/zdock_MB_all'

docking_score_cutoff = -1


interaction2sasa = defaultdict(lambda: defaultdict(lambda: ([], []))) #uniprot -> pdb -> list of chains(list(SASA values))
pdb2uniprots = defaultdict(set)  #store all uniprots seen in each PDB

print 'parsing raw dsasa file...'
for e in parse_dictionary_list(docking_dsasa_file):
	
	interaction = (e['UniProtA'], e['UniProtB'])
	dock_pair = (e['SubunitA'], e['SubunitB'])
	zdock_score = float(e['ZDOCK_Score'])
	
	if zdock_score < docking_score_cutoff:
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
		interaction2sasa[interaction][dock_pair][0].append(dsasa_max)
		interaction2sasa[interaction][dock_pair][1].append(dsasa_max)
	
	else:	
		interaction2sasa[interaction][dock_pair][0].append(dsasas1)
		interaction2sasa[interaction][dock_pair][1].append(dsasas2)



print 'writing feature files...'

output_avg = open(output_prefix + '_avg.txt', 'w')
output_max = open(output_prefix + '_max.txt', 'w')
output_top1 = open(output_prefix + '_top1.txt', 'w')

for i in interaction2sasa:
	p1, p2 = i
	
	p1_dsasas, p2_dsasas = [], []
	
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
	
	mdat = np.ma.masked_array(p2_dsasas, np.isnan(p2_dsasas))
	p2_means = np.mean(mdat, axis=0)
	p2_means = [p2_means.data[r] if p2_means.mask[r]==False else np.nan for r in range(len(p2_means.data))]
	p2_max = np.max(mdat, axis=0)
	p2_max = [p2_max.data[r] if p2_max.mask[r]==False else np.nan for r in range(len(p2_max.data))]
	
	output_avg.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in p1_means]), ';'.join(['%.2f' %(res) for res in p2_means])))
	output_max.write('%s\t%s\t%s\t%s\n' %(p1, p2, ';'.join(['%.2f' %(res) for res in p1_max]), ';'.join(['%.2f' %(res) for res in p2_max])))
	
output_avg.close()
output_max.close()
output_top1.close()
