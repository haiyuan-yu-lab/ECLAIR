from mjm_parsers import unzip_res_range, parse_dictionary_list
import glob, os

'''Write known interface residues to a file. 1=known, 0=not known, nan=unobserved residue'''

interactomes_dir = '../interactomes'
uniprot_info_file = '../uniprot/uniprot_info.txt'
interface_ppi_file = 'ires_perppi_alltax.txt'

min_coverage = 0.50
output_file = '/home/mjm659/ires_ml/features/per_feature/ires50.txt'



###

interactomes = set()
for f in glob.glob(os.path.join(interactomes_dir, '*_binary_hq.txt')):
	print 'parsing', f, '...'
	interactomes.update(  set([tuple(sorted(l.strip().split('\t')[:2])) for l in open(f)][1:])  )

###

uniprot2info = dict((e['id'], e) for e in parse_dictionary_list(uniprot_info_file))

output = open(output_file, 'w')
for e in parse_dictionary_list(interface_ppi_file):

	uniprotA = e['UniProtA']
	uniprotB = e['UniProtB']
	
	assert uniprotA <= uniprotB  #should be sorted already
	
	if uniprotA not in uniprot2info or uniprotB not in uniprot2info:
		continue
	
	if (uniprotA, uniprotB) not in interactomes:
		continue
	
	iresA = set(int(r) for r in unzip_res_range(e['UniProtIresA']))
	iresB = set(int(r) for r in unzip_res_range(e['UniProtIresB']))
	upresA = set(int(r) for r in unzip_res_range(e['UniProtAllResA']))
	upresB = set(int(r) for r in unzip_res_range(e['UniProtAllResB']))
	
	resA_string = []
	for res in range(1, int(uniprot2info[uniprotA]['length'])+1):
		if res in iresA:
			resA_string.append('1')
		elif res in upresA:
			resA_string.append('0')
		else:
			resA_string.append('nan')
	
	resB_string = []
	for res in range(1, int(uniprot2info[uniprotB]['length'])+1):
		if res in iresB:
			resB_string.append('1')
		elif res in upresB:
			resB_string.append('0')
		else:
			resB_string.append('nan')
	
	if resA_string.count('nan')/float(len(resA_string)) > 1.0 - min_coverage:
		continue
	if resB_string.count('nan')/float(len(resB_string)) > 1.0 - min_coverage:
		continue
	
	output.write('%s\t%s\t%s\t%s\n' %(uniprotA, uniprotB, ';'.join(resA_string), ';'.join(resB_string)))

output.close()
	
