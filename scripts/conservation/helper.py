import subprocess

def parse_dictionary_list(filename, delim='\t'):
	'''Uses header row in file as keys for dictionaries'''
	handle = open(filename, 'r')
	header_keys = handle.readline().strip().split(delim)
	
	data = []
	for l in handle:
		line = l.replace('\n', '').split(delim)
		cur_dict = {}
		for i, key in enumerate(header_keys):
			cur_dict[key] = line[i]
		data.append(cur_dict)
	
	return data


def set_from_ires_list(ires_list):
	ires_set = set()
	for elt in ires_list:
		try:
			ires_set.add(int(elt))
		except ValueError:
			try:
				[beg,end] = elt.split('-')
				ires_set.update(range(int(beg), int(end)+1))
			except ValueError:
				continue
	return ires_set


def calc_int3d_ires(structure_file, SEQ_BEGIN1, SEQ_BEGIN2):
	
	PDB_BEGIN1 = None; PDB_END1 = None; PDB_BEGIN2 = None; PDB_END2 = None
	for l in open(structure_file):
		
		if l[:3] not in ['ATO', 'TER']:
			continue
		index = int(l[22:26])
		
		if PDB_BEGIN1 == None and l[:4] == 'ATOM':
			PDB_BEGIN1 = index
			continue
		if PDB_END1 != None and PDB_BEGIN2 == None:
			PDB_BEGIN2 = index
			continue
		if l[:3] == 'TER' and PDB_END1 == None:
			PDB_END1 = index
			continue
		if l[:3] == 'TER' and PDB_END1 != None:
			PDB_END2 = index
			continue

	int_res, err = subprocess.Popen(['python', 'irescalc.py', structure_file], stdout=subprocess.PIPE).communicate()
	
	ires1, ires2 = int_res.split('\n')[:2]
	
	ires1_pdb = [int(r) for r in ires1.split(',') if r != '']
	ires1_unp = [r - PDB_BEGIN1 + SEQ_BEGIN1 for r in ires1_pdb]

	ires2_pdb = [int(r) for r in ires2.split(',') if r != '']
	ires2_unp = [r - PDB_BEGIN2 + SEQ_BEGIN2 for r in ires2_pdb]
	
	pdb1_range = '[%i-%i]' %(PDB_BEGIN1, PDB_END1)
	pdb2_range = '[%i-%i]' %(PDB_BEGIN2, PDB_END2)
	
	return ires1_unp, ires2_unp, pdb1_range, pdb2_range


def sorted_pdb_for_uniprot(pdb_list):
	def num_mappable_res(entry):
		ranges = entry['MappableResInPDBChainOnUniprotBasis'][1:-1].split(',')
		count = 0
		for range in ranges:
			if '-' in range:
				[start, end] = range.split('-')
				count += int(end) - int(start) + 1
			else:
				count += 1
		return count

	def entry_to_tuple(entry):
		return (entry['PDB'], entry['Chain'], str(num_mappable_res(entry)))

	sorted_list = sorted(pdb_list, key=num_mappable_res, reverse=True)

	return map(entry_to_tuple, sorted_list)


def best_pdb_without_homolog(pcs1, homologs2, pdb_uni_mapping):
	# pcs1 must be sorted by best pdb
	for [pdb, chain] in pcs1:
		other_unis = []
		for c in pdb_uni_mapping[pdb]:
			if c != chain:
				other_unis.append(pdb_uni_mapping[pdb][c])
		homolog_set = set(homologs2)
		homolog_set.intersection_update(other_unis)

		if len(homolog_set) == 0:
			return pdb, chain
	return None, None
