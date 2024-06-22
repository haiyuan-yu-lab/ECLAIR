import glob
from mjm_parsers import parse_dictionary_list

interactomes = glob.glob('../interactomes/*_binary_hq.txt')
print interactomes

all_uniprots = set()

for f in interactomes:
	#~ if '9606' in f: continue   #already downloaded all human proteins from modbase
	
	data = parse_dictionary_list(f)
	for e in data:
		all_uniprots.add(e['Uniprot_A'])
		all_uniprots.add(e['Uniprot_B'])

output = open('all_uniprots.txt', 'w')
for u in sorted(all_uniprots):
	output.write(u+'\n')
output.close()




