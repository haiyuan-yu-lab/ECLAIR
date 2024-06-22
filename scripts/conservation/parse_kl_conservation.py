import os, glob

kl_files = sorted(glob.glob('sca_cons/*.kl'))

uniprot2len = dict([(l.split('\t')[0], int(l.strip().split('\t')[1])) for l in open('../interactome/uniprot_human2length.txt')])

output = open('../../features/per_feature/KLconserv.txt', 'w')

for f in kl_files:
	uniprot = os.path.basename(f).split('.')[0]
	cons = [l.strip() for l in open(f)]
	
	if uniprot not in uniprot2len:
		print uniprot, 'not a valid UniProt ID'
		continue
	
	if len(cons) != uniprot2len[uniprot]:
		print uniprot, 'invalid number of residues'
		continue
	
	output.write('%s\t%s\n' %(uniprot, ';'.join(cons)))

output.close()
