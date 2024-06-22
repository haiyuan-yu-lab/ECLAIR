uniprot2len = dict([(l.split('\t')[0],int(l.strip().split('\t')[1])) for l in open('/home/resources/uniprot/parsed_files/uniprot_human2length.txt')])

#~ output = open('../../features/residue_position.txt', 'w')
#~ for u in uniprot2len:
	#~ output.write('%s\t%s\n' %(u, ';'.join(['%.4f' %(float(i)/uniprot2len[u]) for i in range(uniprot2len[u])])))
#~ output.close()

output = open('../../features/residue_terminal_distance.txt', 'w')
for u in uniprot2len:
	output.write('%s\t%s\n' %(u, ';'.join(['%i' %min(i, uniprot2len[u]-i-1) for i in range(uniprot2len[u])])))
output.close()






