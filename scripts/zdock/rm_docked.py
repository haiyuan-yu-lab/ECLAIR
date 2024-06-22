import os, glob
from mjm_parsers import parse_dictionary_list

#docking_set = 'modbase_docking_set.txt'
#docking_dir = 'modbase_docked_models'

docking_set = 'pdb_docking_set.txt'
docking_dir = 'pdb_docked_models'


to_dock = set()
to_dock.update(set([(e['SubA'].upper(), e['SubB'].upper()) for e in parse_dictionary_list(docking_set)]))
print len(to_dock), 'models to dock'
print list(to_dock)[:3]
to_dock.update(set([(e['SubB'].upper(), e['SubA'].upper()) for e in parse_dictionary_list(docking_set)]))

docking_files = glob.glob(os.path.join(docking_dir, '*.out'))

models = set()
for d in docking_files:
	receptor, ligand = os.path.basename(d)[:-11].split('--')
	models.add((receptor, ligand))
	
print len(models), 'models computed'
print list(models)[:3]
models_to_remove = models - to_dock
print len(models_to_remove), 'models to remove'


for m in models_to_remove:
	command = 'rm %s' %(os.path.join(docking_dir, '%s--%s--ZDOCK*' %(m[0], m[1])))
	#print command
	os.system(command)





