import os, time, sys

output_dir = '/mnt/data-rs-wihy299nas/adr66/data/sca/dca'
#dca_path = '/home/mjm659/mjm_path/'
dca_path = '/home/adr66/eclair/data/sca/'

#interaction = 'P68135_Q28372'
interaction = 'P27348_Q9BSD8'
#interaction = 'P48061_P51679'
#interaction = 'Q4G0N4_Q9H1B7'
msa = '/home/adr66/eclair/data/sca/joined_msas/%s.msa' % interaction

if sys.argv[1] == 'matlab':
	output_file = os.path.join(output_dir,'%s_lowram.dca'%interaction)
	command = '''matlab -singleCompThread -nodisplay -nosplash -r "addpath('%s');lowram_dca('%s', '%s');exit"'''  %(dca_path, msa, output_file)
else:
	output_file = os.path.join(output_dir,'%s_python1.dca'%interaction)
	command = 'python dca.py %s %s' % (msa,output_file)

starttime = time.time()
os.system(command)
interval = time.time() - starttime
print interval
