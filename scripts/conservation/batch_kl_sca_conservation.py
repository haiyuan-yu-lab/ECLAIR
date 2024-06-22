from mjm_tools import TimeoutCmd
from mjm_parsers import parse_fasta, write_fasta
import os, glob

sca_path = '/home/mjm659/ires_ml/data/sca/SCA5_forDist/sca5/'
msa_dir = '/home/mjm659/ires_ml/data/conservation/slimmed_msas/'
output_dir = '/home/mjm659/ires_ml/data/conservation/sca_cons/'


slimmed_MSAs = glob.glob(os.path.join(msa_dir, '*.msa'))
print "calculating KL conservation for %s uniprots" %(len(slimmed_MSAs))

for i, msa in enumerate(slimmed_MSAs):
	print '%s/%s %s' %(i, len(slimmed_MSAs), msa)
	
	protein = os.path.basename(msa).split('.')[0]
	
	sca_output = os.path.join(output_dir, '%s.kl' %(protein))
	
	sca_command = '''matlab -nodisplay -nosplash -r "addpath('%s');run_sca_cons('%s', '%s');exit"'''  %(sca_path, msa, sca_output)

	print sca_command
	
	myCmd = TimeoutCmd(sca_command, 3600)   # 3600 = 1 hour, 7200 = 2 hours, 14400 = 4 hours
	myCmd.Run()

	print myCmd.out
