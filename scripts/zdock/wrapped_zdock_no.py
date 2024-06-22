#!/usr/bin/env python
"""..."""

import sys, subprocess, os, re, argparse, glob
from tempfile import mkdtemp
from shutil import rmtree, copyfile

sys.path.append('/data/web-vhosts/instruct2_pipeline/tools')
from mjm_tools import extract_atoms

zdock_path = 'zdock/zdock3.0.2'
createpl_path = 'zdock/create.pl'
create_lig_path = 'zdock/create_lig'

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Runs ZDOCK and returns top models as PDB-formatted files.')
parser.add_argument('receptor', help='File name of receptor PDB file (.pdb or .ent), or gzipped pdb file (.gz), or 4-letter pdb ID.')
parser.add_argument('ligand', help='File name of ligand PDB file (.pdb or .ent), or gzipped pdb file (.gz), or 4-letter pdb ID.')

parser.add_argument('-RC', '-receptor_chain', help='Chain of Receptor file to dock. _ for no chain ID (modbase).', default='A', required=False)
parser.add_argument('-LC', '-ligand_chain', help='Chain of Ligand file to dock. _ for no chain ID (modbase).', default='A', required=False)

parser.add_argument('-N', '--num_models', help='Number of ZDOCK models to calculate.', type=int, default=2000, required=False)
parser.add_argument('-T', '--top_models', help='Number of top ZDOCK models to produce as PDB files.', type=int, default=10, required=False)
parser.add_argument('-S', '--randomization_seed', help='Set the seed for the randomization of the starting PDBs (default to no randomization).', type=int, default=-1, required=False)
parser.add_argument('-F', '--fix_receptor', help='Fix the receptor, preventing its rotation and switching with ligand during execution.', action='store_true', default=False)
parser.add_argument('--fix_modbase', help='Remove fields 73-80 in model file (tend to occur in modbase files, messing up docking). Replace with single element ID in field 77.', action='store_true', default=False)

parser.add_argument('-name', help='Run prefix attached to all output files. (default: a prefix will be created from the given subunits and chains)', required=False)
parser.add_argument('-v', '--verbose', help='Show steps and warnings.', action='store_true', default=False)
parser.add_argument('-o', '--outdir', help='Directory to save output.', default='.', required=False)

args = parser.parse_args()

#---------------------------------- ISOLATE CHAINS ----------------------------------

scratchDir = mkdtemp()

if args.verbose: print 'writing preliminary data to %s...' %(scratchDir)

if args.RC == '_': args.RC = ' '
if args.LC == '_': args.LC = ' '

receptor_basename = os.path.splitext(os.path.basename(args.receptor))[0].upper()
if args.RC != ' ': receptor_basename += '_' + args.RC
cleaned_receptor = os.path.join(scratchDir, '%s.R' %(receptor_basename))
extract_atoms(args.receptor, cleaned_receptor, chain_dict = {args.RC:set()}, chain_map = {args.RC: 'A'}, fix_modbase=args.fix_modbase)

ligand_basename = os.path.splitext(os.path.basename(args.ligand))[0].upper()
if args.LC != ' ': ligand_basename += '_' + args.LC
cleaned_ligand = os.path.join(scratchDir, '%s.L' %(ligand_basename))
extract_atoms(args.ligand, cleaned_ligand, chain_dict = {args.LC:set()}, chain_map = {args.LC: 'B'}, fix_modbase=args.fix_modbase)

if args.name == None:
	run_id = '%s--%s--ZDOCK' %(receptor_basename, ligand_basename)
else:
	run_id = args.name

#---------------------------------- ZDOCK ----------------------------------

cur_dir = os.getcwd()
os.chdir(scratchDir)

zdock_output_file = os.path.join(scratchDir, 'zdock.out')

zdock_params = [zdock_path, '-R', cleaned_receptor, '-L', cleaned_ligand, '-N', '%i' %(args.num_models)]

if args.fix_receptor:
	zdock_params += ['-F']

if args.randomization_seed != -1:
	zdock_params += ['-S', '%i' %(args.randomization_seed)]

zdock_params += ['-o',  zdock_output_file]

if args.verbose: print ' '.join(zdock_params)

os.chdir(cur_dir)
out, err = subprocess.Popen(zdock_params, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

#---------------------------------- MODEL CREATION ----------------------------------

copyfile(createpl_path, os.path.join(scratchDir, 'create.pl'))
os.system('chmod 777 %s' %(os.path.join(scratchDir, 'create.pl')))

copyfile(create_lig_path, os.path.join(scratchDir, 'create_lig'))
os.system('chmod 777 %s' %(os.path.join(scratchDir, 'create_lig')))

top_models_file = zdock_output_file + '.top'

os.chdir(scratchDir)
os.system('head -%i %s > %s' %(4+args.top_models, zdock_output_file, top_models_file))
os.system('./create.pl %s' %(top_models_file))
os.chdir(cur_dir)

#---------------------------------- OUTPUT FILES ----------------------------------

new_zdock_outfile = os.path.join(args.outdir, '%s.out' %(run_id))
copyfile(zdock_output_file, new_zdock_outfile)

#remove tmp paths from receptor and ligand file names in zdock output file:
replace_command = "grep -rl '%s/' %s | xargs sed -i 's/%s\///g'" %(scratchDir, new_zdock_outfile, scratchDir.replace('/', '\/'))
os.system(replace_command)

copyfile(cleaned_receptor, os.path.join(args.outdir, os.path.basename(cleaned_receptor)))
copyfile(cleaned_ligand, os.path.join(args.outdir, os.path.basename(cleaned_ligand)))

complex_files = glob.glob(os.path.join(scratchDir, '*.pdb'))
for f in complex_files:
	num = os.path.basename(f).split('.')[1]
	if len(num) == 1 : num = '0' + num
	new_file = os.path.join(args.outdir, '%s-%s.pdb' %(run_id, num))
	copyfile(f, new_file)
	

#cleanup temporary directory
rmtree(scratchDir)
