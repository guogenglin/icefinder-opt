# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 11:36:36 2025

@author: Genglin Guo
@e-mail: 2019207025.njau.edu.cn

This is an optimized version of ICEfinder2. I have rewritten some scripts and updated the logic.
The original developers are Meng Wang and Hong-Yu OU from the School of Life Sciences & Biotechnology, 
Shanghai Jiao Tong University. You can find the original website and their contact information here: 
https://tool2-mml.sjtu.edu.cn/ICEberg3/ICEfinder.php
"""

import pathlib
import argparse
import subprocess
import sys
import shutil
import multiprocessing
from Bio import SeqIO
# this warning ignore set is for biopython, the higher version of biopython may output some warning, 
# but it dosen't affect the prediction.
import warnings
from Bio import BiopythonDeprecationWarning
warnings.simplefilter("ignore", BiopythonDeprecationWarning)

__version__ = '1.0'

workdir = pathlib.Path.cwd()
rootdir = pathlib.Path(__file__).resolve().parent
tmp_dir = workdir / 'tmp'
fa_dir = tmp_dir / 'fasta'
gb_dir = tmp_dir / 'gbk'

def check_dependencies():
    # Checks dependencies are available
    dependencies = ['kraken', 'defense-finder', 'blastp', 'blastn', 'prodigal', 'prokka', 
                    'macsyfinder', 'hmmsearch', 'vmatch', 'abricate']
    for i in dependencies:
        try:
            subprocess.check_call(['which', i], stdout = subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print(f"Error: could not find {i} tool", file = sys.stderr)
            sys.exit(1)

def get_argument():
    # Parsers
    parser = argparse.ArgumentParser(description = 'icefinder-opt', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_group_1 = parser.add_argument_group('Input and Output')
    parser_group_2 = parser.add_argument_group('Parameters')

    # Input and output
    parser_group_1.add_argument('-i', '--input', required = True, nargs = '+', type = str, 
                                help = 'FASTA/Genbank format file, Genbank format file accepted only for single genome.')
    parser_group_1.add_argument('-t', '--type', required = False, type = str, default = 'single', 
                                help = 'Genome Type: Single/Metagenome')
    parser_group_1.add_argument('-o', '--output', required = False, type = str, default = 'ICEfinder_result',
                              help = 'Output dir')

    # Parameters
    parser_group_2.add_argument('-c', '--threads', required = False, type = int, default = min(multiprocessing.cpu_count(), 4), 
                        help = 'Threads to use for BLAST searches')
    parser_group_2.add_argument('-j', '--json', action='store_true', help = 'output the json based genetic map')
    parser_group_2.add_argument('-v', '--version', action = 'version', version = 'ICEfinder-optimal v' + __version__, 
                        help = 'Show version number and exit')
    
    return parser

def pending_file_type(filename):
# check the file type is, only fasta and gbk is available
	filetype = ''
    # try to open it by SeqIO
	with open(filename, 'rt') as handle1:
		fasta = SeqIO.parse(handle1, 'fasta')
		if any(fasta):
			filetype = 'fa'
	with open(filename, 'rt') as handle2:
		gbk = SeqIO.parse(handle2, 'gb')
		if any(gbk):
			filetype = 'gb'
	return filetype

def main():
    # Initialize
    args = get_argument().parse_args()
    
    print(
        'The original authors of ICEfinder2 are Meng Wang and Hong-Yu Ou.\n'
        'For more information, please visit: https://tool2-mml.sjtu.edu.cn/ICEberg3/ICEfinder.php\n\n'
        'This "optimal" version is a modified version based on their script.\n'
        'For detailed information about the ICEfinder2, please contact the original authors.\n'
        'If you have any questions and suggestions regarding this modified version, I would be very happy to hear.\n'
        'E-mail:2019207025@njau.edu.cn\n'
    )
    
    print(
        'Please note that ICE predictions in draft genomes may be biased.\n'
        'Experimental validation or additional sequencing is recommended to confirm the results.\n'
    )
    
    check_dependencies()

    intype = args.type
    input_files = args.input
    
    # create some work dir
    tmp_dir.mkdir(exist_ok = True)
    gb_dir.mkdir(exist_ok = True)
    fa_dir.mkdir(exist_ok = True)
    pathlib.Path(args.output).mkdir(exist_ok = True)
    
    for input_file in input_files:
        input_path = pathlib.Path(input_file)
        infile_name = input_path.stem
    
        filetype = pending_file_type(input_file)
    
        if intype == 'single':
            # input_file is the file, filetype is fa(fasta) or gb(genbank)
            from script.single import _single
            _single(infile_name, input_file, filetype, rootdir, args.output, args.json, args.threads)
        elif intype == 'metagenome':
            from script.metaICE import _meta
            _meta(infile_name, input_file, args.threads)
        else:
            print('Error: unrecongnized file type.')
            sys.exit(1)
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    if pathlib.Path(workdir / 'defense-finder-tmp').exists():
        shutil.rmtree(workdir / 'defense-finder-tmp')

if __name__ == '__main__':
    main()
