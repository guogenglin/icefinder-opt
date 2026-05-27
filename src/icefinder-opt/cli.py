# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 11:36:36 2025

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn

This is an optimized version of ICEfinder2. I have rewritten some scripts and updated the logic.
The original developers are Meng Wang and Hong-Yu OU from the School of Life Sciences & Biotechnology, 
Shanghai Jiao Tong University. You can find the original website and their contact information here: 
https://tool2-mml.sjtu.edu.cn/ICEberg3/ICEfinder.php
"""

from pathlib import Path
import argparse
import sys
import tempfile

# this warning ignore set is for biopython, the higher version of biopython may output some warning, 
# but it dosen't affect the prediction.
import warnings
from Bio import BiopythonDeprecationWarning
warnings.simplefilter("ignore", BiopythonDeprecationWarning)

from .log import formatted_description, bold, ICEfinderError, log
from .utils import check_cpus, check_programs_shutil
from .single import assembly_prediction
from .metaICE import meta_prediction

__version__ = '1.1.0'

# ─────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────
def parse_args(a) -> argparse.Namespace:
    '''
    Parse command-line arguments.
    '''
    parser = argparse.ArgumentParser(prog = 'icefinder-opt', description = formatted_description, add_help = False,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    
    subparsers = parser.add_subparsers(title=bold('Command'), dest='subparser_name', metavar="")

    assembly_subparser(subparsers)
    metagenome_subparser(subparsers)

    opts = parser.add_argument_group(bold('Other options'), '')
    other_opts(opts)

    if len(a) == 0:  # No arguments, print help message
        parser.print_help(sys.stderr)
        raise ICEfinderError(f'Please specify a command; choose from {{assembly, metagenome}}')
    if any(x in a for x in {'-v', '--version'}):  # Version message
        print(__version__)
        sys.exit(0)
    if subparser := subparsers.choices.get(a[0], None):  # Check if the first argument is a subparser
        if len(a) == 1:  # Subparser help message
            subparser.print_help(sys.stderr)
            raise ICEfinderError(f'Insufficient arguments')
        if any(x in a[1:] for x in {'-h', '--help'}):  # Subparser help message
            subparser.print_help(sys.stderr)
            sys.exit(0)
    elif any(x in a for x in {'-h', '--help'}):  # Help message
        parser.print_help(sys.stderr)
        sys.exit(0)
    else:  # Unknown command
        parser.print_help(sys.stderr)
        raise ICEfinderError(f'Unknown command "{a[0]}"; choose from {{assembly, metagenome}}')
    
    return parser.parse_args(a)

def assembly_subparser(subparsers):
    '''
    Create the argument parser for the assembly subcommand.
    '''
    assembly_parser = subparsers.add_parser('assembly', description = formatted_description, add_help = False, 
                                            formatter_class = argparse.RawTextHelpFormatter, 
                                            help = 'Predict ICE/IME in an assembled genome file.', 
                                            usage = 'icefinder-opt assembly [options]')
    io_group = assembly_parser.add_argument_group(bold('Input and Output'), '')
    in_out_opts(io_group)
    # Parameters group
    param_group = assembly_parser.add_argument_group(bold('Parameters'), '')
    parameter_opts(param_group)
    opts = assembly_parser.add_argument_group(bold('Other options'), "")
    other_opts(opts)

def metagenome_subparser(subparsers):
    '''
    Create the argument parser for the metagenome subcommand.
    '''
    metagenome_parser = subparsers.add_parser('metagenome', description = formatted_description, add_help = False, 
                                            formatter_class = argparse.RawTextHelpFormatter, 
                                            help = 'Predict ICE/IME in a metagenome file.', 
                                            usage = 'icefinder-opt metagenome [options]')
    # Input and Output group
    io_group = metagenome_parser.add_argument_group(bold('Input and Output'), '')
    in_out_opts(io_group)
    # Parameters group
    param_group = metagenome_parser.add_argument_group(bold('Parameters'), '')
    parameter_opts(param_group)
    opts = metagenome_parser.add_argument_group(bold('Other options'), "")
    other_opts(opts)

def in_out_opts(opts: argparse.ArgumentParser):
    '''
    Add input and output options to the argument parser.
    '''
    opts.add_argument('-i', '--input', required = True, nargs = '+', type = str,
                        help = 'FASTA/Genbank format file, Genbank format file accepted only for single genome.')
    opts.add_argument('-o', '--output', type = str, default = 'ICEfinder_result',
                        help = 'The name of the output dict (default: ICEfinder_result)')


def parameter_opts(opts: argparse.ArgumentParser):
    '''
    Add parameter options to the argument parser.
    '''
    opts.add_argument('-t', '--threads', type = int, default = check_cpus(),
                        help = 'Number of threads (default: all available cpus)')
    # opts.add_argument('-s', '--max_space', type = int, help = 'Max gene space for scanning ICE by Macsyfinder(default: 30)')
    opts.add_argument('-j', '--json', action='store_true', 
                        help = 'output the json based genetic map')



def other_opts(opts: argparse._ArgumentGroup):
    '''
    Add other options to the argument parser.
    '''
    opts.add_argument('-V', '--verbose', action='store_true', 
                      help='Keep all output files and print debug messages to stderr.')
    opts.add_argument('-v', '--version', help='Show version number and exit', action='version')
    opts.add_argument('-h', '--help', help='Show this help message and exit', metavar='')


def main(argv: list[str] | None = None):
    '''
    Main method for ICEfinder-opt.
    '''
    argv = sys.argv[1:] if argv is None else argv
    args = parse_args(argv)
    
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
    
    output_dir = Path.cwd() / args.output
    output_dir.mkdir(exist_ok = True)
    
    check_programs_shutil(['defense-finder', 'blastp', 'blastn', 'prodigal', 'prokka', 'macsyfinder', 'hmmsearch', 
                           'vmatch', 'abricate', 'mob_recon']
                             + (['kraken'] if args.subparser_name == 'metagenome' else []), 
                           verbose = args.verbose)
    
    rootdir = Path(__file__).parents[0]
    
    if args.subparser_name == 'assembly':
        for input_file in args.input:
            with tempfile.TemporaryDirectory() as temp_dir:
                assembly_prediction(Path(input_file), args.json, args.threads, rootdir, output_dir, Path(temp_dir), args.verbose)
    elif args.subparser_name == 'metagenome':
        for input_file in args.input:
            with tempfile.TemporaryDirectory() as temp_dir:
                meta_prediction(Path(input_file), args.json, args.threads, rootdir, output_dir, Path(temp_dir), args.verbose)
    else:
        raise ICEfinderError('Invalid mode. Please choose either "assembly" or "metagenome".')
    
    log('ICEfinder completed successfully.', verbose = args.verbose)


if __name__ == '__main__':
    raise SystemExit(main())