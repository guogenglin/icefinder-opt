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
import subprocess

workdir = pathlib.Path.cwd()

tmp_dir = workdir / 'tmp'
gb_dir = tmp_dir / 'gbk'

def getdf(infile_name, rootdir, threads):
    
    locus_tag_faa = tmp_dir / infile_name / f'{infile_name}.locus_tag.faa'
    if not locus_tag_faa.exists():
        infaa = gb_dir / f'{infile_name}.faa'
    else:
        infaa = locus_tag_faa
    
    #defense-finder used the same way of macsyfinder
    model_dir = rootdir / 'data' / 'macsydata'
    defcmd = ['defense-finder', 'run', '-w', str(threads), '--models-dir', str(model_dir), str(infaa)]
    
    process = subprocess.run(defcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out = process.stdout.decode()
    
    dfdict = {}
    for line in out.strip().split('\n'):
        lines = line.strip().split('\t')
        if lines[0] and lines[0] != 'replicon':
            dfdict[lines[1]] = lines[2].replace('__',',')   # hit_id : gene_name
    
    return dfdict

def run_blast(fasta_file, database, threads, evalue = '0.0001'):

    command = ['blastp', '-query', fasta_file, '-db', database, '-evalue', evalue, '-num_threads', 
               str(threads), '-max_hsps', '1', '-num_descriptions', '1', '-num_alignments', '1', 
               '-outfmt', '6 std slen stitle']
    process = subprocess.run(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out = process.stdout.decode()
    
    return out

def havalue(threshold, out):
    # pending if the gene reached the threshold
    blast_filter = {}
    for line in out.strip().split('\n'):
        lines = line.strip().split('\t')
        coverage = int(lines[3]) / int(lines[12])
        identity = float(lines[2]) / 100
        havalue = coverage * identity
        if havalue >= float(threshold):   # lines[0] : eg. TMPID_00001
            blast_filter[lines[0]] = lines[1].split('|')[1]   # line[1] : eg. resfinder|tet(B)|2
    return blast_filter

def abricate(input_file, db):
    
    command = ['abricate', '-db', db, input_file]
    process = subprocess.run(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out = process.stdout.decode()
    blast_filter = {}
    gene_list = []
    for line in out.strip().split('\n'):
        if not line.startswith('#'):
            lines = line.strip().split('\t')
            blast_filter[lines[1]] = lines[5]
            gene_list.append(lines[5].rsplit('_', 1)[0])
    return blast_filter, gene_list

def getblast(infile_name, rootdir, threads):
    
    # define inputfile
    locus_tag_faa = tmp_dir / infile_name / f'{infile_name}.locus_tag.faa'
    if not locus_tag_faa.exists():
        infaa = pathlib.Path(gb_dir) / f'{infile_name}.faa'
    else:
        infaa = tmp_dir / infile_name / f'{infile_name}.locus_tag.faa'
    
    # a dict set the config of every blast
    blast_config = {
        'IS':    {'db': str(rootdir / 'data' / 'transposase')},
        'METAL': {'db': str(rootdir / 'data' / 'metal')},
        'POP':   {'db': str(rootdir / 'data' / 'degradation')},
        'SYM':   {'db': str(rootdir / 'data' / 'symbiosis')},
        }
    
    # run blast and get output
    for db_name, config in blast_config.items():
        config['out'] = run_blast(infaa, config['db'], threads)
    
    # pending the result and collect the match
    isdict = havalue('0.64', blast_config['IS']['out'])
    vfdict, _ = abricate(infaa, 'vfdb')
    argdict, _ = abricate(infaa, 'resfinder')
    metaldict = havalue('0.64', blast_config['METAL']['out'])
    popdict = havalue('0.64', blast_config['POP']['out'])
    symdict = havalue('0.64', blast_config['SYM']['out'])
    
    dfdict = getdf(infile_name, rootdir, threads)
    
    return argdict, vfdict, isdict, dfdict, metaldict, popdict, symdict
