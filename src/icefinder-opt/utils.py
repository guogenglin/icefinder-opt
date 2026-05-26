# -*- coding: utf-8 -*-
"""
Created on Thu May  7 16:10:06 2026

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

from __future__ import annotations
from typing import TYPE_CHECKING

import os
import shutil
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

from typing import Any
from .log import log, ICEfinderError
if TYPE_CHECKING:
    from .alignment import TypingResults
    from .assembly import Assembly, Gene


def check_cpus(cpus: Any = None, max_cpus: int = 32, verbose: bool = False) -> int:
    '''
    Check the number of CPUs to use for parallel processing. If cpus is None, use all available.
    '''
    avail_cpus = os.cpu_count() or max_cpus
    if isinstance(cpus, str):
        cpus = int(cpus) if cpus.isdigit() else avail_cpus
    elif isinstance(cpus, float):
        cpus = int(cpus)
    else:
        cpus = avail_cpus
    cpus = min(cpus, avail_cpus, max_cpus)
    log(f'Using {cpus=}', verbose)
    return cpus


def check_file(input_file: Path):
    '''
    check the file type, only fasta and gbk are available
    '''

    try:
        fasta = SeqIO.parse(input_file, 'fasta')
        next(fasta)
        return 'fasta'
    except Exception:
        pass

    try:
        gbk = SeqIO.parse(input_file, 'gb')
        next(gbk)
        return 'gbk'
    except Exception:
        pass

    raise ICEfinderError(
        'Wrong input file type, only "fasta" or "genbank" format are accepted.'
    )


def check_programs_shutil(progs: list[str], verbose: bool = False):
    '''
    Check if programs are installed and executable using shutil.which.
    '''
    missing = []
    for prog in progs:
        path = shutil.which(prog)
        if path:
            print(f'{prog}: {path}')
        else:
            missing.append(prog)
    if not missing:
        log(f'All required programs found: {", ".join(progs)}', verbose = verbose)
    if missing:
        raise ICEfinderError(f'Missing programs: {", ".join(missing)}')


def format_sequence(sequence: str) -> str:
    '''
    Format a sequence into lines of 60 characters.
    '''
    return '\n'.join(sequence[i:i+60] for i in range(0, len(sequence), 60))


def gene_key(prefix: str, idx: int) -> str: return f'{prefix}_{str(idx).zfill(5)}'

def get_DR(sample_path: Path, temp_dir: Path, verbose: bool = False):
    '''
    Get the direct repeats (DRs) from the assembly.
    '''
    DRindex = temp_dir / f'{sample_path.stem}_DR'
    # fasta and gbk can both be the input of mkvtree
    # Build k-mer index for the input sequence
    mkvtree_cmd = ['mkvtree', '-db', str(sample_path), '-indexname', str(DRindex), '-dna', '-pl', '-allout']
    subprocess.run(mkvtree_cmd, check = True)

    vmatch_cmd = ['vmatch', '-l', '15', str(DRindex)]
    result = subprocess.run(vmatch_cmd, capture_output = True, text = True, check = True)
    out = result.stdout
    
    # DRout : 22    0 1228868   D    22    0 1736928   0    8.87e-02     44   100.00
    # length_left, seq_num_i, post_left, strand, length_right, seq_num_j, pos_right, distance, E-value, score, identity
    DRlist = []
    for line in out.strip().split('\n'):
        if not line.startswith('#'):
            lines = line.strip().split()
            if not lines:
                continue
            
            left_start  = int(lines[2])
            left_end    = int(lines[2]) + int(lines[0])
            right_start = int(lines[6])
            right_end   = int(lines[6]) + int(lines[4])
            
            # Filter DRs that are too long or too short
            if (right_end - left_start) > 500000 or (right_end - left_start) < 5000:
                continue
            DRlist.append(f'{left_start},{left_end},{right_start},{right_end}')

    log('Direct repeats (DRs) extracted.', verbose = verbose)
    return DRlist 

def getfa(fasta_file: Path, contig: str, start: int, end: int) -> Seq:
    '''Get the sequence from the fasta file based on the contig, start and end positions.'''
    sequence = Seq('')
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id == contig:
            sequence = record.seq[int(start) - 1 : int(end)]
            break
    
    return sequence

def process_gene(gene_id: str, assembly_info: Assembly, MGE_result: TypingResults, feature_default = 'Flank', 
                 is_ICE = False) -> Gene:
    '''
    Process a gene and return a Gene object with the appropriate feature and product information.
    '''
    gene = assembly_info.contigs[MGE_result.contig].genes[gene_id]
    feature = feature_default

    if is_ICE and min([hit.locus_num for hit in MGE_result.elements]) <= gene.locus_num <= max([hit.locus_num for hit in MGE_result.elements]) :   # if the gene within ICE
        if gene_id in [hit.hit_id for hit in MGE_result.elements]:
            gene_hit = next((hit for hit in MGE_result.elements if hit.hit_id == gene_id))
            feature, prot4 = gene_hit.tags.split('@')   # tags = {category}@{tag}
        else:
            feature,prot4 = '', ''   # if the gene within ICE but not belong to T4SS
        if prot4:   # if it is belong to T4SS
            if gene.product == 'hypothetical protein':
                gene.product = prot4
            else:
                gene.product = f'{gene.product}, {prot4}'

    gene.get_args(assembly_info, feature)

    return gene

def oritseq(sample_path: Path, regi: str, contig: str, start: int, end: int, rootdir: Path, 
            temp_dir: Path, threads: int) -> str:
    '''
    Get the oriT sequence of the predicted ICE by blasting against the oriT database.
    '''
    oritseq = '-'
    fafile = temp_dir / f'{regi}_fororit.fa'   # whole length sequence of ICE
    with open(fafile,'w') as orif:
        seq = getfa(sample_path, contig, start, end)
        orif.write('>fororit\n')
        orif.write(str(seq))
        
    oriT_Database = rootdir / 'data' / 'oriT_db'
    # search oriT
    blast_cmd = ['blastn', '-db', str(oriT_Database), '-query', str(fafile), '-evalue', '0.01', 
                 '-num_threads', str(threads), '-word_size', '11', '-outfmt', '6 std slen stitle', 
                 '-num_alignments', '1']

    process = subprocess.run(blast_cmd, capture_output = True, text = True, check = True)
    
    for line in process.stdout.strip().split('\n'):
        lines = line.strip().split('\t')
        if lines[0]:   # if there is a matched oriT
            coverage = int(lines[3]) / int(lines[12])
            identity = float(lines[2]) / 100
            havalue = coverage * identity
            if havalue > 0.49:
                oritseq = getfa(fafile, 'fororit', int(lines[6]) - 1, int(lines[7]))
                # lines[6] matched start, lines[7] matched end
                break
    return str(oritseq)

def abricate(faa_file: Path, db: str) -> tuple[dict[str, str], list[str]]:
    '''
    Run abricate to predict AMR genes and return a dictionary of gene_id to AMR type and a list of gene_ids.
    '''
    command = ['abricate', '-db', db, str(faa_file)]
    process = subprocess.run(command, capture_output = True, text = True, check = True)

    blast_filter = {}
    gene_list = []
    for line in process.stdout.strip().split('\n'):
        if not line.startswith('#'):
            lines = line.strip().split('\t')
            blast_filter[lines[1]] = lines[5]
            gene_list.append(lines[5].rsplit('_', 1)[0])
    return blast_filter, gene_list