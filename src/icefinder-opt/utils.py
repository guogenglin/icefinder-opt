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
import warnings
from copy import copy
from pathlib import Path
from Bio import SeqIO, BiopythonWarning
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from typing import Any
from .log import log, ICEfinderError
from .commands import run_command
if TYPE_CHECKING:
    from .alignment import TypingResults
    from .assembly import Assembly, Gene


def load_gbk_records(path: Path) -> dict[str, SeqRecord]:
    '''
    Load the annotated contigs once for a sample.
    '''
    with open(path) as handle:
        return SeqIO.to_dict(SeqIO.parse(handle, 'genbank'))


def extract_gbk_region(records: dict[str, SeqRecord], contig: str, region: slice,
                       record_id: str, reverse: bool = False) -> SeqRecord:
    '''
    Slice an annotated contig in its output orientation, retaining enclosed features.
    '''
    if contig not in records:
        raise ValueError(f'Contig {contig!r} is missing from the source GenBank file.')
    source = records[contig]
    oriented = source.reverse_complement(id = True, name = True, description = True) if reverse else source
    result = oriented[region]
    result.id = record_id
    # Use the same MGE name in FASTA, LOCUS and ACCESSION.
    result.name = record_id
    start, stop, _ = region.indices(len(source))
    result.description = f'Region of {contig}, oriented slice [{start}:{stop}], strand {"-" if reverse else "+"}'
    result.annotations = {
        key: source.annotations[key]
        for key in ('organism', 'taxonomy', 'source') if key in source.annotations
    }
    result.annotations.update(molecule_type = 'DNA', topology = 'linear')
    return result


def write_matching_gbk(records: list[SeqRecord], fasta_path: Path) -> None:
    '''
    Verify record order and sequences before writing the matching GenBank file.
    '''
    with open(fasta_path) as handle:
        fasta_records = list(SeqIO.parse(handle, 'fasta'))
    if len(records) != len(fasta_records):
        raise ValueError(f'GenBank/FASTA record count mismatch for {fasta_path}.')
    for annotated, fasta in zip(records, fasta_records):
        if annotated.id != fasta.id or annotated.seq != fasta.seq:
            raise ValueError(f'GenBank/FASTA mismatch for {fasta_path}: {fasta.id}.')
    # Long MGE identifiers are intentional and must match FASTA and ACCESSION.
    # Suppress only the writer's expected long-LOCUS notice; retain all other warnings.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message = (r'^Increasing length of locus line to allow long name\. '
                       r'This will result in fields that are not in usual positions\.$'),
            category = BiopythonWarning, module = r'^Bio\.SeqIO\.InsdcIO$')
        SeqIO.write(records, fasta_path.with_suffix('.gbk'), 'genbank')


def check_cpus(cpus: Any = None, max_cpus: int = 32, verbose: bool = False) -> int:
    '''
    Check the number of CPUs to use for parallel processing. If cpus is None, use all available.
    '''
    avail_cpus = os.cpu_count() or max_cpus
    if isinstance(cpus, str):
        cpus = int(cpus) if cpus.isdigit() else avail_cpus
    elif isinstance(cpus, (int, float)):
        cpus = int(cpus)
    else:
        cpus = avail_cpus
    if cpus < 1:
        raise ICEfinderError('The number of threads must be at least 1.')
    cpus = min(cpus, avail_cpus, max_cpus)
    log(f'Using {cpus=}', verbose)
    return cpus


def check_file(input_file: Path) -> str:
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


def check_programs_shutil(progs: list[str], verbose: bool = False) -> None:
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
    print('\n')
    if not missing:
        log(f'All required programs found: {", ".join(progs)}', verbose = verbose)
    if missing:
        raise ICEfinderError(f'Missing programs: {", ".join(missing)}')


def format_sequence(sequence: str) -> str:
    '''
    Format a sequence into lines of 60 characters.
    '''
    return '\n'.join(sequence[i:i+60] for i in range(0, len(sequence), 60))


def gene_key(prefix: str, idx: int) -> str:
    '''
    Build an internal gene identifier with a five-digit locus number.
    '''
    return f'{prefix}_{str(idx).zfill(5)}'

def get_DR(sample_path: Path, temp_dir: Path, verbose: bool = False) -> list[str]:
    '''
    Get the direct repeats (DRs) from the assembly.
    '''
    DRindex = temp_dir / f'{sample_path.stem}_DR'
    # fasta and gbk can both be the input of mkvtree
    # Build k-mer index for the input sequence
    mkvtree_cmd = ['mkvtree', '-db', str(sample_path), '-indexname', str(DRindex), '-dna', '-pl', '-allout']
    run_command(mkvtree_cmd, capture_output = False)

    vmatch_cmd = ['vmatch', '-l', '15', str(DRindex)]
    result = run_command(vmatch_cmd)
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
    '''
    Get the sequence from the fasta file based on the contig, start and end positions.
    '''
    sequence = Seq('')
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id == contig:
                sequence = record.seq[int(start) - 1 : int(end)]
                break
    
    return sequence

def process_gene(gene_id: str, assembly_info: Assembly, MGE_result: TypingResults, feature_default: str = 'Flank', 
                 is_ICE: bool = False) -> Gene:
    '''
    Process a gene and return a Gene object with the appropriate feature and product information.
    '''
    # Rendering a map must not change the original annotation or a later map's description.
    gene = copy(assembly_info.contigs[MGE_result.contig].genes[gene_id])
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
            temp_dir: Path, threads: int, *, hits: list[dict[str, int]] | None = None,
            segment_lengths: list[int] | None = None) -> str:
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

    process = run_command(blast_cmd)
    
    boundaries = []
    offset = 0
    for length in segment_lengths or [len(seq)]:
        boundaries.append((offset, offset + length))
        offset += length
    for line in process.stdout.strip().split('\n'):
        lines = line.strip().split('\t')
        if lines[0]:   # if there is a matched oriT
            coverage = int(lines[3]) / int(lines[12])
            identity = float(lines[2]) / 100
            havalue = coverage * identity
            if havalue > 0.49:
                qstart, qend = int(lines[6]), int(lines[7])
                left, right = min(qstart, qend) - 1, max(qstart, qend)
                if not any(a <= left < right <= b for a, b in boundaries):
                    continue
                strand = 1 if (qend - qstart) * (int(lines[9]) - int(lines[8])) >= 0 else -1
                matched = seq[left:right]
                oritseq = matched if strand == 1 else matched.reverse_complement()
                if hits is not None:
                    hits.append({'start': left, 'end': right, 'strand': strand})
                break
    return str(oritseq)

def abricate(faa_file: Path, db: str) -> tuple[dict[str, str], list[str]]:
    '''
    Run abricate to predict AMR genes and return a dictionary of gene_id to AMR type and a list of gene_ids.
    '''
    command = ['abricate', '-db', db, str(faa_file)]
    process = run_command(command)

    blast_filter = {}
    gene_list = []
    for line in process.stdout.strip().split('\n'):
        if line.strip() and not line.startswith('#'):
            lines = line.strip().split('\t')
            blast_filter[lines[1]] = lines[5]
            gene_list.append(lines[5].rsplit('_', 1)[0])
    return blast_filter, gene_list
