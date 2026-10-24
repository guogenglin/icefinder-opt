# -*- coding: utf-8 -*-
"""
Created on Fri May 23 17:23:58 2026

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from Bio import SeqIO

from .utils import write_matching_gbk
from .assembly import Assembly
from .alignment import TypingResults
from .recovery import Recovery
from .log import log, ICEfinderError
from .annotation import resistance_names


DETAIL_HEADER = ['Isolate', 'MGE', 'Location', 'Length', 'GC', 'Relaxase_Type', 'Systems',
                 'oriT', 'attL', 'attR', 'att_seq', 'tRNA', 'AMR', 'metal', 'middle']
SUMMARY_HEADER = ['Isolate', 'ICE(rICE_included)', 'ICE(rICE_excluded)', 'IME', 'AICE', 'rICE']


def get_color(region: str) -> str:
    '''
    Get the color for the region type in the visualization.
    '''
    coldict = {
        'T4SS-type ICE' : 'fill:rgba(0, 128, 164,0.9)',
        'IME':'fill:rgba(41,76,166,0.9)', 
        'AICE':'fill:rgba(255, 84,0,0.9)'
        }
    
    return coldict[region]

def generate_comp_output(assembly_info: Assembly, MGE_results: dict[str, TypingResults], final_dir: Path, 
                         js_dir: Path, rootdir: Path, verbose: bool = False) -> None:
    '''
    Generate the output files for the visualization of the predicted ICEs/IMEs/AICEs.
    '''
    jsback = rootdir / 'data' / 'js'
    
    i = 1 
    ICEsumlist = []
    homelist = []
    if MGE_results:
        for MGE_name, MGE_result in MGE_results.items():
            regi = assembly_info.sample_name + '_' + MGE_result.sys_id   # eg. aHPS7_ICE2
            if 'IME' in regi:
                typeIE = 'IME'
            elif 'AICE' in regi:
                typeIE = 'AICE'
            else:
                typeIE = 'T4SS-type ICE'
            ICEs = {
                'region': 'Region'+ str(i),
                'location': f'{str(MGE_result.DR1)}..{str(MGE_result.DR4)}',
                'length': str(MGE_result.length),
                'gc': MGE_result.gc,
                'type': typeIE,
                'detail': regi
            }
            homedict = {
                'start' : MGE_result.DR1,
                'end' : MGE_result.DR4,
                'color' : get_color(typeIE),
                'info' : '_'.join([str(MGE_result.length), MGE_result.gc, typeIE]),
                'text' : 'Region'+ str(i)
               }
            homelist.append(homedict)
            ICEsumlist.append(ICEs)
            i += 1
    
    ICEsum = final_dir / f'{assembly_info.sample_name}_ICEsum.json'
    with open(ICEsum,'w') as ice_file:
        json.dump(ICEsumlist, ice_file, indent = 4)
    
    shutil.copytree(jsback, js_dir, dirs_exist_ok = True)

def ensure_output_table(path: Path, columns: list[str]) -> None:
    '''
    Initialize a table or reject an incompatible header before appending sample results.
    '''
    header = '\t'.join(columns) + '\n'
    if path.exists() and path.stat().st_size:
        with open(path) as handle:
            if handle.readline() != header:
                raise ICEfinderError(f'Incompatible output table header: {path}. Use a fresh output directory.')
    else:
        with open(path, 'w') as handle:
            handle.write(header)


def detail_row(assembly_info: Assembly, result: TypingResults | Recovery, location: str,
               att_left: str, att_right: str, amr: str, metal: str, middle: str = '-') -> list[str]:
    '''
    Serialize the shared annotation fields of a single-contig MGE or a recovered ICE.
    '''
    return [assembly_info.sample_name, result.ICE_name, location, str(result.length), result.gc,
            ','.join(sorted(result.mob)) or '-', ','.join(sorted(result.mpf)) or '-', result.oritseqs,
            att_left, att_right, result.DR_seq, result.trnaout, amr, metal, middle]


def export_prediction(result: TypingResults | Recovery, final_dir: Path) -> None:
    '''
    Write already prepared annotation records through the same FASTA and GenBank validation path.
    '''
    if not result.output_ready or not result.output_records:
        raise ICEfinderError(f'Output records are not ready for {result.ICE_name or "unnamed prediction"}.')
    final_dir.mkdir(parents = True, exist_ok = True)
    outfa = final_dir / f'{result.ICE_name}.fasta'
    with open(outfa, 'w') as handle:
        SeqIO.write(result.output_records, handle, 'fasta')
    write_matching_gbk(result.output_records, outfa)

def generate_output(assembly_info: Assembly, MGE_results: dict[str, TypingResults],
                    recovery_ICEs: dict[str, Recovery], output_dir: Path,
                    verbose: bool = False) -> None:
    '''
    Generate the output files for the predicted ICEs/IMEs/AICEs and recovery ICEs.
    '''
    assembly_info.require_resistance()
    for result in [*MGE_results.values(), *recovery_ICEs.values()]:
        if not result.output_ready:
            raise ICEfinderError('Prepare all prediction records before generating output tables.')
    output_dir.mkdir(parents = True, exist_ok = True)
    final_dir = output_dir / assembly_info.sample_name
    out_detail = output_dir / 'ICE_details.tsv'
    out_summary = output_dir / 'ICE_summary.tsv'
    ensure_output_table(out_detail, DETAIL_HEADER)
    ensure_output_table(out_summary, SUMMARY_HEADER)
    detail_rows: list[list[str]] = []
    
    ICE_count, IME_count, AICE_count = 0, 0, 0
    
    if MGE_results:
        for MGE_name, MGE_result in MGE_results.items():
            if MGE_result.category == 'IME':
                IME_count += 1
            elif MGE_result.category == 'AICE':
                AICE_count += 1
            else:
                ICE_count += 1
                
            attL, attR = '-', '-'
            if MGE_result.DR2:
                attL = f'{MGE_result.contig}:{MGE_result.DR1}..{MGE_result.DR2}'
                attR = f'{MGE_result.contig}:{MGE_result.DR3}..{MGE_result.DR4}'

            export_prediction(MGE_result, final_dir)

            # trna eg. tRNA-Leu,contig1:1421705..1421790[+]
            ICE_amr, ICE_metal = resistance_names(assembly_info.resistance_hits, MGE_result.output_regions)
            detail_rows.append(detail_row(assembly_info, MGE_result,
                f'{MGE_result.contig}:({MGE_result.DR1}..{MGE_result.DR4})', attL, attR, ICE_amr, ICE_metal))
    
    rICE_count = 0
    if recovery_ICEs:
        for rICE_name, ICE_result in recovery_ICEs.items():
            # info: [DR1, DR2, DR3, DR4, final_left, final_right, end0_last, end1_last, trnalist, left_contig, right_contig]
            rICE_count += 1
            export_prediction(ICE_result, final_dir)
            attL, attR = '-', '-'
            if ICE_result.DR2:
                attL = f'{next(iter(ICE_result.left_contig.keys()))}:{ICE_result.DR1}..{ICE_result.DR2}'
                attR = f'{next(iter(ICE_result.right_contig.keys()))}:{ICE_result.DR3}..{ICE_result.DR4}'

            if not final_dir.exists():
                final_dir.mkdir(parents = True)
            ICE_amr, ICE_metal = resistance_names(assembly_info.resistance_hits, ICE_result.output_regions)

            detail_rows.append(detail_row(assembly_info, ICE_result, ';'.join(ICE_result.ICE_location),
                attL, attR, ICE_amr, ICE_metal, ';'.join(ICE_result.all_middle) or '-'))
            
    summary_content = [assembly_info.sample_name, str(ICE_count + rICE_count), str(ICE_count), str(IME_count), str(AICE_count), str(rICE_count)]
        
    # Append rows only after every exported region in this sample has been generated successfully.
    with open(out_detail, 'a') as detail, open(out_summary, 'a') as summary:
        detail.writelines('\t'.join(row) + '\n' for row in detail_rows)
        summary.write('\t'.join(summary_content) + '\n')
    
    print(f'{assembly_info.sample_name} : {str(ICE_count + rICE_count)} ICE, {str(IME_count)} IME')

    log('Output files generated.', verbose = verbose)
