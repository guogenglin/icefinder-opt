# -*- coding: utf-8 -*-
"""
Created on Fri May 23 17:23:58 2026

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

import json
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .utils import abricate
from .assembly import Assembly
from .alignment import TypingResults
from .recovery import Recovery
from .log import log


def get_color(region: str) -> str:
    '''Get the color for the region type in the visualization.'''
    coldict = {
        'T4SS-type ICE' : 'fill:rgba(0, 128, 164,0.9)',
        'IME':'fill:rgba(41,76,166,0.9)', 
        'AICE':'fill:rgba(255, 84,0,0.9)'
        }
    
    return coldict[region]

def generate_comp_output(assembly_info: Assembly, MGE_results: dict[str, TypingResults], final_dir: Path, 
                         js_dir: Path, rootdir: Path, verbose = False):
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

def final_refine(recovery_ICEs: dict[str, Recovery], MGE_results: dict[str, TypingResults]) -> dict[str, Recovery]:
    remove_keys = set()
    for recovery_name, recovery_ICE in recovery_ICEs.items():
        middle_contigs = [middle_contig for middle_contig in recovery_ICE.middle_contig]
        for MGE_result in MGE_results.values():
            if MGE_result.contig in middle_contigs:
                remove_keys.add(recovery_name)
                continue
    for key in remove_keys:
        recovery_ICEs.pop(key)

    return recovery_ICEs

def generate_output(assembly_info: Assembly, MGE_results: dict[str, TypingResults], temp_dir: Path, 
                    recovery_ICEs: dict[str, Recovery], rootdir: Path, output_dir: Path, threads: int, 
                    verbose: bool = False):
    '''
    Generate the output files for the predicted ICEs/IMEs/AICEs and recovery ICEs.
    '''
    final_dir = output_dir / assembly_info.sample_name
    out_detail = output_dir / 'ICE_details.tsv'
    out_summary = output_dir / 'ICE_summary.tsv'

    if not out_detail.exists():
        
        detail_header = ['Isolate', 'MGE', 'Location', 'Length', 'GC', 'Relaxase_Type', 
                         'Systems', 'oriT', 'attL', 'attR', 'att_seq', 'tRNA', 'AMR', 'middle']
        summary_header = ['Isolate', 'ICE(rICE_included)', 'ICE(rICE_excluded)', 'IME', 'AICE', 'rICE']
        
        with open(out_detail, 'w') as detail, open(out_summary, 'w') as summary:
            detail.write(('\t'.join(detail_header)))
            detail.write('\n')
            summary.write(('\t'.join(summary_header)))
            summary.write('\n')
    
    out_detail_file = open(out_detail, 'a')
    out_summary_file = open(out_summary, 'a')
    
    ICE_count, IME_count, AICE_count = 0, 0, 0
    
    if MGE_results:
        for MGE_name, MGE_result in MGE_results.items():
            if 'IME' in MGE_name:
                IME_count += 1
            elif 'AICE' in MGE_name:
                AICE_count += 1
            else:
                ICE_count += 1
                
            relaxase = ','.join(MGE_result.mob) if MGE_result.mob else '-'
            system = ','.join(MGE_result.mpf) if MGE_result.mpf else '-'
            
            attL, attR = '-', '-'
            if MGE_result.DR2:
                attL = f'{MGE_result.contig}:{MGE_result.DR1}..{MGE_result.DR2}'
                attR = f'{MGE_result.contig}:{MGE_result.DR3}..{MGE_result.DR4}'

            if not final_dir.exists():
                final_dir.mkdir(parents = True)

            outfa = final_dir / f'{MGE_result.ICE_name}.fasta'

            with open(outfa, 'w') as out_handle:
                SeqIO.write(SeqRecord(MGE_result.seq, id = MGE_result.ICE_ID, description=''), out_handle, 'fasta')

            # trna eg. tRNA-Leu,contig1:1421705..1421790[+]
            _, amr_list = abricate(outfa, 'resfinder')
            ICE_amr = (';').join(amr_list) if amr_list else '-'
            middle_contig = '-'
            
            detail_content = [assembly_info.sample_name, MGE_result.ICE_name, 
                              f'{MGE_result.contig}:({MGE_result.DR1}..{MGE_result.DR4})', str(MGE_result.length), 
                              MGE_result.gc, relaxase, system, MGE_result.oritseqs, attL, attR, MGE_result.DR_seq, 
                              MGE_result.trnaout, ICE_amr, middle_contig]
            
            out_detail_file.write(('\t'.join(detail_content)))
            out_detail_file.write('\n')
    
    rICE_count = 0
    if recovery_ICEs:
        recovery_ICEs = final_refine(recovery_ICEs, MGE_results)
        for rICE_name, ICE_result in recovery_ICEs.items():
            # info: [DR1, DR2, DR3, DR4, final_left, final_right, end0_last, end1_last, trnalist, left_contig, right_contig]
            rICE_count += 1
            ICE_result.getrfasta(assembly_info, rICE_count, output_dir, temp_dir, rootdir, threads, verbose)
            relaxase = ','.join(ICE_result.mob) if ICE_result.mob else '-'
            system = ','.join(ICE_result.mpf) if ICE_result.mpf else '-'

            attL, attR = '-', '-'
            if ICE_result.DR2:
                attL = f'{next(iter(ICE_result.left_contig.keys()))}:{ICE_result.DR1}..{ICE_result.DR2}'
                attR = f'{next(iter(ICE_result.right_contig.keys()))}:{ICE_result.DR3}..{ICE_result.DR4}'

            if not final_dir.exists():
                final_dir.mkdir(parents = True)
            _, amr_list = abricate(final_dir / f'{ICE_result.ICE_name}.fasta', 'resfinder')
            ICE_amr = (';').join(amr_list) if amr_list else '-'        

            if not ICE_result.all_middle:
                ICE_result.all_middle = ['-']
            detail_content = [assembly_info.sample_name, ICE_result.ICE_name, (';').join(ICE_result.ICE_location), 
                              str(ICE_result.length), ICE_result.gc, relaxase, system, ICE_result.oritseqs, attL, attR, 
                              ICE_result.DR_seq, ICE_result.trnaout, ICE_amr, (';').join(ICE_result.all_middle)]
            
            out_detail_file.write(('\t'.join(detail_content)))
            out_detail_file.write('\n')
            
    summary_content = [assembly_info.sample_name, str(ICE_count + rICE_count), str(ICE_count), str(IME_count), str(AICE_count), str(rICE_count)]
        
    out_summary_file.write(('\t'.join(summary_content)))
    out_summary_file.write('\n')
    
    out_detail_file.close()
    out_summary_file.close()
    
    print(f'{assembly_info.sample_name} : {str(ICE_count + rICE_count)} ICE, {str(IME_count)} IME')

    log('Output files generated.', verbose = verbose)