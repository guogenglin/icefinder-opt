# -*- coding: utf-8 -*-
"""
Created on Fri May 22 10:34:31 2026

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

import json
from pathlib import Path
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq

from .assembly import Assembly
from .alignment import TypingResults
from .utils import gene_key, process_gene
from .log import log


def getcolor(feature: str, product: list[str]) -> tuple[str, str]:
    
    coldict = {
        'DR': 'black', 'Gene': '#C0C0C0', 'Hyp': '#DCDCDC',
        'Integrase': 'blue', 'Transposase': 'yellow',
        'T4SS': 'lightpink', 'T4CP': 'orange', 'Relaxase': 'brown',
        'AR': 'red', 'tRNA': 'black', 'Flank': 'gray', 'VF': '#ba8448',
        'Defense': '#00B050', 'Metal': '#03A89E', 'Degradation': '#640B0F',
        'Symbiosis': '#FFFFCD', 'Rep': 'black', 'Tra': 'black'
        }
    
    namedict = {
        'Hyp': 'Hypothetical protein', 'Gene': 'Other gene',
        'AR': 'Antibiotic resistance gene', 'VF': 'Virulence factor',
        'Metal': 'Metal resistance', 'Flank': 'Flank region',
        'Defense': 'Defense system', 'Transposase': 'Transposase',
        'Relaxase': 'Relaxase', 'T4CP': 'T4CP', 'T4SS': 'T4SS',
        'Integrase': 'Integrase', 'Degradation': 'Degradation',
        'Symbiosis': 'Symbiosis', 'Rep': 'Rep', 'Tra': 'Tra'
        }
    
    keyword_map = {
        'Integrase': 'Integrase',
        'T4SS': 'T4SS',
        'T4CP': 'T4CP',
        'Relaxase': 'Relaxase',
        'Rep': 'Rep',
        'Tra': 'Tra',
        'IS': 'Transposase',
        'VF': 'VF',
        'AR': 'AR',
        'Defense': 'Defense',
        'Metal': 'Metal',
        'Degradation': 'Degradation',
        'Symbiosis': 'Symbiosis'
        }
    
    for key, value in keyword_map.items():
        if key in feature:
            feature = value
            break
    else:   # only if for loop haven't be break
        if feature == 'Flank':
            feature = 'Flank'
        elif not feature:
            feature = 'Hyp' if 'hypothetical protein' in product else 'Gene'
        else:
            feature = 'Gene'
    
    return coldict[feature], namedict[feature]
        
def calculate_gc(seq: Seq, window_size = 500, step_size = 50) -> dict:

    gc_contents = []
    pos = []

    for i in range(0, len(seq) - window_size + 1, step_size):
        window = seq[i:i + window_size]
        gc_content = gc_fraction(window) * 100
        gc_contents.append(round(gc_content, 2))
        # midpoint of window
        pos.append(i + window_size // 2)

    gcdict = {
        'xData': pos,
        'datasets': [{
            'name': 'GC content',
            'data': gc_contents,
            'unit': '%',
            'type': 'line',
            'valueDecimals': 1
        }]
    }

    return gcdict

def get_map(assembly_info: Assembly, MGE_results: dict[str, TypingResults], rootdir: Path, temp_dir: Path, 
            output_dir: Path, threads: int, verbose: bool = False) -> tuple[Path, Path]:
    
    final_dir = output_dir / assembly_info.sample_name
    final_dir.mkdir(parents = True, exist_ok = True)
    js_dir = final_dir / 'js'
    js_dir.mkdir(parents = True, exist_ok = True)
    viewfile = rootdir / 'data' / 'js' / 'view.html'
    gcmap = rootdir / 'data' / 'js' / 'gcmap.js'
    
    assembly_info.getblast(rootdir, temp_dir, threads, verbose)
    
    for MGE_result in MGE_results.values():
        genelist = []
        genefile = final_dir / f'{MGE_result.ICE_name}_gene.json'
        infofile = final_dir / f'{MGE_result.ICE_name}_info.json'
        gcjson = js_dir / f'{MGE_result.sys_id}_gc.js'
        mapfile = js_dir / f'{MGE_result.sys_id}.js'
        htmlfile = final_dir / f'{MGE_result.ICE_name}.html'

        max_in_contig = max(gene.locus_num for gene in assembly_info.contigs[MGE_result.contig].genes.values())
        min_in_contig = min(gene.locus_num for gene in assembly_info.contigs[MGE_result.contig].genes.values())

        # scan 5 gene before the ICE
        for mov in range(max(MGE_result.start_point - 5, min_in_contig), MGE_result.start_point):
            gene_id = gene_key(assembly_info.temp_id, mov)
            genelist.append(process_gene(gene_id, assembly_info, MGE_result))

        # genelist = [{'gene': locus_tag, 'pos': pos, 'prod': product, 'featu': feature}]
        # scan all gene within the ICE
        for mov in range(MGE_result.start_point, MGE_result.end_point + 1):
            gene_id = gene_key(assembly_info.temp_id, mov)
            genelist.append(process_gene(gene_id, assembly_info, MGE_result, feature_default = '', is_ICE = True))
        
        # scan 5 gene after the ICE
        for mov in range(MGE_result.end_point + 1, min(MGE_result.end_point + 6, max_in_contig)):
            gene_id = gene_key(assembly_info.temp_id, mov)
            genelist.append(process_gene(gene_id, assembly_info, MGE_result))
            
        with open(genefile, 'w') as gene_file:
            for gene in genelist:
                gene_info = {
                    'gene': gene.origin_id,
                    'pos': f'{gene.start}..{gene.end} [{gene.strand}], {int(gene.end) - int(gene.start) + 1}',
                    'prod': gene.product,
                    'featu': gene.feature
                }
                gene_file.write(json.dumps(gene_info) + '\n')

        if 'IME' in MGE_result.ICE_name:
            typeIE = 'IME'
        elif 'AICE' in MGE_result.ICE_name:
            typeIE = 'AICE'
        else:
            typeIE = 'T4SS-type ICE'
            
        ICEinfo = {
            'Type' : typeIE,
            'Location (nt)' : str(MGE_result.DR1) + '..' + str(MGE_result.DR4),
            'Length (bp)' : str(MGE_result.length),
            'GC Content (%)' : MGE_result.gc,
            'oriT seq' : MGE_result.oritseqs,
            'DRs' : MGE_result.DRw,
            'Relaxase Type' : ','.join(MGE_result.mob),
            'Mating pair formation systems' : ','.join(MGE_result.mpf),
            'Close to tRNA' : MGE_result.trnaout
        }
        with open(infofile, 'w') as info_file:
            json.dump(ICEinfo, info_file, indent = 4)   # output json file

        i = 1
        mapforlist = []   # forward strand
        maprevlist = []   # reverse strand
        for gene in genelist:
            color, name = getcolor(gene.feature, gene.product_map)
            anno = {
                'start' : gene.start,
                'end' : gene.end,
                'strand' : gene.strand,
                'locus_tag' : 'M'+ str(i),
                'type' : 'others',
                'color' : color,
                'description' : 
                    f'Location: {gene.start}..{gene.end} ({gene.end - gene.start + 1} bp)<br>Type: {name}<br>Detail: {gene.product_map}'
                }

            if gene.strand == '+':   # pending if the gene located in forward or reverse strand
                mapforlist.append(anno)
            else:
                maprevlist.append(anno)
            i += 1
        
        head = 'var borders = [];\nvar tta_codons = [];\nvar orfs ='
        start = str(genelist[0].start)
        end = str(genelist[-1].end)
        
        gcdict = calculate_gc(MGE_result.seq)
        
        with open(gcmap, 'rt') as original_file:   # this is a template file to draw a gcmap
            original_content = original_file.read()
        with open(gcjson,'wt') as gein2:   # write to a file for draw gcmap for this genome
            gein2t = 'var jsonData = ' + str(gcdict)+';'
            gein2.write(gein2t)
            gein2.write(original_content)	
        
        maps = str(mapforlist)+';\nvar orfs2 ='+str(maprevlist)+';\nvar clusterf2 = { start: '+start+', end: '+ \
				  end+', idx: 1, orfs: orfs, borders: borders, tta_codons:tta_codons,\
				  label: \'\', unordered: true };\nvar clusterr2 = { start: '+ start+', end: '+ \
				  end+', idx: 2, orfs: orfs2, borders: borders, tta_codons:tta_codons,\
				  label: \'\', unordered: true };\nsvgene.drawClusters("'+MGE_result.sys_id+'", [clusterf2, clusterr2], 50, 920);'
        
        with open(mapfile,'wt') as map_file:
            map_file.write(head + maps)
            
        with open(viewfile, 'rt') as file:
            file_content = file.read()
        new_content = file_content.replace('XXXX', MGE_result.sys_id)
        with open(htmlfile, 'wt') as file:
            file.write(new_content)


    log('The genetic map of ICEs are drawn', verbose = verbose)
    return final_dir, js_dir