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
import subprocess
import shutil
from collections import defaultdict
from itertools import combinations

from .utils import check_file
from .assembly import Assembly
from .alignment import TypingResults, MacSyHit
from .recovery import UnionFind, Recovery, is_connected, validate_cluster
from .output import generate_comp_output, generate_output
from .get_map import get_map
from .log import log


# ─────────────────────────────────────────────
# Core methods
# ─────────────────────────────────────────────
def prokkanno(assembly_info: Assembly, temp_dir: Path, threads: int, verbose: bool = False):
    '''
    Annotate the fasta sequence by prokka.
    '''
    cmd = ['prokka', assembly_info.sample_path.name, '--force', '--cpus', str(threads), '--cdsrnaolap', '--prefix', 
           assembly_info.sample_name, '--locustag', assembly_info.temp_id, '-o', 
           str(temp_dir / assembly_info.sample_name)]
    subprocess.run(cmd, check = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

    assembly_info.gff_file = temp_dir / assembly_info.sample_name / f'{assembly_info.sample_name}.gff'
    assembly_info.gbk_file = temp_dir / assembly_info.sample_name / f'{assembly_info.sample_name}.gbk'
    assembly_info.faa_file = temp_dir / assembly_info.sample_name / f'{assembly_info.sample_name}.faa'
    assembly_info.fna_file = temp_dir / assembly_info.sample_name / f'{assembly_info.sample_name}.ffn'

    log('Prokka annotation completed.', verbose = verbose)

def ICE_filter(ICE_res: Path) -> list[str]:
    '''
    Among results, there might be some incomplete IMEs, or some IMEs are part of ICE, those will be removed
    from the results, and finally return a list with the name of all candidates MGEs (the name contains the 
    mark of the MGE type).
    '''
    # Dictionaries to track IME genes and positions
    IME_gen_dict = {}    # {sys_id: [gene_name1, gene_name2, ...]}
    IME_pos_dict = {}    # {sys_id: [gene_hit1, gene_hit2, ...]}
    ICE_pos_dict = {}    # {sys_id: [gene_hit1, ...]}

    # Lists to track ICE and AICE systems
    ICE_list = []
    AICE_list = []

    # Open the ICE result file
    with open(ICE_res) as file:
        for line in file:
            if 'Chromosome' not in line:   # XML results are located in ./ICEscan/Chromosome
                continue

            lines = line.strip().split('\t')
            
            # Filter out loner genes not belonging to a cluster
            if lines[7] != '1':   # locus num != 1 indicates loner genes, e.g., phage integrase
                continue

            hit_id, gene_name = lines[1], lines[2]  # lines[2] is gene_name
            sys_id = lines[5]                        # system ID, e.g., aHPS7_T4SS_typeG_2, aHPS7_IME_1

            if 'IME' in sys_id:
                IME_gen_dict.setdefault(sys_id, []).append(gene_name)
                IME_pos_dict.setdefault(sys_id, []).append(hit_id)
            elif 'AICE' in sys_id:
                if sys_id not in AICE_list:
                    AICE_list.append(sys_id)
            else:  # Real ICE, e.g., T4SS-containing system
                ICE_pos_dict.setdefault(sys_id, []).append(hit_id)
                if sys_id not in ICE_list:
                    ICE_list.append(sys_id)

    # Keywords to detect mob and integrase genes in IME
    mob_keywords = ('Relaxase_', 'T4SS_MOB')
    int_keywords = {
        'Phage_integrase', 'UPF0236', 'Recombinase', 'rve',
        'TIGR02224', 'TIGR02249', 'TIGR02225', 'PB001819'
    }

    # Filter out incomplete IMEs that lack both mob and integrase
    if IME_gen_dict:
        for pre_IME, genes in IME_gen_dict.items():
            has_mob = any(any(k in g for k in mob_keywords) for g in genes)
            has_int = any(g in int_keywords for g in genes)
            if not (has_mob and has_int):
                IME_pos_dict.pop(pre_IME)

    # Determine IMEs that are not subregions of any ICE
    IME_list = []
    if IME_pos_dict:
        for pre_IME, pos_hits in IME_pos_dict.items():
            is_subregion = any(set(pos_hits).issubset(set(ice_hits))
                               for ice_hits in ICE_pos_dict.values())
            if not is_subregion:
                IME_list.append(pre_IME)

    # Return a combined list of all ICE, IME, and AICE system IDs
    return ICE_list + IME_list + AICE_list



def get_MGE(assembly_info: Assembly, rootdir: Path, temp_dir: Path, verbose: bool = False) -> dict[str, TypingResults]: 
    '''
    Detect all MGE (Mobile Genetic Elements) in the assembly using MacSyFinder.
    '''
    assembly_info.ICE_dir = temp_dir / f'{assembly_info.sample_name}_MacsyFinder'
    ICE_res = assembly_info.ICE_dir / 'all_systems.tsv'

    # Remove the ICE output folder if it already exists
    if assembly_info.ICE_dir.exists():
        shutil.rmtree(assembly_info.ICE_dir)
    
    assembly_info.ICE_dir.parent.mkdir(parents = True, exist_ok = True)

    # Run MacSyFinder to detect ICE systems
    model_dir = rootdir / 'data' / 'macsydata'
    ICE_cmd = ['macsyfinder', '--db-type', 'ordered_replicon', '--models-dir', str(model_dir), 
               '--models', 'ICEscan', 'all', '--replicon-topology', 'linear', '--coverage-profile', 
               '0.3', '--sequence-db', assembly_info.faa_file, '-o', str(assembly_info.ICE_dir)]
    
    subprocess.run(ICE_cmd, check = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    # ICE_res columns: replicon, hit_id, gene_name, hit_pos, model_fqn, sys_id, sys_loci, locus_num
    ftag = ICE_filter(ICE_res)  # List containing sys_id of all predicted ICE, AICE, and IME

    MGE_results = {}

    if ftag:
        with open(ICE_res) as file:
            for line in file:
                if 'Chromosome' not in line:
                    continue
                content = line.strip().split('\t')
    
                # Filter out genes not belonging to a cluster
                if content[7] != '1':  # locus num != 1 means loner genes like phage integrase
                    continue
                if content[5] not in ftag:  # sys_id must be included in ftag
                    continue
                
                # Determine category based on sys_id
                if 'IME' in content[5]:
                    category = 'IME'
                elif 'AICE' in content[5]:
                    category = 'AICE'
                else:
                    category = 'ICE'

                # Store gene tags
                if content[5] not in MGE_results:
                    MGE_results[content[5]] = TypingResults(content[5], category)
                MGE_results[content[5]].add_hit(content[2], content[4], content[1], assembly_info.contigs)
    
    # Sort each value in MGE_results by element.locus_num
    for v in MGE_results.values():
        v.elements.sort(key = lambda i: i.locus_num)

    log('The primary MGE detection by MacSyFinder is completed.', verbose = verbose)
    return MGE_results


def check_recoveries(recovery_MGEs: dict[str, list[MacSyHit]]) -> dict[str, list[MacSyHit]]:
    '''
    Check the recovery candidates and filter out those that do not meet the criteria.
    '''
    checked_recovery_MGEs = {}
    counters = {}
    for MGE_type, MGE_hits in recovery_MGEs.items():
        if not MGE_hits:
            continue

        n = len(MGE_hits)
        uf = UnionFind(n)
        
        # build graph
        for i in range(n):
            for j in range(i + 1, n):
                if is_connected(MGE_hits[i], MGE_hits[j], 10 if 'T4SS' not in MGE_type else 30):
                    uf.union(i, j)

        # group clusters
        clusters = defaultdict(list)
        for i in range(n):
            root = uf.find(i)
            clusters[root].append(MGE_hits[i])
        if MGE_type == 'IME':
            # validate clusters
            for cluster in clusters.values():
                if validate_cluster(cluster, MGE_type, 2):
                    counters[MGE_type] = counters.get(MGE_type, 0) + 1
                    checked_recovery_MGEs[f'{MGE_type}_{counters[MGE_type]}'] = cluster
        elif MGE_type == 'AICE':
            # validate clusters
            for cluster in clusters.values():
                if validate_cluster(cluster, MGE_type, 3):
                    counters[MGE_type] = counters.get(MGE_type, 0) + 1
                    checked_recovery_MGEs[f'{MGE_type}_{counters[MGE_type]}'] = cluster
        else:
            # validate clusters
            for cluster in clusters.values():
                if validate_cluster(cluster, MGE_type, 8):
                    counters[MGE_type] = counters.get(MGE_type, 0) + 1
                    checked_recovery_MGEs[f'{MGE_type}_{counters[MGE_type]}'] = cluster

    return checked_recovery_MGEs


def ICE_recovery(assembly_info: Assembly, verbose: bool = False) -> dict[str, list[MacSyHit]]:
    '''
    Recover ICEs from rejected candidates.
    '''
    # For draft genomes, MacSyFinder may miss some ICEs. This function recovers ICEs 
    # from rejected results.
    ICE_rej = assembly_info.ICE_dir / 'rejected_candidates.tsv'

    recovery_MGEs = {}    # MGE type : MacSyHit
    # Only 2 genes are required for an IME, so recovery is usually not necessary, but included here.

    # Parse rejected candidates
    with open(ICE_rej, 'rt') as file:
        for line in file:
            if 'Chromosome' not in line:  # XML results are located in ./ICEscan/Chromosome
                continue

            lines = line.strip().split('\t')
            MGE_type, hit_id, gene_name = lines[2].rsplit('/', 1)[-1], lines[4], lines[6]  
            # MGE_type e.g., ICEscan/Chromosome/T4SS_typeB; lines[6] is gene_name

            if MGE_type not in recovery_MGEs:
                recovery_MGEs[MGE_type] = []
            recovery_MGEs[MGE_type].append(MacSyHit.from_macsy_line(hit_id, gene_name, assembly_info.contigs))

    if not recovery_MGEs:
        return {}

    checked_recovery_MGEs = check_recoveries(recovery_MGEs)
    # Sort elements inside each object of the three lists
    for results in checked_recovery_MGEs.values():
        results.sort(key = lambda i: i.locus_num)

    log('Collected candidates for ICE recovery.', verbose = verbose)
    return checked_recovery_MGEs


def remove_duplicate(assembly_info: Assembly, MGE_results: dict, recovery_MGEs: dict, verbose: bool = False
                     ) -> tuple[dict[str, TypingResults], dict[str, list[MacSyHit]]]:
    '''
    Remove duplicate MGE results from the recovery list.
    '''
    # IME is very short, so it do not necessary need to be recovered, however, we will still check for duplicates
    # But in the rest of the script, we will only consider ICEs for further analysis
    new_recovery_MGEs = {}

    for mge_name, results in recovery_MGEs.items():
        if 'IME' in mge_name:
            if any('ICE' in key for key in MGE_results) or any('IME' in key for key in MGE_results):
                mge_locus_sets = [{element.locus_num for element in tr.elements}
                                  for tr in MGE_results.values()]
                if not any({hit.locus_num for hit in results}.issubset(locus_set)
                           for locus_set in mge_locus_sets):
                    new_recovery_MGEs[mge_name] = results
            else:
                new_recovery_MGEs[mge_name] = results
        elif 'AICE' in mge_name:
            if any('AICE' in key for key in MGE_results):
                mge_locus_sets = [{element.locus_num for element in tr.elements}
                                  for tr in MGE_results.values()]
                if not any({hit.locus_num for hit in results}.issubset(locus_set)
                           for locus_set in mge_locus_sets):
                    new_recovery_MGEs[mge_name] = results
            else:
                new_recovery_MGEs[mge_name] = results
        else:
            # The rest is ICEs, and the name should be 'T4SS_typeB' or similar
            if any('ICE' in key for key in MGE_results):
                mge_locus_sets = [{element.locus_num for element in tr.elements}
                                for tr in MGE_results.values()]
                if not any({hit.locus_num for hit in results}.issubset(locus_set)
                        for locus_set in mge_locus_sets):
                    new_recovery_MGEs[mge_name] = results
            else:
                new_recovery_MGEs[mge_name] = results

    if any('T4SS' in key for key in recovery_MGEs) and (
        any('IME' in key for key in MGE_results)):

        recovery_locus_sets = [{hit.locus_num for hit in ice_result} 
                               for ice_name, ice_result in recovery_MGEs.items() if 'T4SS' in ice_name]

        remove_keys = [
            sys_id for sys_id, tr in MGE_results.items()
            if any(
                {element.locus_num for element in tr.elements}.issubset(recovery_set)
                for recovery_set in recovery_locus_sets
            )
        ]
        for key in remove_keys:
            MGE_results.pop(key)

    remove_keys = set()
    for MGE_name, MGE_result in MGE_results.items():
        hit_id = set([element.locus_num for element in MGE_result.elements])
        for contig_id, contig in assembly_info.contigs.items():
            if contig.plasmid != 'T' and hit_id.issubset({gene.locus_num for gene in contig.genes.values()}):
                MGE_result.contig = contig_id
                break
        if MGE_result.category == 'ICE' and not MGE_result.contig:
            for mpf in MGE_result.mpf:
                new_tag = f'T4SS_{next(iter(mpf))}'
                if new_tag not in new_recovery_MGEs:
                    new_recovery_MGEs[new_tag] = []
                for i in MGE_result.elements:
                    new_recovery_MGEs[new_tag].append(i)
                remove_keys.add(MGE_name)
        elif not MGE_result.contig:
            remove_keys.add(MGE_name)
    for key in remove_keys:
            MGE_results.pop(key)

    remove_keys = set()
    for (MGE_name1, MGE_result1), (MGE_name2, MGE_result2) in combinations(MGE_results.items(), 2):

        if MGE_name1 in remove_keys or MGE_name2 in remove_keys:
            continue
        loci1 = {element.locus_num for element in MGE_result1.elements}
        loci2 = {element.locus_num for element in MGE_result2.elements}

        if loci1 == loci2:
            MGE_result1.mob.update(MGE_result2.mob)
            MGE_result1.mpf.update(MGE_result2.mpf)

            remove_keys.add(MGE_name2)

    for key in remove_keys:
        MGE_results.pop(key)

    log('Duplicate MGE results removed.', verbose = verbose)
    return MGE_results, new_recovery_MGEs

def MGE_reorder(MGE_results: dict[str, TypingResults], verbose: bool = False) -> dict[str, TypingResults]:
    '''
    Reorder MGE results based on their boundary positions.
    '''
    remove_keys = []

    for name1, MGE1 in MGE_results.items():
        for name2, MGE2 in MGE_results.items():
            if name2 == name1:
                continue
            if MGE2.DR1 <= MGE1.DR1 <= MGE2.DR4 and MGE2.DR1 <= MGE1.DR4 <= MGE2.DR4:
                remove_keys.append(name1)

    for k in remove_keys:
        MGE_results.pop(k)

    new_MGE_results = {}
    counters = {}

    for sys_id, MGE_result in MGE_results.items():
        counters[MGE_result.category] = counters.get(MGE_result.category, 0) + 1          # count
        new_key = f'{MGE_result.category}_{counters[MGE_result.category]}'                # IME_1、ICE_2
        MGE_result.sys_id = new_key
        new_MGE_results[new_key] = MGE_result

    log('MGE results reordered.', verbose = verbose)
    return new_MGE_results

def boundary_of_rICE(assembly_info: Assembly, recovery_MGEs: dict[str, list[MacSyHit]], temp_dir: Path, 
                     verbose: bool = False, EXPAND = 20) -> dict[str, Recovery]:

    recovery_ICEs = {}   # eg. {T4SS_typeG : Recovery}
    recount = 0
    for mge_name, mge_result in recovery_MGEs.items():
        if not 'T4SS' in mge_name:
            continue
        recount += 1
        recovery_ICE = Recovery().parse_raw_recovery(mge_name, mge_result, assembly_info, EXPAND)
        
        if recovery_ICE:
            recovery_ICE.merge_tRNA(assembly_info, temp_dir, EXPAND)
            recovery_ICEs[f'rICE_{recount}'] = recovery_ICE

    for key in list(recovery_ICEs.keys()):
        if not recovery_ICEs[key]:
            recovery_ICEs.pop(key)
    
    remove_keys = set()

    for (rICE_name1, rICE_result1), (rICE_name2, rICE_result2) in combinations(recovery_ICEs.items(), 2):

        if rICE_name1 in remove_keys or rICE_name2 in remove_keys:
            continue

        genes1 = rICE_result1.all_gene_num
        genes2 = rICE_result2.all_gene_num

        overlap = len(genes1 & genes2)
        overlap_ratio = overlap / min(len(genes1), len(genes2))
        # sometimes the edge will have one or two unmatched gene, so we cannot just use issubset to pending if they are 
        # the same ICE
        if overlap_ratio >= 0.8:

            # keep the larger one
            if len(genes1) >= len(genes2):
                rICE_result1.mob.update(rICE_result2.mob)
                rICE_result1.mpf.update(rICE_result2.mpf)

                remove_keys.add(rICE_name2)

            else:
                rICE_result2.mob.update(rICE_result1.mob)
                rICE_result2.mpf.update(rICE_result1.mpf)

                remove_keys.add(rICE_name1)

    for key in remove_keys:
        recovery_ICEs.pop(key, None)
    
    log('Recovery ICE results merged with tRNA information.', verbose = verbose)
    return recovery_ICEs    

def assembly_prediction(input_file: Path, json_set: bool, threads: int, rootdir: Path, output_dir: Path, temp_dir: Path, 
            verbose: bool):
    '''
    The main function for single genome analysis.
    '''
    # Despite I already changed the ICE distance from 20 to 40, but the inner pending expand still be 20 to aviod noise.
    # ── Step 1: Annotate the assembly and get gene information ────────────────────────────────
    assembly_info = Assembly(input_file, input_file.stem, check_file(input_file), 'TMPID')
    # treat fasta and genbank differerntly to collect the gene details
    if assembly_info.file_type == 'fasta':
        prokkanno(assembly_info, temp_dir, threads, verbose)
        assembly_info.get_gene_info_from_gff(verbose)
    else:
        assembly_info.get_gene_info_from_gbk(temp_dir, verbose)
    assembly_info.information_filter_marker(temp_dir, threads, verbose)
    # ── Step 2: Detect MGEs and recover ICEs ───────────────────────────────────────────────
    MGE_results = get_MGE(assembly_info, rootdir, temp_dir, verbose)
    recovery_MGEs = ICE_recovery(assembly_info, verbose)
    MGE_results, recovery_MGEs = remove_duplicate(assembly_info, MGE_results, recovery_MGEs, verbose)
    # ── Step 3: Reorder MGEs and determine the boundary of rICEs ─────────────────────────────
    for MGE_result in MGE_results.values():
        MGE_result.merge_tRNA(assembly_info, temp_dir, EXPAND = 20)
    log('MGE results merged with tRNA information.', verbose = verbose)
    MGE_results = MGE_reorder(MGE_results, verbose)
    recovery_ICEs = boundary_of_rICE(assembly_info, recovery_MGEs, temp_dir, verbose, EXPAND = 20)
    # ── Step 4: Generate output files ───────────────────────────────────────────────
    if MGE_results:
        for MGE_name, MGE_result in MGE_results.items():
            MGE_result.basic_info(assembly_info, rootdir, temp_dir, threads)
        if json_set:
            final_dir, js_dir = get_map(assembly_info, MGE_results, rootdir, temp_dir, output_dir, threads, verbose)
            generate_comp_output(assembly_info, MGE_results, final_dir, js_dir, rootdir, verbose)
    generate_output(assembly_info, MGE_results, temp_dir, recovery_ICEs, rootdir, output_dir, threads, verbose)