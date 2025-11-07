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
import json
import subprocess
import shutil
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqFeature import CompoundLocation
from script.function import getblast, abricate
from collections import defaultdict
from itertools import combinations, chain


workdir = pathlib.Path.cwd()

tmp_dir = workdir / 'tmp'
fa_dir = tmp_dir / 'fasta'
gb_dir = tmp_dir / 'gbk'

def prokkanno(infile_name, input_file, prefix, threads):
    # annotate the fasta sequence by prokka
    cmd = ['prokka', input_file, '--force', '--cpus', str(threads), '--cdsrnaolap', '--prefix', infile_name, 
        '--locustag', prefix, '-o', gb_dir]
    subprocess.run(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

def getgff1(infile_name):
    
    gffile = gb_dir / f'{infile_name}.gff'   # open the gff file

    genome_info = {}   # contain three dict, for every contig, trna pos, and position of every CDS.

    with open(gffile, 'rt') as file:
        for line in file:
            if 'ID=' not in line:
                continue  # skip lines without ID immediately

            info = line.strip().split('\t')
            contig = info[0]

            if contig not in genome_info:
                genome_info[contig] = {
                    'locusdict': {},
                    'trnadict': {},
                    'posdict': {},
                    'plasmid': 'F',
                }

            attributes = dict(
                item.split('=') for item in info[8].split(';') if '=' in item
            )

            ids = attributes.get('ID')
            if not ids:
                continue  # skip if no ID

            genome_info[contig]['locusdict'][ids] = ids
            product = attributes.get('product', '')  # default empty if no product
            pos = [info[3], info[4], info[6], product]

            if info[2] in ('tRNA', 'tmRNA'):
                genome_info[contig]['trnadict'][ids] = pos

            genome_info[contig]['posdict'][ids] = pos

    return genome_info, gb_dir / f'{infile_name}.gbk'

def format_sequence(sequence):
    # Split sequence into lines of 60 characters and join with newline
    return '\n'.join(sequence[i:i+60] for i in range(0, len(sequence), 60))

def getgff(infile_name, input_file, prefix):
    gbfile = pathlib.Path(input_file)
    faafile = gb_dir / f'{infile_name}.faa'
    fnafile = gb_dir / f'{infile_name}.fna'

    records = SeqIO.parse(gbfile, 'genbank')
    genome_info = {}
    i = 1

    with open(faafile, 'wt') as faa, open(fnafile, 'wt') as fna:
        for record in records:
            genome_info[record.id] = {
                'locusdict': {},
                'trnadict': {},
                'posdict': {},
                'plasmid': 'F',
            }

            for feature in record.features:
                if feature.type not in ['CDS', 'rRNA', 'tRNA', 'tmRNA']:
                    continue

                if 'locus_tag' not in feature.qualifiers:
                    continue

                locus_tag = feature.qualifiers['locus_tag'][0]

                # Determine start, end, and strand
                ori = ''
                # if compoundlocation exists, find the start and end, or use freature.location
                if isinstance(feature.location, CompoundLocation):
                    loc_list = sorted(feature.location.parts, key=lambda part: part.start)
                    loc = loc_list[-1]
                    if int(loc.end) == len(record.seq):
                        start = loc_list[-1].start
                        end = loc_list[0].end
                        ori = 'T'
                    else:
                        start = str(int(loc.start))
                        end = str(int(loc.end))
                else:
                    loc = feature.location
                    start = str(int(loc.start))
                    end = str(int(loc.end))
                #.strand is a number, -1 or 1, switch to + or -
                strand = '+' if loc.strand == 1 else '-'

                # Determine product
                if feature.type in ['CDS', 'rRNA']:
                    product = feature.qualifiers.get('product', ['-'])[0]
                elif feature.type == 'tRNA':
                    product = feature.qualifiers.get('product', ['-'])[0]
                else:  # tmRNA
                    product = 'tmRNA'

                newid = gene_key(prefix, i)
                genome_info[record.id]['locusdict'][newid] = locus_tag
                genome_info[record.id]['posdict'][newid] = [start, end, strand, product]

                # Write sequences for CDS/rRNA
                if feature.type in ['CDS', 'rRNA'] and 'translation' in feature.qualifiers:
                    aa_sequence = str(feature.qualifiers['translation'][0])
                    faa.write(f'>{newid} {product}\n')
                    faa.write(f'{format_sequence(aa_sequence)}\n')

                    # write new id, information, and na sequence into fna
                    if not ori:
                        cds_sequence = str(record.seq[feature.location.start:feature.location.end])
                    else:
                        cds_sequence = str(record.seq[start:len(record.seq)]) + str(record.seq[0:end])
                    fna.write(f'>{newid} {product}\n')
                    fna.write(f'{format_sequence(cds_sequence)}\n')

                # Record tRNA
                if feature.type in ['tRNA', 'tmRNA']:
                    genome_info[record.id]['trnadict'][newid] = [start, end, strand, product]

                i += 1

    return genome_info, gbfile

def infomation_filter_marker(genome_info, input_file):
    # Remove contigs with empty data
    genome_info = {k: v for k, v in genome_info.items()
                   if v != {'locusdict': {}, 'trnadict': {}, 'posdict': {}, 'plasmid': 'F'}}

    # Run abricate plasmidfinder
    command = ['abricate', '-db', 'plasmidfinder', str(input_file)]
    process = subprocess.run(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
    
    for line in process.stdout.strip().split('\n'):
        if line.startswith('#') or not line.strip():
            continue
        fields = line.strip().split('\t')
        contig_id = fields[1]
        if contig_id in genome_info:
            genome_info[contig_id]['plasmid'] = 'T'

    return genome_info

def ICE_filter(ICE_res):
    # Dictionaries to track IME genes and positions
    IME_gen_dict = {}    # {sys_id: [gene_name1, gene_name2, ...]}
    IME_pos_dict = {}    # {sys_id: [gene_hit1, gene_hit2, ...]}
    ICE_pos_dict = {}    # {sys_id: [gene_hit1, ...]}

    # Lists to track ICE and AICE systems
    ICE_list = []
    AICE_list = []

    # Open the ICE result file
    with open(ICE_res, 'rt') as file:
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

def get_feature(feature):
    # sort the gene_names as the rule
    featuredict = {
		'Phage_integrase':'Integrase','UPF0236':'Integrase',
        'Recombinase':'Integrase','rve':'Integrase',
        'TIGR02224':'Integrase','TIGR02249':'Integrase',
        'TIGR02225':'Integrase','PB001819':'Integrase',
		'RepSAv2':'Rep','DUF3631':'Rep','Prim-Pol':'Rep',
		'FtsK_SpoIIIE':'Tra'}
    
    if feature in featuredict:
        return f'{featuredict[feature]}@{feature}'
    
    patterns = [
        ('T4SS_MOB', 'Relaxase'), 
        ('Relaxase_', 'Relaxase'), 
        ('t4cp', 'T4CP'), 
        ('tcpA', 'T4CP'), 
        ('FATA_', 'T4SS'), 
        ('FA_', 'T4SS'), 
        ('T4SS_', 'T4SS'), 
    ]
    
    for key, category in patterns:
        if key in feature:
            parts = feature.split('_', 1)
            tag = parts[1] if len(parts) > 1 else feature
            return f'{category}@{tag}'
    
    return f'T4SS@{feature}'   # this return may never be used

def get_MGE(infile_name, infile, genome_info, rootdir): 
    
    ICE_dir = tmp_dir / infile_name / f'{infile_name}_ICE'
    ICE_res = ICE_dir / 'all_systems.tsv'

    # Remove the ICE output folder if it already exists
    if ICE_dir.exists():
        shutil.rmtree(ICE_dir)
    
    ICE_dir.parent.mkdir(parents=True, exist_ok=True)
    
    anno_fa = gb_dir / f'{infile_name}.faa'

    # Run MacSyFinder to detect ICE systems
    model_dir = rootdir / 'data' / 'macsydata'
    ICE_cmd = ['macsyfinder', '--db-type', 'ordered_replicon', '--models-dir', str(model_dir), 
               '--models', 'ICEscan', 'all', '--replicon-topology', 'linear', '--coverage-profile', 
               '0.3', '--sequence-db', str(anno_fa), '-o', str(ICE_dir)]
    
    subprocess.run(ICE_cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    # ICE_res columns: replicon, hit_id, gene_name, hit_pos, model_fqn, sys_id, sys_loci, locus_num
    ftag = ICE_filter(ICE_res)  # List containing sys_id of all predicted ICE, AICE, and IME

    MGEdict = {}
    infodict = {}

    if ftag:
        with open(ICE_res, 'rt') as file:
            for line in file:
                if 'Chromosome' not in line:
                    continue

                lines = line.strip().split('\t')
                
                # Filter out genes not belonging to a cluster
                if lines[7] != '1':  # locus num != 1 means loner genes like phage integrase
                    continue
                elif lines[5] not in ftag:  # sys_id must be included in ftag
                    continue
                else:
                    tags = get_feature(lines[2])  # lines[2] is gene_name; tags = {category}@{tag}
                    
                    # Extract T4SS type if present
                    mpf = lines[4].split('/')[-1].split('_')[1] if 'T4SS' in lines[4] else ''
                    # Extract relaxase information
                    mob = tags.split('@')[1] if 'Relaxase@' in tags else ''    
                    
                    # Determine category based on sys_id
                    if 'IME' in lines[5]:
                        category = 'IME'
                    elif 'AICE' in lines[5]:
                        category = 'AICE'
                    else:
                        category = 'ICE'
                    
                    orig_num = lines[5].split('_')[-1]

                    # Store gene tags
                    MGEdict.setdefault(f'{category}{orig_num}', {})[lines[1]] = tags
                    # Store accessory info: relaxase and T4SS
                    infodict.setdefault(f'{category}{orig_num}', {'mob': set(), 'mpf': set()})
                    if mob:
                        infodict[f'{category}{orig_num}']['mob'].add(mob)
                    if mpf:
                        infodict[f'{category}{orig_num}']['mpf'].add(mpf)

    return MGEdict, infodict

def remove_loner(MGE_dict, threshold = 20):
    # Exclude loner genes from MGE candidates
    new_MGE_dict = {}
    for sys_id, sys_data in MGE_dict.items():
        hit_ids = list(sys_data['hit_id'].keys())
        if len(hit_ids) < 2:
            continue
            # Keep only hits that are close to at least one other hit within the threshold
        filtered_hit_id = sorted([
            x for x in hit_ids
            if any(
                abs(int(x.split('_')[1]) - int(y.split('_')[1])) <= threshold
                and int(x.split('_')[1]) != int(y.split('_')[1])
                for y in hit_ids
            )
        ])
        if filtered_hit_id:
            # Retain only the filtered hits
            new_MGE_dict[sys_id] = sys_data.copy()
            new_MGE_dict[sys_id]['hit_id'] = {
                k: v for k, v in sys_data['hit_id'].items() if k in filtered_hit_id
            }

    # Remove duplicate MGEs based on identical hit positions
    ICE_type = ''
    seen = set()
    unique_dict = {}

    for sys_id in new_MGE_dict.keys():
        # Detect when a new ICE type starts
        if ICE_type == '' or sys_id.rsplit('_', 1)[0] != ICE_type:
            ICE_type = sys_id.rsplit('_', 1)[0]
            seen = set()
        t = tuple(new_MGE_dict[sys_id]['hit_id'].keys())
        if t not in seen:
            unique_dict[sys_id] = new_MGE_dict[sys_id]
            seen.add(t)

    return unique_dict

def validate_hits(MGE_dict, genome_info):
    # Verify if all hits of an ICE are located on the same contig and not at the edges.
    # Remove false matches. Only ICEs at contig edges (potential breaks) or covering the whole contig are retained.
    
    validated_MGE_dict = {}
    for MGE, info in MGE_dict.items():  # MGE_dict[sys_id] = {'gene_name': [], 'hit_id': []}
        record_ids = []
        hit_status = {}
        for hit_id in info['hit_id']:
            for record_id, ginfo in genome_info.items():
                tmpids = list(ginfo['locusdict'].keys())
                for i, tmpid in enumerate(tmpids):
                    if hit_id == tmpid:
                        record_ids.append(record_id)
                        # Mark if the gene is not at the contig edge
                        if 20 <= i <= (len(tmpids) - 20):
                            hit_status[hit_id] = 'V'  # valid
                        else:
                            hit_status[hit_id] = 'P'  # potential edge
                        break

        # Determine if the ICE might be fake
        break_status = set(hit_status.values())
        if len(set(record_ids)) == 1:  # all hits are on a single contig
            contig_id = next(iter(record_ids))
            if genome_info[contig_id]['plasmid'] != 'T' and 'P' in break_status:
                new_entry = info.copy()
                new_entry['hit_id'] = hit_status
                validated_MGE_dict[MGE] = new_entry
                continue  # keep ICE at edge of chromosome
                # discard ICE on plasmid or fully internal
                # hits are scattered across multiple contigs → discard

    return validated_MGE_dict
    
def check_ICE_genes(ICE_dict, genome_info):
    # This function attempts to recover potential ICEs based on mandatory and accessory genes.
    # Mandatory genes must be present in the combined end regions of ICEs.
    # Accessory genes increase confidence and are used to incorporate middle ICE regions.

    mandatory_genes = {
      'T4SS_virb4' : ['T4SS_virb4', 'T4SS_I_traU'], 
      'T4SS_t4cp1' : ['T4SS_tcpA', 'T4SS_t4cp1', 'T4SS_t4cp2'], 
      # despite 1 extra gene tcpA exists in T4SS_typeFA, this gene will not exist in the result of 
      # other type of ICE, so it won't affect the result of others
      'T4SS_MOBV' : ['T4SS_MOBV', 'T4SS_MOBB', 'T4SS_MOBT', 'T4SS_MOBQ', 'T4SS_MOBP3', 
                     'T4SS_MOBP2', 'T4SS_MOBP1', 'T4SS_MOBH', 'T4SS_MOBF', 'T4SS_MOBC', 
                     'Relaxase_firmi_MOBL', 'Relaxase_firmi_Rep_2', 'Relaxase_firmi_Viral_Rep_A', 
                     'Relaxase_firmi_Viral_Rep_B1', 'Relaxase_firmi_Viral_Rep_B2', 
                     'Relaxase_PHA_IME_A1', 'Relaxase_PHA_IME_B', 'Relaxase_profile_MOBT'], 
      'Phage_integrase' : ['Phage_integrase', 'UPF0236', 'Recombinase', 'rve', 'TIGR02224', 
                           'TIGR02249', 'TIGR02225', 'PB001819'], 
      }
    
    accessory_genes = {
        'T4SS_typeB' : ['T4SS_B_traE', 'T4SS_B_traF', 'T4SS_B_traH', 'T4SS_B_traI', 'T4SS_B_traJ', 
                   'T4SS_B_traK', 'T4SS_B_traL', 'T4SS_B_traM', 'T4SS_B_traN', 'T4SS_B_traO', 
                   'T4SS_B_traP', 'T4SS_B_traQ'], 
        'T4SS_typeC' : ['T4SS_C_alr7204', 'T4SS_C_alr7205', 'T4SS_C_alr7207', 'T4SS_C_alr7208', 
                       'T4SS_C_alr7209', 'T4SS_C_alr7210', 'T4SS_C_alr7211', 'T4SS_C_alr7212'], 
        'T4SS_typeF' : ['T4SS_F_traB', 'T4SS_F_traE', 'T4SS_F_traF', 'T4SS_F_traH', 'T4SS_F_traK', 
                       'T4SS_F_traL', 'T4SS_F_traN', 'T4SS_F_traU', 'T4SS_F_traV', 'T4SS_F_traW'], 
        'T4SS_typeFA' : ['FA_orf13', 'FA_orf14', 'FA_orf15', 'FA_orf17a', 'FA_orf17b', 'FA_orf19', 
                       'FA_orf23'], 
        'T4SS_typeFATA' : ['FATA_cd411', 'FATA_cd419_1', 'FATA_cd419a', 'FATA_cd419b', 'FATA_cd424', 
                       'FATA_gbs1346', 'FATA_gbs1347', 'FATA_gbs1350', 'FATA_gbs1354', 'FATA_gbs1365', 
                       'FATA_gbs1369', 'FATA_prgB', 'FATA_prgC', 'FATA_prgF', 'FATA_prgHa', 
                       'FATA_prgHb', 'FATA_prgIa', 'FATA_prgIb', 'FATA_prgIc', 'FATA_prgK', 
                       'FATA_prgL', 'FATA_trsC', 'FATA_trsD', 'FATA_trsF', 'FATA_trsG', 
                       'FATA_trsJ', 'FATA_trsL'], 
        'T4SS_typeG' : ['T4SS_G_tfc7', 'T4SS_G_tfc8', 'T4SS_G_tfc9', 'T4SS_G_tfc10', 'T4SS_G_tfc11', 
                       'T4SS_G_tfc12', 'T4SS_G_tfc13', 'T4SS_G_tfc14', 'T4SS_G_tfc15', 'T4SS_G_tfc17', 
                       'T4SS_G_tfc18', 'T4SS_G_tfc19', 'T4SS_G_tfc22', 'T4SS_G_tfc23', 'T4SS_G_tfc24', 
                       'T4SS_G_tfc2', 'T4SS_G_tfc3', 'T4SS_G_tfc5', 'T4SS_I_traE'], 
        'T4SS_typeI' : ['T4SS_I_traI', 'T4SS_I_traK', 'T4SS_I_traL', 'T4SS_I_traM', 'T4SS_I_traN', 
                       'T4SS_I_traO', 'T4SS_I_traP', 'T4SS_I_traQ', 'T4SS_I_traR', 'T4SS_I_traW', 
                       'T4SS_I_traY', 'T4SS_I_trbA', 'T4SS_I_trbB'], 
        'T4SS_typeT' : ['T4SS_T_virB1', 'T4SS_T_virB10', 'T4SS_T_virB11', 'T4SS_T_virB2', 
                       'T4SS_T_virB3', 'T4SS_T_virB5', 'T4SS_T_virB6', 'T4SS_T_virB8', 
                       'T4SS_T_virB9'], 
        }
    
    ICE_types = sorted({'_'.join(ice.split('_')[-3:-1]) for ice in ICE_dict.keys()})    
    
    recovery_ICE = defaultdict(lambda: {'end': [], 'middle': [], 'middle_importance' : ''})
    ICE_count = 1
    for ICE_type in ICE_types:
        end = []
        middle = []
        for ICE, info in ICE_dict.items():
            if ICE.rsplit('_', 2)[1] != ICE_type:
                continue
            break_markers = list(info['hit_id'].values())
            if all(x == 'P' for x in break_markers):
                hit_ids = list(info['hit_id'].keys())
                hit_numbers = [int(hit.split('_')[1]) for hit in hit_ids]
                for contig_info  in genome_info.values():
                    tmpids = list(contig_info ['locusdict'].keys())
                    if hit_ids in tmpids:
                        record_numbers = [int(tmpid.split('_')[1]) for tmpid in tmpids]
                        if min(hit_numbers) - min(record_numbers) <= 20 and max(record_numbers) - max(hit_numbers) <= 20:
                            middle.append(ICE)
                        else:
                            end.append(ICE)
                        break
            else:
                end.append(ICE)
        # a ICE could only have 2 end, but how many middle exist is a problem
        # also, may be some end will be predicted as middle, this situation haven't be included for now
        if len(end) >= 2:
        
            for pair in list(combinations(end, 2)):
                all_genes = set(chain.from_iterable(ICE_dict[p]['gene_name'] for p in pair))
                core_present = all(any(gene in all_genes for gene in genes) for genes in mandatory_genes.values())
                # if there are one group in mandatory haven't be found
                if not core_present:
                    continue
                count = sum(1 for gene in accessory_genes[ICE_type] if gene in all_genes)
                if middle:   # if there are a potential middle of ICE
                    all_middles = list(chain.from_iterable(combinations(middle, r) for r in range(1, len(middle) + 1)))
                    best_middle = tuple()
                    importance_mark = ''
                    for one_middle in all_middles:
                        middle_gene = set(chain.from_iterable(ICE_dict[m]['gene_name'] for m in one_middle))
                        new_count = count + sum(1 for gene in accessory_genes[ICE_type] if gene in middle_gene)
                        if count >= 4 and new_count >= 4:
                            best_middle = one_middle
                            importance_mark = 'L'
                        elif count < 4 and new_count >= 4:
                            best_middle = one_middle
                            importance_mark = 'H'
                        else:
                            continue
                    if best_middle:
                        recovery_ICE[f'rICE{ICE_count}']['end'] = list(pair)
                        recovery_ICE[f'rICE{ICE_count}']['middle'] = list(best_middle)
                        recovery_ICE[f'rICE{ICE_count}']['middle_importance'] = importance_mark
                        ICE_count += 1
                else:
                    if count >= 4:
                        recovery_ICE[f'rICE{ICE_count}']['end'] = list(pair)
                        ICE_count += 1
                    else:
                        continue

    return recovery_ICE

def ICE_recovery(infile_name, genome_info):
    # For draft genomes, MacSyFinder may miss some ICEs. This function recovers ICEs 
    # from rejected results.
    ICE_dir = tmp_dir / infile_name / f'{infile_name}_ICE'
    ICE_rej = ICE_dir / 'rejected_candidates.tsv'

    # Only 2 genes are required for an IME, so recovery is usually not necessary, but included here.
    IME_dict = defaultdict(lambda: {'gene_name': [], 'hit_id': {}})
    ICE_dict = defaultdict(lambda: {'gene_name': [], 'hit_id': {}})

    # Parse rejected candidates
    with open(ICE_rej, 'rt') as file:
        for line in file:
            if 'Chromosome' not in line:  # XML results are located in ./ICEscan/Chromosome
                continue

            lines = line.strip().split('\t')
            sys_id, hit_id, gene_name = lines[0], lines[4], lines[6]  
            # sys_id e.g., aHPS7_T4SS_typeG_2, aHPS7_IME_1; lines[6] is gene_name

            if 'IME' in sys_id:
                IME_dict[sys_id]['gene_name'].append(gene_name)
                IME_dict[sys_id]['hit_id'][hit_id] = ''
            elif 'T4SS' in sys_id:
                ICE_dict[sys_id]['gene_name'].append(gene_name)
                ICE_dict[sys_id]['hit_id'][hit_id] = ''
    
    # Remove loner IMEs and ICEs
    IME_dict = remove_loner(IME_dict) if IME_dict else {}
    ICE_dict = remove_loner(ICE_dict) if ICE_dict else {}

    # Validate hits against genome information
    IME_dict = validate_hits(IME_dict, genome_info) if IME_dict else {}
    ICE_dict = validate_hits(ICE_dict, genome_info) if ICE_dict else {}

    # Check for ICEs with complete gene sets
    recovery_ICE = check_ICE_genes(ICE_dict, genome_info) if ICE_dict else {}

    # List of end positions to consider when screening for DR (ignore middle positions)
    partial = [x for v in recovery_ICE.values() for x in v['end']]

    rICEdict = {}
    rinfodict = {}

    if recovery_ICE:
        # Re-scan rejected candidates to fill recovery dictionaries
        with open(ICE_rej, 'rt') as file:
            for line in file:
                if 'Chromosome' not in line:
                    continue

                lines = line.strip().split('\t')
                if lines[0] not in partial:
                    continue

                tags = get_feature(lines[6])  # lines[6] is gene_name, tags = {category}@{tag}

                # Extract T4SS type if present
                mpf = ''
                if 'T4SS' in lines[2]:  # lines[2] is model_fqn, e.g., ICEscan/Chromosome/T4SS_typeG
                    mpf = lines[2].split('/')[-1].split('_')[1]

                # Extract relaxase if present
                mob = ''
                if 'Relaxase@' in tags:
                    mob = tags.split('@')[1]

                for key, value in recovery_ICE.items():
                    if lines[0] in value['end'] or lines[0] in value['middle']:
                        if lines[4] not in list(ICE_dict[lines[0]]['hit_id'].keys()):
                            continue
                        rICEdict.setdefault(key, {}).setdefault(lines[0], {})[lines[4]] = tags
                        rinfodict.setdefault(key, {'mob': [], 'mpf': []})

                        if mob and mob not in rinfodict[key]['mob']:
                            rinfodict[key]['mob'].append(mob)
                        if mpf and mpf not in rinfodict[key]['mpf']:
                            rinfodict[key]['mpf'].append(mpf)

    return recovery_ICE, ICE_dict, rICEdict, rinfodict

def get_DR(infile_name, input_file, reverse = False):
    
    DRindex = tmp_dir / infile_name / f'{infile_name}_DR'
    # fasta and gbk can both be the input of mkvtree
    # Build k-mer index for the input sequence
    mkvtree_cmd = [
        'mkvtree', '-db', str(input_file), '-indexname', str(DRindex), 
        '-dna', '-pl', '-allout'
    ]
    subprocess.run(mkvtree_cmd, check=True)
    vmatch_cmd = ['vmatch', '-l', '15', str(DRindex)]
    process = subprocess.run(vmatch_cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out = process.stdout.decode()
    
    # Optionally search the reverse complement
    if reverse:
        vmatch_rev_cmd = ['vmatch', '-p', '-l', '15', str(DRindex)]
        process = subprocess.run(vmatch_rev_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out += process.stdout.decode()
    # DRout : 22    0 1228868   D    22    0 1736928   0    8.87e-02     44   100.00
    # length_left, seq_num_i, post_left, strand, length_right, seq_num_j, pos_right, distance, E-value, score, identity
    DRlist = []
    for line in out.strip().split('\n'):
        if not line.startswith('#'):
            lines = line.strip().split()
            if not lines:
                continue
            
            left_start = int(lines[2]) + 1
            left_end   = int(lines[2]) + int(lines[0])
            right_start = int(lines[6]) + 1
            right_end   = int(lines[6]) + int(lines[4])
            
            # Filter DRs that are too long or too short
            if (right_end - left_start) > 500000 or (right_end - left_start) < 5000:
                continue
            DRlist.append(f'{left_start},{left_end},{right_start},{right_end}')

    return DRlist 

def remove_duplicate(MGEdict, rICEdict):
    # maybe some IME predicted by get_MGE is part of the recovered ICE.
    remove_keys = []
    for MGE in MGEdict:
        hit_list = list(MGEdict[MGE].keys())
        for rICE in rICEdict.values():
            rhit_list = []
            for i in rICE.values():
                rhit_list += list(i.keys())
            if all(gene in rhit_list for gene in hit_list):
                remove_keys.append(MGE)

    for k in remove_keys:
        MGEdict.pop(k, None)
    return MGEdict

def find_max_distance(numbers):
	max_distance = -1
	max_distance_index = -1

	for i in range(len(numbers) - 1):
		distance = abs(numbers[i] - numbers[i+1])
		if distance > max_distance:
			max_distance = distance
			max_distance_index = i

	if max_distance_index == -1:
		return None

	return numbers[max_distance_index], numbers[max_distance_index + 1]

def pos_tag(pos, posdict, edge_point, final_edge, max_in_contig, dirtag):

    # posdict[newid] = [start, end, strand, product]
	for hit_id, info in posdict.items():
		vstart, vend = int(info[0]), int(info[1])
		if vend >= int(pos):   # find whether gene is in the end of ICE. if end > DR
			if dirtag == 's':
				edge_point = int(hit_id.split('_')[1])
				final_edge = max(max_in_contig, edge_point - 5)
			else:
				if vstart > int(pos):   # if the start of the gene also < DR
					edge_point = int(hit_id.split('_')[1]) - 1   # which means the last gene is the end
				else:
					edge_point = int(hit_id.split('_')[1])   # or, the current gene is the end
				final_edge = min(max_in_contig, edge_point + 5)
			break   # when found the gene, the search is over.
	return edge_point, final_edge

def merge_tRNA(info, DRlist, prefix, genome_info):
    
    start_point = int(next(iter(info)).split('_')[1])
    end_point = int(next(reversed(info)).split('_')[1])   # python 3.8+ needed
    
    # ori problem could be solved by recovery, so we don't need to consider it here.
    # The authors defined that if a tRNA gene is located near the predicted ICE region, it would be 
    # considered as a new boundary. However, they set the threshold at only five genes beyond the 
    # boundary, which seems insufficient. In their original version, a region of 250 kb around the 
    # integrase was considered. Although the algorithm differs, based on my personal experience, 
    # using only five genes is not enough. The current algorithm relies on ConjScan as its core, 
    # which only detects T4SS-related structures, while many other types of genes can also be present 
    # within an ICE.
    # In one of my cases, a tRNA gene was found 16 genes away from the predicted ICE, and its 
    # excision/integration form was experimentally verified by PCR. Therefore, I will temporarily 
    # set the threshold to 20 genes to test the accuracy.
    
    record = None
    # Sometimes, at this step, an ICE may be predicted to be located on a broken contig. However, 
    # this is incorrect, because the input for macsyFinder is an FAA file, which does not record 
    # contig information and instead concatenates proteins from different contigs in a disordered 
    # manner. Therefore, such results need to be filtered out. 
    # If an ICE truly lies on a broken contig, it will be detected later during the recovery step.
    
    for key, value in genome_info.items():
        hit_list = list(value['locusdict'].keys())
        if all(gene in hit_list for gene in list(info.keys())):
            record = key
            max_in_contig = int(hit_list[-1].split('_')[1])
            min_in_contig = int(hit_list[0].split('_')[1])
            break
    
    if record is None:
        return None
    
    expanded_start = max(min_in_contig, start_point - 20)
    expanded_end = min(max_in_contig, end_point + 20)
        
    ICE_region = [expanded_start, expanded_end]
    
    trnalist = []
    for gene_id, trna_info in genome_info[record]['trnadict'].items():   # trnadict[newid] = [start, end, strand, product]
        if expanded_start <= int(gene_id.split('_')[1].lstrip('0')) <= expanded_end:   # if trna within ICE
            ICE_region.append(int(gene_id.split('_')[1].lstrip('0')))
            trnalist.append(trna_info)   # eg. trna_info = ['18651', '18726', '+', 'tRNA-Glu(ttc)']
            
    ICE_region.sort()
    finalstart, finalend = find_max_distance(ICE_region)   # find the two genes with largest gap
    
    # here has a prerequisite, which is trna will not located within the ICE
    # the background is: normally ICE will insert into the genome nearby the tRNA
    DR1 = genome_info[record]['posdict'][gene_key(prefix, start_point)][0]   # start pos of ICE  posdict[newid] = [start, end, strand, product]
    DR2 = ''
    DR3 = ''
    DR4 = genome_info[record]['posdict'][gene_key(prefix, end_point)][1]   # end pos of ICE
    
    if trnalist:
        if finalend == expanded_end:   # right end hasn't change
            start_point = finalstart   # maybe add of tRNA will change the left end, re-define it
            finalstart =  max(min_in_contig, finalstart - 5)
            start_edge_gene = genome_info[record]['posdict'][gene_key(prefix, start_point)]
            start, end = int(start_edge_gene[0]), int(start_edge_gene[1])
            DR1 = start   # start_point may has been changed
            for DR in DRlist:
                DRs = DR.split(',')
                # if the left DR is inside the boundary gene
                if start < int(DRs[0]) < end:
                    end_point, finalend = pos_tag(DRs[3], genome_info[record]['posdict'], end_point, finalend, max_in_contig, 'e')
                    DR1 = DRs[0]
                    DR2 = DRs[1]
                    DR3 = DRs[2]
                    DR4 = DRs[3]
                    break
                            
        elif finalstart == expanded_start:
            end_point = finalend
            finalend = min(max_in_contig, finalend + 5)
            end_edge_gene = genome_info[record]['posdict'][gene_key(prefix, end_point)]
            start, end = int(end_edge_gene[0]), int(end_edge_gene[1])
            DR4 = start			
            for DR in DRlist:
                DRs = DR.split(',')
                if start < int(DRs[3]) < end:
                    start_point, finalstart = pos_tag(DRs[0], genome_info[record]['posdict'], start_point, finalstart, max_in_contig, 's')
                    DR1 = DRs[0]
                    DR2 = DRs[1]
                    DR3 = DRs[2]
                    DR4 = DRs[3]
                    break
                
    return DictMGE([DR1, DR2, DR3, DR4, start_point, end_point, finalstart, finalend, trnalist, record])

class DictMGE(object):
    
    def __init__(self, info):
        self.DR1 = info[0]
        self.DR2 = info[1]
        self.DR3 = info[2]
        self.DR4 = info[3]
        self.start_point = info[4]
        self.end_point = info[5]
        self.finalstart = info[6]
        self.finalend = info[7]
        self.trnalist = info[8]
        self.record = info[9]

def MGE_reorder(MGEdict):

    new_MGEdict = {}
    counters = {}
    
    for key, value in MGEdict.items():
        prefix = ''.join([c for c in key if not c.isdigit()])  # eg. 'IME'
        counters[prefix] = counters.get(prefix, 0) + 1          # count
        new_key = f"{prefix}{counters[prefix]}"                 # IME1、ICE2
        new_MGEdict[new_key] = value
    
    new_MGEdict
    
    return new_MGEdict

def search_trna(contig, start_gene, direction, genome_info, trnalist):

    locus_ids = list(genome_info[contig]['locusdict'])
    min_in_contig = int(locus_ids[0].split('_')[1])
    max_in_contig = int(locus_ids[-1].split('_')[1])

    if direction == 'left':
        expanded = max(min_in_contig, start_gene - 20)
        in_range = lambda gid: expanded <= int(gid.split('_')[1]) <= start_gene
    else:
        expanded = min(max_in_contig, start_gene + 20)
        in_range = lambda gid: start_gene <= int(gid.split('_')[1]) <= expanded

    final_gene = start_gene
    for gene_id, info in genome_info[contig]['trnadict'].items():
        if in_range(gene_id):
            trnainfo = [contig] + info
            trnalist.append(trnainfo)
            final_gene = gene_id.split('_')[1]
            
    return final_gene

def gene_key(prefix, idx): return f'{prefix}_{str(idx).zfill(5)}'

def boundary_of_rICE(infile_name, gbfile, recovery_ICE, ICE_dict, genome_info, rICEdict, prefix):
    
    rdictICE = {}
    for ICE, info in recovery_ICE.items():
        
        end0 = ICE_dict[info['end'][0]]
        end1 = ICE_dict[info['end'][1]]
        
        end0_last = list(end0['hit_id'].values())[-1]
        end1_last = list(end1['hit_id'].values())[-1]

        left_contig = find_located_contig(ICE_dict[info['end'][0]]['hit_id'].keys(), genome_info)
        right_contig = find_located_contig(ICE_dict[info['end'][1]]['hit_id'].keys(), genome_info)

        trnalist = []
        if end0_last == 'P' and end1_last == 'P':
            left_gene = int(next(iter(rICEdict[ICE][info['end'][0]])).split('_')[1])
            right_gene = int(next(iter(rICEdict[ICE][info['end'][1]])).split('_')[1])
        
            final_left = search_trna(left_contig, left_gene, 'left', genome_info, trnalist)
            final_right = search_trna(right_contig, right_gene, 'left', genome_info, trnalist)
        
        elif end0_last != 'P' and end1_last != 'P':
            left_gene = int(next(reversed(rICEdict[ICE][info['end'][0]])).split('_')[1])
            right_gene = int(next(reversed(rICEdict[ICE][info['end'][1]])).split('_')[1])
        
            final_left = search_trna(left_contig, left_gene, 'right', genome_info, trnalist)
            final_right = search_trna(right_contig, right_gene, 'right', genome_info, trnalist)
        
        else:
            if end0_last == 'P':
                left_gene = int(next(iter(rICEdict[ICE][info['end'][0]])).split('_')[1])
                right_gene = int(next(reversed(rICEdict[ICE][info['end'][1]])).split('_')[1])
        
                final_left = search_trna(left_contig, left_gene, 'left', genome_info, trnalist)
                final_right = search_trna(right_contig, right_gene, 'right', genome_info, trnalist)
            else:
                left_gene = int(next(iter(rICEdict[ICE][info['end'][1]])).split('_')[1])
                right_gene = int(next(reversed(rICEdict[ICE][info['end'][0]])).split('_')[1])
        
                final_left = search_trna(right_contig, left_gene, 'left', genome_info, trnalist)
                final_right = search_trna(left_contig, right_gene, 'right', genome_info, trnalist)

        DR1, DR2, DR3, DR4 = '', '', '', ''
        if end0_last == 'P' and end1_last == 'P':
            DR1 = genome_info[left_contig]['posdict'][gene_key(prefix, final_left)][0]   # start pos of ICE  posdict[newid] = [start, end, strand, product]
            DR4 = genome_info[right_contig]['posdict'][gene_key(prefix, final_right)][0]
        elif end0_last != 'P' and end1_last != 'P':
            DR1 = genome_info[left_contig]['posdict'][gene_key(prefix, final_left)][1]
            DR4 = genome_info[right_contig]['posdict'][gene_key(prefix, final_right)][1]
        else:
            if end0_last == 'P':
                DR1 = genome_info[left_contig]['posdict'][gene_key(prefix, final_left)][0]
                DR4 = genome_info[right_contig]['posdict'][gene_key(prefix, final_right)][1]
            else:
                DR1 = genome_info[right_contig]['posdict'][gene_key(prefix, final_left)][0]
                DR4 = genome_info[left_contig]['posdict'][gene_key(prefix, final_right)][1]
        
        records = SeqIO.parse(gbfile, 'genbank')
        genome_sequence = {}
        for record in records:
            genome_sequence[record.id] = record.seq
        
        if end0_last == 'P' and end1_last == 'P':
            start_seq = str(genome_sequence[left_contig])
            end_seq = str(genome_sequence[right_contig].reverse_complement())
        elif end0_last != 'P' and end1_last != 'P':
            start_seq = str(genome_sequence[left_contig].reverse_complement())
            end_seq = str(genome_sequence[right_contig])
        else:
            if end0_last == 'P':
                start_seq = str(genome_sequence[left_contig])
                end_seq = str(genome_sequence[right_contig])
            else:
                start_seq = str(genome_sequence[right_contig])
                end_seq = str(genome_sequence[left_contig])
                
        len_left, len_right = len(start_seq), len(end_seq)
        
        if trnalist:
            temp_for_DR = tmp_dir/ infile_name / f'refinder_{ICE}'
            
            with open(temp_for_DR, 'wt') as file:   # write the start and end of the ICE to a file, rescreen the DR
                file.write('>temp_DR\n')
                file.write(start_seq)
                file.write(end_seq)
            
            rDRlist = get_DR(infile_name, temp_for_DR, reverse = True)   # in case the direct of contig is wrong

            left_pos = genome_info[left_contig]['posdict'][gene_key(prefix, final_left)]
            right_pos = genome_info[right_contig]['posdict'][gene_key(prefix, final_right)]
            left_range = lambda x: int(left_pos[0]) < x < int(left_pos[1])
            right_range = lambda x: int(right_pos[0]) < x < int(right_pos[1])
            
            if final_right == right_gene or final_left == left_gene:
                for DR in rDRlist:
                    DRs = DR.split(',')
                    d0, d1, d2, d3 = map(int, DRs)
            
                    if final_right == right_gene:
                        if end0_last == 'P' and end1_last == 'P':
                            cond = left_range(d0)
                            shift = len_left
                        elif end0_last != 'P' and end1_last != 'P':
                            cond = left_range(len_left - d1)
                            shift = len_left
                        elif end0_last == 'P':  # end1_last != 'P'
                            cond = left_range(d0)
                            shift = len_left
                        else:  # end0_last != 'P' and end1_last == 'P'
                            cond = right_range(d0)
                            shift = len_right
            
                    else:  # final_left == left_gene
                        if end0_last == 'P' and end1_last == 'P':
                            cond = right_range(len_right - d2 + len_left)
                            shift = len_left
                        elif end0_last != 'P' and end1_last != 'P':
                            cond = right_range(d3 - len_left)
                            shift = len_left
                        elif end0_last == 'P':  # end1_last != 'P'
                            cond = right_range(d3 - len_left)
                            shift = len_left
                        else:  # end0_last != 'P' and end1_last == 'P'
                            cond = left_range(len_left - d2 + len_right)
                            shift = len_right
            
                    if cond:
                        DR1, DR2 = DRs[0], DRs[1]
                        DR3, DR4 = str(d2 - shift), str(d3 - shift)
                        break
                        
        rdictICE[ICE] = RdictICE([DR1, DR2, DR3, DR4, final_left, final_right, end0_last, end1_last, 
                                 trnalist, left_contig, right_contig, len_left, len_right])
    
    return rdictICE    

class RdictICE(object):
    
    def __init__(self, info):
        self.DR1 = info[0]
        self.DR2 = info[1]
        self.DR3 = info[2]
        self.DR4 = info[3]
        self.final_left = info[4]
        self.final_right = info[5]
        self.end0_last = info[6]
        self.end1_last = info[7]
        self.trnalist = info[8]
        self.left_contig = info[9]
        self.right_contig = info[10]
        self.len_left = info[11]
        self.len_right = info[12]

def find_located_contig(hit_list, genome_info):
    
    for contig_id, info in genome_info.items():
        tmpids = list(info['locusdict'].keys())
        if all(gene in tmpids for gene in list(hit_list)):
            contig = contig_id
    
    return contig

def get_args(argdict, vfdict, isdict, dfdict, metaldict, popdict, symdict, gene, feature, product):
    # search if genes belong to these dict, if so, add correspond feature and product info
    feature = [feature]
    product = [product]
    
    dict_label_map  = {
        'AR' : argdict, 
        'VF' : vfdict, 
        'IS' : isdict, 
        'Defense' : dfdict, 
        'Metal' : metaldict, 
        'Degradation' : popdict, 
        'Symbiosis' : symdict, 
        }

    for label, dic in dict_label_map.items():
        if gene in dic:
            feature.append(label)
            product.append(dic[gene])
    
    feature = '; '.join(list(filter(None, feature)))
    product = '; '.join(list(filter(None, product)))
    
    return feature, product

def process_gene(gene, posdict, ICEdict, argdict, vfdict, isdict, dfdict, metaldict, popdict, 
                 symdict, locusdict, feature_default = 'Flank', is_ICE = False, ICEdict_tag = None):
    
    s, e, strand, pro = posdict[gene]   # posdict[newid] = [start, end, strand, product]
    pos = f'{s}..{e} [{strand}], {int(e) - int(s) + 1}'   # eg. 1945..2222 [+], 278

    feature = feature_default
    product = pro

    if is_ICE and gene in ICEdict[ICEdict_tag]:   # if the gene within ICE
        if gene in ICEdict[ICEdict_tag]:
            feature, pro11 = ICEdict[ICEdict_tag][gene].split('@')   # tags = {category}@{tag}
        else:
            feature,pro11 = '', ''   # if the gene within ICE but not belong to T4SS
        if pro11:   # if it is belong to T4SS
            if pro == 'hypothetical protein':
                product = pro11
            else:
                product = f'{pro}, {pro11}'

    feature, product = get_args(argdict, vfdict, isdict, dfdict, metaldict, popdict, symdict, gene, feature, product)

    if 'hypothetical protein;' in product:
        product = product.replace('hypothetical protein;', '')

    return {
        'gene': locusdict[gene],
        'pos': pos,
        'prod': product,
        'featu': feature
    }

def gc(input_file, contig, start, end, filetype):
    # the original writer haven't consider the filetype here, however, now draft genome still not be considered
    if filetype == 'fa':
        records = SeqIO.parse(input_file, 'fasta')
    else:
        records = SeqIO.parse(input_file, 'gb')
    for record in records:
        if record.id == contig:
            if start == 0:
                start = 1
            sequence = record.seq[int(start) - 1 : int(end)]
        
    gcs = f'{gc_fraction(sequence) * 100:.2f}'
    
    return gcs

def getfa(fasta_file, start, end, contig, filetype):
    # the original writer haven't consider the filetype here, however, now draft genome still not be considered
    if filetype == 'fa':
        records = SeqIO.parse(fasta_file, 'fasta')
    else:
        records = SeqIO.parse(fasta_file, 'gb')
    
    for record in records:
        if record.id == contig:
        	sequence = record.seq[int(start) - 1 : int(end)]
    
    return str(sequence)

def oritseq(infile_name, regi, input_file, contig, start, end, filetype, rootdir, threads):
    
    oritseq = '-'
    fafile = tmp_dir / infile_name / f'{regi}_fororit.fa'   # whole length sequence of ICE
    with open(fafile,'wt') as orif:
        seq = getfa(input_file, start, end, contig, filetype)
        orif.write('>fororit\n')
        orif.write(seq)
        
    oriT_Database = rootdir / 'data' / 'oriT_db'
    # search oriT
    blast_cmd = ['blastn', '-db', str(oriT_Database), '-query', str(fafile), '-evalue', '0.01', 
                 '-num_threads', str(threads), '-word_size', '11', '-outfmt', '6 std slen stitle', 
                 '-num_alignments', '1']

    process = subprocess.run(blast_cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out = process.stdout.decode()
    
    for line in out.strip().split('\n'):
        lines = line.strip().split('\t')
        if lines[0]:   # if there is a matched oriT
            coverage = int(lines[3]) / int(lines[12])
            identity = float(lines[2]) / 100
            havalue = coverage * identity
            if havalue > 0.49:
                oritseq = getfa(fafile, str(int(lines[6])-1), lines[7], 'fororit', 'fa')
                # lines[6] matched start, lines[7] matched end
                break
    
    return oritseq

def getcolor(feature, product):
    
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
            feature = 'Hyp' if product == 'hypothetical protein' else 'Gene'
        else:
            feature = 'Gene'
    
    return coldict[feature], namedict[feature]

def calculate_gc(fasta_file, contig, filetype, start, end, window_size = 500, step_size = 50):
    
    if filetype == 'fa':
        records = SeqIO.parse(fasta_file, 'fasta')
    else:
        records = SeqIO.parse(fasta_file, 'gb')
    
    for record in records:
        if record.id == contig:
        	if start == 0:
        		start = 1
        	sequence = record.seq[start-1:end]

	# windows = []
    gc_contents = []
    pos = []
    j = start/1000 + 0.025
    for i in range(0, len(sequence) - window_size + 1, step_size):
        window = sequence[i:i+window_size]
        gc_content = (gc_fraction(window) * 100)
        gc_contents.append(gc_content)
        pos.append(round(j, 4))
        j += 0.05
    gcdict = {
		        'xData':pos,
		        'datasets':[{
		            'name':'',
		            'data':gc_contents,
		            'unit':'%',
		            'type':'line',
		            "valueDecimals": 1
		        }]
	    }
    return gcdict

def get_map(infile_name, input_file, MGEdict, infodict, genome_info, dictMGE, filetype, prefix, rootdir, out_dir, threads):
    
    final_dir = workdir / out_dir / infile_name
    if not final_dir.exists():
        final_dir.mkdir(parents = True)
    js_dir = workdir / out_dir / infile_name / 'js'
    if not js_dir.exists():
        js_dir.mkdir(parents = True)
    viewfile = rootdir / 'script' / 'js' / 'view.html'
    gcmap = rootdir / 'script' / 'js' / 'gcmap.js'
    
    argdict, vfdict, isdict, dfdict, metaldict, popdict, symdict = getblast(infile_name, rootdir, threads)
    
    ICEss = {}
    for ICEtag, mge in dictMGE.items():
        genelist = []
        regi = infile_name + '_' + ICEtag   # eg. aHPS7_ICE2
        regijs = ICEtag
        genefile = final_dir / f'{regi}_gene.json'
        infofile = final_dir / f'{regi}_info.json'
        gcjson = js_dir / f'{regijs}_gc.js'
        mapfile = js_dir / f'{regijs}.js'
        htmlfile = final_dir / f'{regi}.html'
        
        contig = find_located_contig(MGEdict[ICEtag].keys(), genome_info)
        
        
        #[myDR1, myDR2, myDR3, myDR4, start_point, end_point, finalstart, finalend, trnalist, record]
        
        # scan 5 gene before the ICE
        for mov in range(mge.finalstart, mge.start_point):
            gene = gene_key(prefix, mov)
            genelist.append(process_gene(
                gene, genome_info[contig]['posdict'], MGEdict, argdict, vfdict, isdict, dfdict, 
                metaldict, popdict, symdict, genome_info[contig]['locusdict'])
                )
        # genelist = [{'gene': locus_tag, 'pos': pos, 'prod': product, 'featu': feature}]
        # scan all gene within the ICE
        for mov in range(mge.start_point, mge.end_point + 1):
            gene = gene_key(prefix, mov)
            genelist.append(process_gene(
                gene, genome_info[contig]['posdict'], MGEdict, argdict, vfdict, isdict, dfdict, 
                metaldict, popdict, symdict, genome_info[contig]['locusdict'], feature_default = '', 
                is_ICE = True, ICEdict_tag = ICEtag)
                )
        
        # scan 5 gene after the ICE
        for mov in range(mge.end_point + 1, mge.finalend + 1):
            gene = gene_key(prefix, mov)
            genelist.append(process_gene(
                gene, genome_info[contig]['posdict'], MGEdict, argdict, vfdict, isdict, dfdict, 
                metaldict, popdict, symdict, genome_info[contig]['locusdict'])
                )
            
        with open(genefile, 'wt') as gene_file:   # output json file
            json.dump(genelist, gene_file, indent = 4)
            
        sgene = gene_key(prefix, mge.start_point)   # the first gene
        egene = gene_key(prefix, mge.end_point)   # the last gene
        s1, e1, strand1, pro1 = genome_info[contig]['posdict'][sgene]   # info of the first gene
        s2, e2, strand2, pro2 = genome_info[contig]['posdict'][egene]   # info of the last gene
        if mge.DR1 == '0':   # incase the DR is located in ori
            mge.DR1 = '1'
        gcc = gc(input_file, contig, int(mge.DR1), int(mge.DR4), filetype)   # gc content, eg. 37.86
        ICEss[regi] = ','.join([contig, mge.DR1, mge.DR4, str(mge.start_point), str(mge.end_point), 
                                gcc])   #startpos, endpos, sgene, egene, gc
        
        if mge.DR2: # pending if there is a attL
            DR1 = getfa(input_file, mge.DR1, mge.DR2, contig, filetype)   # attL
            DR2 = getfa(input_file, mge.DR3, mge.DR4, contig, filetype)   # attR
            DRw = 'attL:' + mge.DR1 + '..' + mge.DR2 + '(' + DR1 + ')  ' + 'attR:' + mge.DR3 + '..' + mge.DR4 + '(' + DR2 + ')'
            # DRw eg. attL:1353270..1353285(CGGATTTTGAATCCG)  attR:1421747..1421762(CGGATTTTGAATCCG)
        else:
            DRw = '-'
        if mge.trnalist:   # trnalist = [[start, end, strand, product]]
            trnaout = mge.trnalist[0][3] + ' (' + mge.trnalist[0][0] + '..' + mge.trnalist[0][1] + ') [' + mge.trnalist[0][2] + ']'
            # trnaout eg. tRNA-Leu (1421705..1421790) [+]
        else:
            trnaout = '-'
            
        oritseqs = oritseq(infile_name, regi, input_file, contig, mge.DR1, mge.DR4, filetype, rootdir, threads)
#		oritdesc = "<br>".join([oritseqs[i:i+63] for i in range(0, len(oritseqs), 63)])

        if 'IME' in regi:
            typeIE = 'IME'
        elif 'AICE' in regi:
            typeIE = 'AICE'
        else:
            typeIE = 'T4SS-type ICE'
            
        ICEinfo = {
			'Type' : typeIE,
			'Location (nt)' : mge.DR1 + '..' + mge.DR4,
			'Length (bp)' : str(int(mge.DR4) - int(mge.DR1) + 1),
			'GC Content (%)' : gcc,
			'oriT seq' : oritseqs,
			'DRs' : DRw,
			'Relaxase Type' : ','.join(infodict[ICEtag]['mob']),
			'Mating pair formation systems' : ','.join(infodict[ICEtag]['mpf']),
			'Close to tRNA' : trnaout
		}
        with open(infofile, 'wt') as info_file:
            json.dump(ICEinfo, info_file, indent = 4)   # output json file
            
        i = 1
        mapforlist = []   # forward strand
        maprevlist = []   # reverse strand
        for gene in genelist:   # genelist = [{'gene': locus_tag, 'pos': pos, 'prod': product, 'featu': feature}]
            color, name = getcolor(gene['featu'], gene['prod'])
            start = gene['pos'].split(' ')[0].split('..')[0]      # pos : eg. 1945..2222 [+], 278
            end = gene['pos'].split(' ')[0].split('..')[1]
            if gene['pos'].split('[')[1].split(']')[0] == '+':
                strand = 1
            else:
                strand = 0
            product = gene['prod']
            
            if product == '':
                product = 'hypothetical protein'
                
            anno = {
					'start' : start,
					'end' : end,
					'strand' : strand,
					'locus_tag' : 'M'+ str(i),
					'type' : 'others',
					'color' : color,
					'description' : 
                        'Location: ' + gene['pos'].split(' ')[0] + ' ('+gene['pos'].split(' ')[2] + 
                        ' bp)<br>Type: ' + name + '<br>Detail: ' + product
				}
            
            if strand == 1:   # pending if the gene located in forward or reverse strand
                mapforlist.append(anno)
            else:
                maprevlist.append(anno)
            i += 1
        
        head = 'var borders = [];\nvar tta_codons = [];\nvar orfs ='
        s = genelist[0]['pos'].split(' ')[0].split('..')[0]
        e = genelist[-1]['pos'].split(' ')[0].split('..')[1]
        
        gcdict = calculate_gc(input_file, contig, filetype, int(s), int(e))
        
        with open(gcmap, 'rt') as original_file:   # this is a template file to draw a gcmap
            original_content = original_file.read()
        with open(gcjson,'wt') as gein2:   # write to a file for draw gcmap for this genome
            gein2t = 'var jsonData = ' + str(gcdict)+';'
            gein2.write(gein2t)
            gein2.write(original_content)	
        
        maps = str(mapforlist)+';\nvar orfs2 ='+str(maprevlist)+';\nvar clusterf2 = { start: '+s+', end: '+ \
				  e+', idx: 1, orfs: orfs, borders: borders, tta_codons:tta_codons,\
				  label: \'\', unordered: true };\nvar clusterr2 = { start: '+ s+', end: '+ \
				  e+', idx: 2, orfs: orfs2, borders: borders, tta_codons:tta_codons,\
				  label: \'\', unordered: true };\nsvgene.drawClusters("'+regijs+'", [clusterf2, clusterr2], 50, 920);'
        
        with open(mapfile,'wt') as map_file:
            map_file.write(head + maps)
            
        with open(viewfile, 'rt') as file:
            file_content = file.read()
        new_content = file_content.replace('XXXX', regijs)
        with open(htmlfile, 'wt') as file:
            file.write(new_content)
            
    return ICEss

def get_color(region):

        coldict = {
            'T4SS-type ICE' : 'fill:rgba(0, 128, 164,0.9)',
            'IME':'fill:rgba(41,76,166,0.9)', 
            'AICE':'fill:rgba(255, 84,0,0.9)'
            }
        
        return coldict[region]

def getfasta(infile_name, input_file, ICE_name, contig, start_pos, end_pos, filetype, final_dir):
    
    if not final_dir.exists():
        final_dir.mkdir(parents = True)
        
    outfa = final_dir / f'{ICE_name}.fa'
    
    # the original writer haven't consider the filetype here, however, now draft genome still not be considered
    if filetype == 'fa':
        records = SeqIO.parse(input_file, 'fasta')
    else:
        records = SeqIO.parse(input_file, 'gb')
    
    with open(outfa, 'wt') as file:
        for record in records:
            if record.id == contig:
                sequence = str(record.seq[int(start_pos) - 1 : int(end_pos)])
        sequence_for_write = format_sequence(sequence)
        ICE_ID = f'>{infile_name}_{contig}_{start_pos}-{end_pos}'
        file.write(ICE_ID)
        file.write('\n')
        file.write(sequence_for_write)
        file.write('\n')
    
    return outfa, ICE_ID[1:]

def getrfasta(infile_name, input_file, Itag, ICE_name, ice, recovery_ICE, ICE_dict, genome_info,
              gbfile, out_dir, rootdir, threads):
    
    final_dir = workdir / out_dir / infile_name
    if not final_dir.exists():
        final_dir.mkdir(parents = True)
    outfa = final_dir / f'{ICE_name}.fa'
    for_gc = tmp_dir / f'{ICE_name}_gc.fa'
    
    records = SeqIO.parse(gbfile, 'genbank')
    genome_sequence = {}
    for record in records:
        genome_sequence[record.id] = record.seq
    # info: [DR1, DR2, DR3, DR4, final_left, final_right, end0_last, end1_last, trnalist, left_contig, right_contig]
    if ice.end0_last == 'P' and ice.end1_last == 'P':
        start_seq = str(genome_sequence[ice.left_contig])[int(ice.DR1) - 1:]
        start_ID = f'>{infile_name}_{ice.left_contig}_{int(ice.DR1) - 1}-{ice.len_left}'
        start_loc = f'{ice.left_contig}:{int(ice.DR1) - 1}..{ice.len_left}'
        end_seq = str(genome_sequence[ice.right_contig].reverse_complement())[:int(ice.DR4)]
        end_ID = f'>{infile_name}_rev_{ice.right_contig}_1-{ice.DR4}'
        end_loc = f'rev_{ice.right_contig}:1..{ice.DR4}'
        att_seq = str(start_seq[: int(ice.DR2) - int(ice.DR1) + 1]) if ice.DR2 else '-'

    elif ice.end0_last != 'P' and ice.end1_last != 'P':
        start_seq = str(genome_sequence[ice.left_contig].reverse_complement())[int(ice.DR1) - 1:]
        start_ID = f'>{infile_name}_rev_{ice.left_contig}_{int(ice.DR1) - 1}-{ice.len_left}'
        start_loc = f'rev_{ice.left_contig}:{int(ice.DR1) - 1}..{ice.len_left}'
        end_seq = str(genome_sequence[ice.right_contig])[:int(ice.DR4)]
        end_ID = f'>{infile_name}_{ice.right_contig}_1-{ice.DR4}'
        end_loc = f'{ice.right_contig}:1..{ice.DR4}'
        att_seq = str(start_seq[: int(ice.DR2) - int(ice.DR1) + 1]) if ice.DR2 else '-'

    else:
        if ice.end0_last == 'P':
            start_seq = str(genome_sequence[ice.left_contig])[int(ice.DR1) - 1:]
            start_ID = f'>{infile_name}_{ice.left_contig}_{int(ice.DR1) - 1}-{ice.len_left}'
            start_loc = f'{ice.left_contig}:{int(ice.DR1) - 1}..{ice.len_left}'
            end_seq = str(genome_sequence[ice.right_contig])[:int(ice.DR4)]
            end_ID = f'>{infile_name}_{ice.right_contig}_1-{ice.DR4}'
            end_loc = f'{ice.right_contig}:1..{ice.DR4}'
            att_seq = str(start_seq[: int(ice.DR2) - int(ice.DR1) + 1]) if ice.DR2 else '-'
        else:
            start_seq = str(genome_sequence[ice.right_contig])[int(ice.DR1) -1:]
            start_ID = f'>{infile_name}_{ice.right_contig}_{ice.DR1}-{ice.len_right}'
            start_loc = f'{ice.right_contig}:{ice.DR1}..{ice.len_right}'
            end_seq = str(genome_sequence[ice.left_contig])[:int(ice.DR4)]
            end_loc = f'{ice.left_contig}:1..{ice.DR4}'
            att_seq = str(start_seq[: int(ice.DR2) - int(ice.DR1) + 1]) if ice.DR2 else '-'

    ICE_location = [start_loc]
    ICE_length = len(start_seq) + len(end_seq)
    middle_contigs = '-'
    with open(outfa, 'wt') as file, open(for_gc, 'wt') as gc_file:
        sequence_for_write = format_sequence(start_seq)
        file.write(start_ID)
        file.write('\n')
        file.write(sequence_for_write)
        file.write('\n')
        
        gc_file.write(f'>{ICE_name}')
        gc_file.write('\n')
        gc_file.write(start_seq)
        
        if recovery_ICE[Itag]['middle']:
            for middle in recovery_ICE[Itag]['middle']:
                middle_contig = find_located_contig(ICE_dict[middle]['hit_id'].keys(), genome_info)
                ICE_location.append(middle_contig)
                middle_contigs.append(middle_contig)
                mid_seq = str(genome_sequence[middle_contig])
                ICE_length += len(mid_seq)
                mid_ID = f'>{infile_name}_{mid_seq}'
                sequence_for_write = format_sequence(mid_seq)
                file.write(mid_ID)
                file.write('\n')
                file.write(sequence_for_write)
                file.write('\n')
                
                gc_file.write(mid_seq)
        
        sequence_for_write = format_sequence(end_seq)
        file.write(end_ID)
        file.write('\n')
        file.write(sequence_for_write)
        file.write('\n')
        
        gc_file.write(end_seq)
        
    ICE_location.append(end_loc)
    ICE_gc = gc(for_gc, ICE_name, 1, ICE_length, 'fa')
    ori = oritseq(infile_name, ICE_name, for_gc, ICE_name, 1, ICE_length, 'fa', rootdir, threads)
    
    return outfa, ICE_location, str(ICE_length), ICE_gc, ori, att_seq, (';').join(middle_contigs) if middle_contigs else '-'

def generate_comp_output(infile_name, input_file, ICEss, filetype, genome_info, out_dir, workdir, rootdir):
    
    final_dir = workdir / out_dir / infile_name
    js_dir = workdir / out_dir / infile_name / 'js'
    jsback = rootdir / 'script' / 'js'
    
    i = 1 
    ICEsumlist = []
    homelist = []
    if ICEss:
        for key, value in ICEss.items():
            [contig, s, e, stag, etag, gc] = value.split(',')
            lengt = int(e) - int(s) + 1
            
            if 'IME' in key:
                typeIE = 'IME'
            elif 'AICE' in key:
                typeIE = 'AICE'
            else:
                typeIE = 'T4SS-type ICE'
            
            ICEs = {
				'region':'Region'+ str(i),
		        'location': s + '..' + e,
		        'length': str(lengt),
		        'gc':gc,
		        'type': typeIE,
		        'detail': key
		    }
            homedict = {
			  	'start' : int(s),
			   	'end' : int(e),
			   	'color' : get_color(typeIE),
			   	'info' : '_'.join([str(lengt), str(gc), typeIE]),
			   	'text' : 'Region'+ str(i)
			   }
            homelist.append(homedict)
            ICEsumlist.append(ICEs)
            i += 1
    
    ICEsum = final_dir / f'{infile_name}_ICEsum.json'
    with open(ICEsum,'wt') as ice_file:
        json.dump(ICEsumlist, ice_file, indent = 4)
    
    shutil.copytree(jsback, js_dir, dirs_exist_ok=True)

def generate_output(out_dir, dictMGE, infile_name, input_file, filetype, infodict, rdictICE, 
                    recovery_ICE, ICE_dict, rinfodict, genome_info, gbfile, rootdir, threads):
    
    final_dir = workdir / out_dir / infile_name
    out_detail = workdir / out_dir / 'ICE_details.tsv'
    out_summary = workdir / out_dir / 'ICE_summary.tsv'
    
    if not out_detail.exists():
        
        detail_header = ['Isolate', 'MGE', 'Location', 'Length', 'GC', 'Relaxase_Type', 
                         'Systems', 'oriT', 'attL', 'attR', 'att_seq', 'tRNA', 'AMR', 'middle', 
                         'middle_probability']
        summary_header = ['Isolate', 'ICE(rICE_included)', 'ICE(rICE_excluded)', 'IME', 'AICE', 'rICE']
        
        with open(out_detail, 'wt') as detail, open(out_summary, 'wt') as summary:
            detail.write(('\t'.join(detail_header)))
            detail.write('\n')
            summary.write(('\t'.join(summary_header)))
            summary.write('\n')
    
    out_detail_file = open(out_detail, 'at')
    out_summary_file = open(out_summary, 'at')
    
    ICE_count, IME_count, AICE_count = 0, 0, 0
    
    if dictMGE:
        for Itag, mge in dictMGE.items():
            if 'IME' in Itag:
                IME_count += 1
            elif 'AICE' in Itag:
                AICE_count += 1
            else:
                ICE_count += 1
                
            ICE_name = f'{infile_name}_{Itag}'
            outfa, ICE_ID = getfasta(infile_name, input_file, ICE_name, mge.record, int(mge.DR1), int(mge.DR4), filetype, final_dir)
            ICE_seq = SeqIO.read(outfa, 'fasta')
            ICE_length = str(len(ICE_seq))
            ICE_gc = gc(outfa, ICE_ID, 0, ICE_length, 'fa')
            relaxase = ','.join(infodict[Itag]['mob']) if infodict[Itag]['mob'] else '-'
            system = ','.join(infodict[Itag]['mpf']) if infodict[Itag]['mpf'] else '-'
            ori = oritseq(infile_name, ICE_name, outfa, ICE_ID, 1, ICE_length, 'fa', rootdir, threads)
            attL, attR, att_seq = '-', '-', '-'
            if mge.DR2:
                attL = f'{mge.DR1}..{mge.DR2}'
                attR = f'{mge.DR3}..{mge.DR4}'
                att_seq = str(ICE_seq.seq[0 : int(mge.DR2) - int(mge.DR1) + 1])
            if mge.trnalist:
                trna = f'{mge.trnalist[0][3]},{mge.record}:{mge.trnalist[0][0]}..{mge.trnalist[0][1]}[{mge.trnalist[0][2]}]'
            else:
                trna = '-'
            # trna eg. tRNA-Leu,contig1:1421705..1421790[+]
            _, amr_list = abricate(outfa, 'resfinder')
            ICE_amr = (';').join(amr_list) if amr_list else '-'
            middle_contig = '-'
            posibility = '-'
            
            detail_content = [infile_name, ICE_name, f'{mge.record}:({mge.DR1}..{mge.DR4})', ICE_length, 
                              ICE_gc, relaxase, system, ori, attL, attR, att_seq, trna, ICE_amr, 
                              middle_contig, posibility]
            
            out_detail_file.write(('\t'.join(detail_content)))
            out_detail_file.write('\n')
    
    rICE_count = 0
    if rdictICE:
        for Itag, ice in rdictICE.items():
            # info: [DR1, DR2, DR3, DR4, final_left, final_right, end0_last, end1_last, trnalist, left_contig, right_contig]
            rICE_count += 1
            ICE_name = f'{infile_name}_{Itag}'
            outfa, ICE_location, ICE_length, ICE_gc, ori, att_seq, middle_contigs = getrfasta(
                infile_name, input_file, Itag, ICE_name, ice, recovery_ICE, ICE_dict, genome_info, gbfile, 
                out_dir, rootdir, threads)
            relaxase = ','.join(rinfodict[Itag]['mob']) if rinfodict[Itag]['mob'] else '-'
            system = ','.join(rinfodict[Itag]['mpf']) if rinfodict[Itag]['mpf'] else '-'
            attL, attR = '-', '-'
            if ice.DR2:
                if ice.end0_last != 'P' and ice.end1_last == 'P':
                    attL = f'{ice.right_contig}:{ice.DR1}..{ice.DR2}'
                    attR = f'{ice.left_contig}:{ice.DR3}..{ice.DR4}'
                else:
                    attL = f'{ice.left_contig}:{ice.DR1}..{ice.DR2}'
                    attR = f'{ice.right_contig}:{ice.DR3}..{ice.DR4}'
            if ice.trnalist:   # eg. trnainfo = ['contig1', '18651', '18726', '+', 'tRNA-Glu(ttc)']
                trna = f'{ice.trnalist[0][4]},{ice.trnalist[0][0]}:{ice.trnalist[0][1]}..{ice.trnalist[0][2]}[{ice.trnalist[0][3]}]'
            else:
                trna = '-'
            _, amr_list = abricate(outfa, 'resfinder')
            ICE_amr = (';').join(amr_list) if amr_list else '-'        
            if recovery_ICE[Itag]['middle_importance']:
                posibility = recovery_ICE[Itag]['middle_importance']
            else:
                posibility = '-'
                
            detail_content = [infile_name, ICE_name, (',').join(ICE_location), ICE_length, 
                              ICE_gc, relaxase, system, ori, attL, attR, att_seq, trna, ICE_amr, 
                              middle_contigs, posibility]
            
            out_detail_file.write(('\t'.join(detail_content)))
            out_detail_file.write('\n')
            
    summary_content = [infile_name, str(ICE_count + rICE_count), str(ICE_count), str(IME_count), str(AICE_count), str(rICE_count)]
        
    out_summary_file.write(('\t'.join(summary_content)))
    out_summary_file.write('\n')
    
    out_detail_file.close()
    out_summary_file.close()
    
    print(f'{infile_name} : {str(ICE_count + rICE_count)} ICE, {str(IME_count)} IME')
    
    
def _single(infile_name, input_file, filetype, rootdir, out_dir, json_out, threads):
    # create the result dir with genome name, which is infile_name
    
    prefix = 'TMPID'
    # treat fasta and genbank differerntly to collect the gene details
    if  filetype == 'fa':
        prokkanno(infile_name, input_file, prefix, threads)
        genome_info, gbfile = getgff1(infile_name)
    else:
        genome_info, gbfile = getgff(infile_name, input_file, prefix)
    
    genome_info = infomation_filter_marker(genome_info, input_file)
    
    MGEdict, infodict = get_MGE(infile_name, input_file, genome_info, rootdir)
    recovery_ICE, ICE_dict, rICEdict, rinfodict = ICE_recovery(infile_name, genome_info)
    
    DRlist = get_DR(infile_name, input_file)
    
    MGEdict = remove_duplicate(MGEdict, rICEdict)
    
    dictMGE = {}
    if MGEdict:
        for MGE, info in MGEdict.items():
            merged = merge_tRNA(info, DRlist, prefix, genome_info)
            if merged is not None:
                dictMGE[MGE] = merged
    
    MGEdict = MGE_reorder(MGEdict)
    infodict = MGE_reorder(infodict)
    dictMGE = MGE_reorder(dictMGE)
    
    rdictICE = boundary_of_rICE(infile_name, gbfile, recovery_ICE, ICE_dict, genome_info, rICEdict, prefix)
    
    if dictMGE and json_out:
        ICEss = get_map(infile_name, input_file, MGEdict, infodict, genome_info, dictMGE, filetype, prefix, rootdir, out_dir, threads)
        generate_comp_output(infile_name, input_file, ICEss, filetype, genome_info, out_dir, workdir, rootdir)
    
    generate_output(out_dir, dictMGE, infile_name, input_file, filetype, infodict, rdictICE, 
                        recovery_ICE, ICE_dict, rinfodict, genome_info, gbfile, rootdir, threads)
    
    tmpfile = tmp_dir / infile_name
    shutil.rmtree(tmpfile)
