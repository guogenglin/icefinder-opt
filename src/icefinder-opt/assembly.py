# -*- coding: utf-8 -*-
"""
Created on Sun May 17 15:31:17 2026

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

from .utils import abricate, format_sequence, gene_key
from .log import log

def run_blast(faa_file: Path, database: str, threads: int, evalue = 1e-4) -> str:
    '''
    Run a BLASTP search against the specified database and return the output as a string.
    '''
    command = ['blastp', '-query', str(faa_file), '-db', database, '-evalue', str(evalue), '-num_threads', 
               str(threads), '-max_hsps', '1', '-num_descriptions', '1', '-num_alignments', '1', 
               '-outfmt', '6 std slen stitle']
    process = subprocess.run(command, capture_output = True, text = True, check = True)

    return process.stdout

def havalue(threshold: float, out: str) -> dict[str, str]:
    '''
    Calculate the havalue (coverage * identity) for each BLAST hit and filter based on the threshold.
    '''
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

def getdf(faa_file: Path, temp_dir: Path, threads: int) -> dict[str, str]:
    '''
    Run defense-finder and process the results to get the defense system information.
    '''
    #defense-finder used the same way of macsyfinder
    defcmd = ['defense-finder', 'run', '-w', str(threads), str(faa_file), '-o', temp_dir / 'def']
    
    process = subprocess.run(defcmd, check = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    
    def_file = temp_dir / 'def' / f'{faa_file.stem}_defense_finder_genes.tsv'
    dfdict = {}

    # replicon	hit_id	gene_name
    # eg. GCA_011044315	GCA_011044315_00151	HEC-02__HEC-02A

    with open(def_file) as defile:
        for line in defile:
            lines = line.strip().split('\t')
            if lines[0] and lines[0] != 'replicon':
                dfdict[lines[1]] = lines[2].replace('__',',')   # hit_id : gene_name
        
        return dfdict
# ─────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────
class Gene:
    '''
    Represents a gene in a genomic sequence.
    '''
    def __init__(self, gene_id : str, origin_id : str, start: int, end: int, strand: str, product: str, 
                 trna: str = 'F', feature: list[str] | None = None, product_map: list[str] | None = None):
        self.gene_id        = gene_id
        self.origin_id      = origin_id
        self.start          = start
        self.end            = end
        self.strand         = strand
        self.product        = product
        self.trna           = trna
        self.feature        = feature or []
        self.product_map    = product_map or []

    @property
    def locus_num(self):
        '''Return the locus number of the gene, which is the numeric part of the gene ID after the last underscore.'''
        return int(self.gene_id.rsplit('_', 1)[1])
    
    def get_args(self, assembly_info: Assembly, feature: str):
        '''
        Get the feature and product information for the gene based on the assembly information and the provided feature.
        '''
        # search if genes belong to these dict, if so, add correspond feature and product info
        features = [feature] if feature else []
        product_map = [self.product]
        
        dict_label_map  = {
            'AR' : assembly_info.argdict, 
            'VF' : assembly_info.vfdict, 
            'IS' : assembly_info.isdict, 
            'Defense' : assembly_info.dfdict, 
            'Metal' : assembly_info.metaldict, 
            'Degradation' : assembly_info.popdict, 
            'Symbiosis' : assembly_info.symdict, 
            }

        for label, dic in dict_label_map.items():
            if self.gene_id in dic:
                features.append(label)
                product_map.append(dic[self.gene_id])
        self.feature = '; '.join(list(filter(None, features)))
        self.product_map = '; '.join(list(filter(None, product_map))).replace('hypothetical protein;', '').replace(';hypothetical protein', '')


    
class Contig:
    '''
    Represents a contig in a genomic assembly.
    '''
    def __init__(
            self, contig_name: str, genes: dict[str, Gene] | None = None, length: int | None = None, 
            plasmid: str = 'F'):
        self.name           = contig_name
        self.genes          = genes or {}
        self.length         = length or int()
        self.plasmid        = plasmid


class Assembly:
    '''
    Represents a genomic assembly.
    '''
    def __init__(self, sample_path: Path, sample_name: str, file_type: str, temp_id: str = 'TMPID', 
                 gff_file: Path | None = None, gbk_file: Path | None = None, faa_file: Path | None = None, 
                 fna_file: Path | None = None, ICE_dir: Path | None = None, 
                 contigs: dict[str, Contig] | None = None, argdict: dict[str, str] | None = None, 
                 vfdict: dict[str, str] | None = None, isdict: dict[str, str] | None = None, 
                 dfdict: dict[str, str] | None = None, metadict: dict[str, str] | None = None, 
                 popdict: dict[str, str] | None = None, symdict: dict[str, str] | None = None):
        self.sample_path        = sample_path
        self.sample_name        = sample_name
        self.file_type          = file_type
        self.temp_id             = temp_id
        self.gff_file           = gff_file or Path()
        self.gbk_file           = gbk_file or Path()
        self.faa_file           = faa_file or Path()
        self.fna_file           = fna_file or Path()
        self.ICE_dir            = ICE_dir or Path()
        self.contigs            = contigs or {}
        self.argdict            = argdict or {}
        self.vfdict             = vfdict or {}
        self.isdict             = isdict or {}
        self.dfdict             = dfdict or {}
        self.metadict           = metadict or {}
        self.popdict            = popdict or {}
        self.symdict            = symdict or {}

    def get_gene_info_from_gff(self, verbose: bool = False):
        '''
        Extract gene information from a GFF file and populate the assembly information.
        '''
        with open(self.gff_file) as gff_handle:
            for line in gff_handle:
                if 'ID=' not in line:
                    continue  # skip lines without ID immediately
                
                # GFF (General Feature Format) field description example
                # CP001321.1	Prodigal:v2.6	CDS	2	1021	.	+	0	ID=GCA_000021885_00001;eC_number=1.2.1.12;Name=gapA;db_xref=COG:COG0057;gene=gapA;inference=ab initio prediction:Prodigal:v2.6,similar to AA sequence:UniProtKB:P0A9B2;locus_tag=GCA_000021885_00001;product=Glyceraldehyde-3-phosphate dehydrogenase A
                # 0. seqid      : Sequence ID / chromosome / contig name
                # 1. source     : Annotation source or software
                # 2. type       : Feature type
                # 3. start      : Start position on sequence (1-based)
                # 4. end        : End position on sequence
                # 5. score      : Feature score, '.' means no score available
                # 6. strand     : Strand information
                # 7. phase      : Reading frame phase for CDS
                # 8. attributes : Semicolon-separated key=value annotations
                content = line.strip().split('\t')
                if content[2] not in ['CDS', 'rRNA', 'tRNA', 'tmRNA']:
                        continue
                contig = content[0]
                if contig not in self.contigs:
                    self.contigs[contig] = Contig(contig)

                attributes = dict(item.split('=') for item in content[8].split(';') if '=' in item)
                # Attribute explanation:
                # ID           : Unique feature identifier
                # eC_number    : Enzyme Commission number
                # Name         : Feature name
                # db_xref      : Cross-reference to external database
                # gene         : Gene symbol/name
                # inference    : Annotation evidence source
                # locus_tag    : Stable locus identifier
                # product      : Predicted protein product/function
                ids = attributes.get('ID')
                if not ids:
                    continue  # skip if no ID
                product = attributes.get('product', '')  # default empty if no product
                self.contigs[contig].genes[ids] = (Gene(ids, ids, int(content[3]) - 1, int(content[4]), 
                                                                content[6], product))

                if content[2] in ['tRNA', 'tmRNA']:
                    self.contigs[contig].genes[ids].trna = 'T'

        for record in SeqIO.parse(self.sample_path, 'fasta'):
            self.contigs[record.id].length = len(record.seq)
        
        log('Gene information extracted from GFF file.', verbose = verbose)

    def get_gene_info_from_gbk(self, temp_dir: Path, verbose: bool = False):
        '''
        Extract gene information from a GenBank file and populate the assembly information.
        '''
        Path(temp_dir / self.sample_name).mkdir(parents = True, exist_ok = True)
        
        self.faa_file    = temp_dir / self.sample_name / f'{self.sample_name}.faa'
        self.fna_file    = temp_dir / self.sample_name / f'{self.sample_name}.ffn'
        self.gbk_file    = self.sample_path
        self.sample_path = temp_dir / self.sample_name / f'{self.sample_name}.fasta'

        i = 1

        with open(self.faa_file, 'w') as faa, open(self.fna_file, 'w') as fna, open(self.sample_path, 'w') as fasta:
            for record in SeqIO.parse(self.gbk_file, 'genbank'):
                SeqIO.write(record, fasta, 'fasta')
                self.contigs[record.id] = Contig(record.id)

                for feature in record.features:
                    if feature.type not in ['CDS', 'rRNA', 'tRNA', 'tmRNA']:
                        continue
                    if 'locus_tag' not in feature.qualifiers:
                        continue
                    locus_tag = feature.qualifiers['locus_tag'][0]

                    # Determine start, end, and strand
                    ori = ''
                    # if compoundlocation exists, find the start and end, or use freature.location
                    # a gene may across the ori of the bacteria sequencing, the annotation may like this
                    # (999800..1000000,1..500)
                    if isinstance(feature.location, CompoundLocation):
                        loc_list = sorted(feature.location.parts, key=lambda part: part.start)
                        loc = loc_list[-1]
                        if int(loc.end) == len(record.seq):
                            start = int(loc_list[-1].start) - 1
                            end = int(loc_list[0].end)
                            ori = 'T'
                        else:
                            start = int(loc.start) - 1
                            end = int(loc.end)
                    else:
                        loc = feature.location
                        start = int(loc.start) - 1
                        end = int(loc.end)
                    #.strand is a number, -1 or 1, switch to + or -
                    strand = '+' if loc.strand == 1 else '-'

                    # Determine product
                    if feature.type in ['CDS', 'rRNA']:
                        product = feature.qualifiers.get('product', ['-'])[0]
                    elif feature.type == 'tRNA':
                        product = feature.qualifiers.get('product', ['-'])[0]
                    else:  # tmRNA
                        product = 'tmRNA'

                    newid = gene_key(self.temp_id, i)
                    self.contigs[record.id].genes[newid] = (Gene(newid, locus_tag, start, end, strand, 
                                                                            product))

                    # Write sequences for CDS/rRNA
                    if feature.type in ['CDS', 'rRNA'] and 'translation' in feature.qualifiers:
                        aa_sequence = str(feature.qualifiers['translation'][0])
                        faa.write(f'>{newid} {product}\n')
                        faa.write(f'{format_sequence(aa_sequence)}\n')

                        # write new id, information, and na sequence into fna
                        if not ori:
                            cds_sequence = record.seq[feature.location.start:feature.location.end]
                        else:
                            cds_sequence = record.seq[start:len(record.seq)] + record.seq[0:end]
                        if strand == '-':
                            cds_sequence = cds_sequence.reverse_complement()
                        fna.write(f'>{newid} {product}\n')
                        fna.write(f'{format_sequence(str(cds_sequence))}\n')

                    # Record tRNA
                    if feature.type == 'tRNA' or feature.type == 'tmRNA':
                        self.contigs[record.id].genes[newid].trna = 'T'

                    i += 1
        
        for record in SeqIO.parse(self.sample_path, 'fasta'):
            self.contigs[record.id].length = len(record.seq)

        log('Gene information extracted from GenBank file.', verbose = verbose)

    def information_filter_marker(self, temp_dir: Path, threads: int, verbose: bool = False):
        '''
        Filter out contigs with no genes and identify plasmid contigs.
        '''
        # Run abricate plasmidfinder
        command = ['mob_recon', '-i', str(self.sample_path), '-o', str(temp_dir / 'mob_recon_results'), 
                '-n', str(threads)]
        subprocess.run(command, check = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
        
        with open(temp_dir / 'mob_recon_results' / 'contig_report.txt') as file:
            next(file)
            for line in file:
                content = line.strip().split('\t')
                if content[1] != 'plasmid':
                    continue
                contig_id = content[4].split()[0]
                self.contigs[contig_id].plasmid = 'T'

        remove_ids = []
        for contig_id, contig in self.contigs.items():
            if not contig.genes:
                remove_ids.append(contig_id)
        for contig_id in remove_ids:
            self.contigs.pop(contig_id)

        log('Plasmid contigs identified.', verbose = verbose)

    def getblast(self, rootdir: Path, temp_dir: Path, threads: int, verbose: bool = False):
        '''
        Run BLAST searches against various databases and process the results to populate the assembly information.
        '''
        # a dict set the config of every blast
        blast_config = {
            'IS':    {'db': str(rootdir / 'data' / 'transposase')},
            'METAL': {'db': str(rootdir / 'data' / 'metal')},
            'POP':   {'db': str(rootdir / 'data' / 'degradation')},
            'SYM':   {'db': str(rootdir / 'data' / 'symbiosis')},
            }
        
        # run blast and get output
        for db_name, config in blast_config.items():
            config['out'] = run_blast(self.faa_file, config['db'], threads)
        
        # pending the result and collect the match
        self.isdict = havalue(0.64, blast_config['IS']['out'])
        self.vfdict, _ = abricate(self.fna_file, 'vfdb')
        self.argdict, _ = abricate(self.fna_file, 'resfinder')
        self.metaldict = havalue(0.64, blast_config['METAL']['out'])
        self.popdict = havalue(0.64, blast_config['POP']['out'])
        self.symdict = havalue(0.64, blast_config['SYM']['out'])
        
        self.dfdict = getdf(self.faa_file, temp_dir, threads)

        log('Blast results processed.', verbose = verbose)



