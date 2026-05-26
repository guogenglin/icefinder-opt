# -*- coding: utf-8 -*-
"""
Created on Sun May 17 15:31:17 2026

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

from .utils import gene_key, getfa, oritseq, get_DR
from .assembly import Assembly, Contig, Gene


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

AICE_genes = {
    'FtsK_SpoIIIE' : ['FtsK_SpoIIIE'], 
    'RepSAv2' : ['RepSAv2', 'DUF3631', 'Prim-Pol'], 
    'Phage_integrase' : ['Recombinase']
}

def get_feature(feature: str) -> str:
    '''
    Map a feature to its category.
    '''
    featuredict = {
		'Phage_integrase':'Integrase', 'UPF0236':'Integrase',
        'Recombinase':'Integrase', 'rve':'Integrase',
        'TIGR02224':'Integrase', 'TIGR02249':'Integrase',
        'TIGR02225':'Integrase', 'PB001819':'Integrase',
		'RepSAv2':'Rep', 'DUF3631':'Rep', 'Prim-Pol':'Rep',
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
        if key in feature:   # eg. T4SS_G_tfc24, T4SS_MOBH
            parts = feature.split('_', 1)
            tag = parts[1] if len(parts) > 1 else feature
            return f'{category}@{tag}'
    
    return f'T4SS@{feature}'   # this return may never be used


def find_max_distance(loc_region: list) -> tuple[int, int]:
    '''
    Find the maximum distance between consecutive elements in a list.
    '''
    max_distance = -1
    max_distance_index = -1
    
    for i in range(len(loc_region) - 1):
        distance = abs(loc_region[i] - loc_region[i + 1])
        if distance > max_distance:
            max_distance = distance
            max_distance_index = i
    
    return loc_region[max_distance_index], loc_region[max_distance_index + 1]


def pos_tag(DR_pos: str, gene_dict: dict[str, Gene], edge_point: int, expanded_edge: int, max_in_contig: int, 
            dirtag: str, EXPAND: int = 20) -> tuple[int, int]:
    '''
    Tag the position of a gene in relation to the direct repeats (DRs) of an ICE.
    '''
    for gene in gene_dict.values():   # the dict must be sorted increasingly.

        if gene.end >= int(DR_pos):   # find which gene is in the end of ICE. if end > DR
            if dirtag == 's':
                return gene.locus_num, max(max_in_contig, edge_point - EXPAND)
            else:
                if gene.start > int(DR_pos):   # if the start of the gene also < DR
                    edge_point = gene.locus_num - 1   # which means the last gene is the end
                else:
                    edge_point = gene.locus_num   # or, the current gene is the end
                return edge_point, min(max_in_contig, edge_point + EXPAND)
            
    return edge_point, expanded_edge


# ─────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────
class MacSyHit:
    '''
    Represents a hit from the MacSyFinder output.
    '''
    def __init__(self, hit_id: str, gene_name: str | None = None, tags: str | None = None, 
                 loc_contig: str | None = None, dis_to_boundary: int = int()):
        self.hit_id          = hit_id or ''
        self.gene_name       = gene_name or ''
        self.tags            = tags or ''
        self.loc_contig      = loc_contig or ''
        self.dis_to_boundary = dis_to_boundary or int()

    @classmethod
    def from_macsy_line(cls, hit_id: str, gene_name: str, assembly_contigs: dict[str, Contig]):
        '''
        Create a MacSyHit instance from a line of MacSyFinder output.
        '''
        loc_contigs, dis_to_boundary = '', 0
        for contig_id, contigs in assembly_contigs.items():
            if hit_id in contigs.genes.keys():
                loc_contigs = contig_id
                gene_id_loc = [int(gene.gene_id.rsplit('_', 1)[1]) for gene in contigs.genes.values()]

                dis_to_boundary = min(
                    abs(int(hit_id.rsplit('_', 1)[1]) - min(gene_id_loc)),
                    abs(int(hit_id.rsplit('_', 1)[1]) - max(gene_id_loc))
                )

        return(cls(hit_id, gene_name, get_feature(gene_name), loc_contigs, dis_to_boundary))
    
    @property
    def gene_family(self):
        '''Determine the family of the gene based on its name.'''
        for family, genes in mandatory_genes.items():
            if self.gene_name in genes:
                return family
        for family, genes in accessory_genes.items():
            if self.gene_name in genes:
                return family
        for family, genes in AICE_genes.items():
            if self.gene_name in genes:
                return family
        return 'neutral'

    @property
    def locus_num(self):
        '''Get the locus number from the hit ID.'''
        return int(self.hit_id.rsplit('_', 1)[1])


class TypingResults:
    '''
    Represents the results of a typing operation.
    '''
    def __init__(self, sys_id: str, category: str | None, mob: set | None = None, mpf: set | None = None, 
                 elements: list[MacSyHit] | None = None, contig: str | None = None, 
                 start_point: int | None = None, end_point: int | None = None, 
                 trnalist: list[Gene] | None = None, expanded_start: int | None = None, 
                 expanded_end: int | None = None, DR1: int | None = int(), DR2: int | None = int(), 
                 DR3: int | None = int(), DR4: int | None = int(), DR_seq: str | None = None, 
                 DRw: str | None = '-', gc: str | None = None, trnaout: str | None = '-', 
                 oritseqs: str | None = '-', length: int | None = None):
        self.sys_id           = sys_id or ''
        self.category         = category or ''
        self.mob              = mob or set()
        self.mpf              = mpf or set()
        self.elements         = elements or []
        self.contig           = contig or ''
        self.start_point      = start_point
        self.end_point        = end_point
        self.trnalist         = trnalist or []
        self.expanded_start   = expanded_start
        self.expanded_end     = expanded_end
        self.DR1              = DR1 or int()
        self.DR2              = DR2 or int()
        self.DR3              = DR3 or int()
        self.DR4              = DR4 or int()
        self.DR_seq           = DR_seq or ''
        self.DRw              = DRw or '-'
        self.gc               = gc or ''
        self.trnaout          = trnaout or '-'
        self.oritseqs         = oritseqs or '-'
        self.length           = length or int()


    def add_hit(self, gene_name: str, model_fqn: str, hit_id: str, assembly_contigs: dict[str, Contig]):
        '''
        Add a new hit to the typing results.
        '''
        tags = get_feature(gene_name)  # lines[2] is gene_name; tags = {category}@{tag}

        # Extract relaxase information
        if 'Relaxase@' in tags:
            self.mob.add(tags.split('@')[1])
        # Extract T4SS type if present
        if 'T4SS' in model_fqn:
            self.mpf.add(model_fqn.rsplit('/', 1)[-1].split('_')[1])   # eg. ICEscan/Chromosome/T4SS_typeG
        
        self.elements.append(MacSyHit.from_macsy_line(hit_id, gene_name, assembly_contigs))

    def merge_tRNA(self, assembly_info: Assembly, temp_dir: Path, EXPAND: int = 20):
        '''
        Merge tRNA information into the typing results.
        '''
        # ori problem could be solved by recovery, so we don't need to consider it here.
        # The origin authors defined that if a tRNA gene is located near the predicted ICE region, it would be 
        # considered as a new boundary. However, they set the threshold at only five genes beyond the 
        # boundary, which seems insufficient. In their original version, a region of 250 kb around the 
        # integrase was considered. Although the algorithm differs, based on my personal experience, 
        # using only five genes is not enough. The current algorithm relies on ConjScan as its core, 
        # which only detects T4SS-related structures, while many other types of genes can also be present 
        # within an ICE.
        # In one of my cases, a tRNA gene was found 16 genes away from the predicted ICE, and its 
        # excision/integration form was experimentally verified by PCR. Therefore, I will temporarily 
        # set the threshold to 20 genes to test the accuracy.
        
        
        # Sometimes, at this step, an ICE may be predicted to be located on a broken contig. However, 
        # this is incorrect, because the input for macsyFinder is an FAA file, which does not record 
        # contig information and instead concatenates proteins from different contigs in a disordered 
        # manner. Therefore, such results need to be filtered out. 
        # If an ICE truly lies on a broken contig, it will be detected later during the recovery step.

        # Meanwhile, the ICE/IMEs were predicted in contigs which were predicted as plasmid will be filtered out.
        
        self.start_point = self.elements[0].locus_num
        self.end_point   = self.elements[-1].locus_num

        max_in_contig = max(gene.locus_num for gene in assembly_info.contigs[self.contig].genes.values())
        min_in_contig = min(gene.locus_num for gene in assembly_info.contigs[self.contig].genes.values())

        self.expanded_start = max(min_in_contig, self.start_point - EXPAND)
        self.expanded_end = min(max_in_contig, self.end_point + EXPAND)
        # We set 20 as expand, there maybe extra unexpected trna within
        candidate_left_trnas = [gene for gene in assembly_info.contigs[self.contig].genes.values() 
                           if gene.trna == 'T' and self.expanded_start <= gene.locus_num <= self.start_point]

        left_trna = min(candidate_left_trnas, key = lambda g: abs(g.locus_num - self.start_point), default = None)
        if left_trna:
            self.trnalist.append(left_trna)
            self.expanded_start = left_trna.locus_num
        candidate_right_trnas = [gene for gene in assembly_info.contigs[self.contig].genes.values() 
                           if gene.trna == 'T' and self.end_point <= gene.locus_num <= self.expanded_end]

        right_trna = min(candidate_right_trnas, key = lambda g: abs(g.locus_num - self.end_point), default = None)
        if right_trna:
            self.trnalist.append(right_trna)
            self.expanded_end = right_trna.locus_num
        # I found a special circumstances, that there are ICE contain more than one trna, one at the start and one 
        # at the end, so the original logic need to be modified.

        # here has a prerequisite, which is trna do not located within the ICE
        # the background is: normally ICE will insert into the genome nearby the tRNA
        self.DR1 = assembly_info.contigs[self.contig].genes[gene_key(assembly_info.temp_id, self.start_point)].start
        self.DR4 = assembly_info.contigs[self.contig].genes[gene_key(assembly_info.temp_id, self.end_point)].end

        if self.trnalist:
            seq_contig = Seq('')
            for record in SeqIO.parse(assembly_info.sample_path, 'fasta'):
                if record.id == self.contig:
                    seq_contig = record.seq
                    break
            temp_DR = temp_dir / 'detect_DR.fasta'

            with open(temp_DR , 'w') as handle:
                handle.write('>temp_DR\n')
                handle.write(str(seq_contig))

            DRlist = get_DR(temp_DR, temp_dir)   # in case the direct of contig is wrong

            left_gene = assembly_info.contigs[self.contig].genes[gene_key(assembly_info.temp_id, self.expanded_start)]
            right_gene = assembly_info.contigs[self.contig].genes[gene_key(assembly_info.temp_id, self.expanded_end)]

            for DR in DRlist:   # here is another prerequisite, the DR within trna must be the attL/R of ICE
                DRs = DR.split(',')   # left_start, left_end, right_start, right_end
                if left_gene.start <= int(DRs[0]) <= left_gene.end:
                    self.end_point, self.expanded_end = pos_tag(
                             DRs[3], assembly_info.contigs[self.contig].genes, self.end_point, self.expanded_end, max_in_contig, 'e', EXPAND)
                    self.start_point = self.expanded_start
                    self.expanded_start = max(self.expanded_start - EXPAND, min_in_contig)
                    self.DR1 = int(DRs[0]) + 1
                    self.DR2 = int(DRs[1])
                    self.DR3 = int(DRs[2]) + 1
                    self.DR4 = int(DRs[3])
                    break
                elif right_gene.start <= int(DRs[3]) <= right_gene.end:
                    self.start_point, self.expanded_start = pos_tag(
                            DRs[0], assembly_info.contigs[self.contig].genes, self.start_point, self.expanded_start, min_in_contig, 's', EXPAND)
                    self.end_point = self.expanded_end
                    self.expanded_end = min(self.expanded_end + EXPAND, max_in_contig)
                    self.DR1 = int(DRs[0]) + 1
                    self.DR2 = int(DRs[1])
                    self.DR3 = int(DRs[2]) + 1
                    self.DR4 = int(DRs[3])
                    break

    def basic_info(self, assembly_info: Assembly, rootdir: Path, temp_dir: Path, threads: int):
        '''
        Extract basic information about the predicted ICE, including sequence, GC content, and direct repeats.
        '''
        self.ICE_name = f'{assembly_info.sample_name}_{self.sys_id}'   # eg. aHPS7_ICE1
        self.ICE_ID = f'>{assembly_info.sample_name}_{self.contig}_{self.DR1}-{self.DR4}'
        self.seq = getfa(assembly_info.sample_path, self.contig, self.DR1, self.DR4)
        self.gc = f'{gc_fraction(self.seq) * 100:.2f}'
        if self.DR2: # pending if there is a attL
            self.DR_seq = str(getfa(assembly_info.sample_path, self.contig, self.DR1, self.DR2))   # attL
            # DRr_seq = str(getfa(assembly_info.sample_path, self.contig, self.DR3, self.DR4))   # attR
            self.DRw = f'attL:{self.DR1}..{self.DR2}({self.DR_seq})  self:{self.DR3}..{self.DR4}({self.DR_seq})'
            # DRw eg. attL:1353270..1353285(CGGATTTTGAATCCG)  attR:1421747..1421762(CGGATTTTGAATCCG)

        if self.trnalist:
            trnaout = []
            for trna in self.trnalist:
                trnaout.append(f'{trna.product} ({trna.start}..{trna.end}) [{trna.strand}]')
            self.trnaout = (';').join(trnaout)
            # trnaout eg. [tRNA-Leu (1421705..1421790) [+]; ]

        self.oritseqs = oritseq(assembly_info.sample_path, self.ICE_name, self.contig, self.DR1, self.DR4, rootdir, temp_dir, threads)
        # oritdesc = "<br>".join([oritseqs[i:i+63] for i in range(0, len(oritseqs), 63)])
        self.length = len(self.seq)