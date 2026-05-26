# -*- coding: utf-8 -*-
"""
Created on Wes May 20 13:19:31 2026

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

from .utils import get_DR, getfa, oritseq
from .alignment import MacSyHit, get_feature
from .assembly import Assembly, Gene
from .log import log, ICEfinderError


def is_connected(a: MacSyHit, b: MacSyHit, max_space: int) -> bool:
    '''
    Determine if two hits are connected based on their contig and locus information, and the distance to the contig boundary.
    '''
    # same contig
    if a.loc_contig == b.loc_contig:
        return abs(a.locus_num - b.locus_num) <= max_space

    # cross contig boundary rule
    return (a.dis_to_boundary + b.dis_to_boundary <= max_space)

def validate_cluster(cluster: list[MacSyHit], MGE_type: str, min_genes_required: int
                     ) -> bool:
    '''
    Validate a cluster of hits based on the MGE type and the required number of genes and families.
    '''
    if len(cluster) < min_genes_required:
        return False

    families_found = set()
    for hit in cluster:
        families_found.add(hit.gene_family)
    if MGE_type == 'IME':
        return set(['T4SS_MOBV', 'Phage_integrase']).issubset(families_found)
    elif MGE_type == 'AICE':
        return set(['FtsK_SpoIIIE', 'Phage_integrase']).issubset(families_found)
    else:
        return set(['T4SS_virb4', 'T4SS_t4cp1', 'T4SS_MOBV', 'Phage_integrase']).issubset(families_found)

def direction_pending(edge_contigs: dict[str, list[MacSyHit]], assembly_info: Assembly, EXPAND: int = 20
                      ) -> tuple[dict[str, list[MacSyHit]], str, Gene, dict[str, list[MacSyHit]], str, Gene]:
    '''
    Determine the direction of the edge contigs and the edge genes based on the distance to the target region.
    '''
    a_contig, b_contig = [{k : v} for k, v in edge_contigs.items()]

    a_target = [gene.locus_num for gene in assembly_info.contigs[next(iter(a_contig))].genes.values()]
    b_target = [gene.locus_num for gene in assembly_info.contigs[next(iter(b_contig))].genes.values()]

    a_locus = [hit.locus_num for hit in next(iter(a_contig.values()))]
    b_locus = [hit.locus_num for hit in next(iter(b_contig.values()))]

    if max(a_target) - max(a_locus) > EXPAND and min(b_locus) - min(b_target) > EXPAND:
        a_gene = next(gene for gene in assembly_info.contigs[next(iter(a_contig))].genes.values() if gene.locus_num == max(a_locus))
        b_gene = next(gene for gene in assembly_info.contigs[next(iter(b_contig))].genes.values() if gene.locus_num == min(b_locus))
        return b_contig, '+', b_gene, a_contig, '+', a_gene
    elif max(b_target) - max(b_locus) > EXPAND and min(a_locus) - min(a_target) > EXPAND:
        a_gene = next(gene for gene in assembly_info.contigs[next(iter(a_contig))].genes.values() if gene.locus_num == min(a_locus))
        b_gene = next(gene for gene in assembly_info.contigs[next(iter(b_contig))].genes.values() if gene.locus_num == max(b_locus))
        return a_contig, '+', a_gene, b_contig, '+', b_gene
    elif max(a_target) - max(a_locus) > EXPAND and max(b_target) - max(b_locus) > EXPAND:
        a_gene = next(gene for gene in assembly_info.contigs[next(iter(a_contig))].genes.values() if gene.locus_num == max(a_locus))
        b_gene = next(gene for gene in assembly_info.contigs[next(iter(b_contig))].genes.values() if gene.locus_num == max(b_locus))
        return a_contig, '-', a_gene, b_contig, '+', b_gene
    elif min(a_locus) - min(a_target) > EXPAND and min(b_locus) - min(b_target) > EXPAND:
        a_gene = next(gene for gene in assembly_info.contigs[next(iter(a_contig))].genes.values() if gene.locus_num == min(a_locus))
        b_gene = next(gene for gene in assembly_info.contigs[next(iter(b_contig))].genes.values() if gene.locus_num == min(b_locus))
        return a_contig, '+', a_gene, b_contig, '-', b_gene
    else:
        raise ICEfinderError('No valid direction found')
    
def search_trna(edge_contig: dict[str, list[MacSyHit]], edge_gene: Gene, edge_direction: str, 
                assembly_info: Assembly, flanking: str, EXPAND: int = 20) -> Gene:
    '''
    Search for tRNA genes in the flanking regions of the ICE and return the list of tRNA genes and the final edge gene after merging with tRNA.
    '''
    trnalist = []

    max_in_contig = max(gene.locus_num for gene in assembly_info.contigs[next(iter(edge_contig))].genes.values())
    min_in_contig = min(gene.locus_num for gene in assembly_info.contigs[next(iter(edge_contig))].genes.values())
    
    if flanking == 'left':
        if edge_direction == '+':
            expanded = max(min_in_contig, edge_gene.locus_num - EXPAND)
            for gene in assembly_info.contigs[next(iter(edge_contig))].genes.values():
                if not gene.trna == 'T':
                    continue
                if expanded <= gene.locus_num <= edge_gene.locus_num:
                    trnalist.append(gene)
            return min(trnalist, key = lambda g: abs(g.locus_num - edge_gene.locus_num), default = None)
            
        else:
            expanded = min(max_in_contig, edge_gene.locus_num + EXPAND)
            for gene in assembly_info.contigs[next(iter(edge_contig))].genes.values():
                if not gene.trna == 'T':
                    continue
                if edge_gene.locus_num <= gene.locus_num <= expanded:
                    trnalist.append(gene)
            return min(trnalist, key = lambda g: abs(g.locus_num - edge_gene.locus_num), default = None)
    elif flanking == 'right':
        if edge_direction == '-':
            expanded = max(min_in_contig, edge_gene.locus_num - EXPAND)
            for gene in assembly_info.contigs[next(iter(edge_contig))].genes.values():
                if not gene.trna == 'T':
                    continue
                if expanded <= gene.locus_num <= edge_gene.locus_num:
                    trnalist.append(gene)
            return min(trnalist, key = lambda g: abs(g.locus_num - edge_gene.locus_num), default = None)
        else:
            expanded = min(max_in_contig, edge_gene.locus_num + EXPAND)
            for gene in assembly_info.contigs[next(iter(edge_contig))].genes.values():
                if not gene.trna == 'T':
                    continue
                if edge_gene.locus_num <= gene.locus_num <= expanded:
                    trnalist.append(gene)
            return min(trnalist, key = lambda g: abs(g.locus_num - edge_gene.locus_num), default = None)
    else:
         raise ICEfinderError('flanking mark should be either "left" or "right".')

# ─────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────
class UnionFind:
    '''Union-Find data structure for clustering hits.'''
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        px, py = self.find(x), self.find(y)

        if px == py:
            return

        if self.rank[px] < self.rank[py]:
            self.parent[px] = py
        elif self.rank[px] > self.rank[py]:
            self.parent[py] = px
        else:
            self.parent[py] = px
            self.rank[px] += 1


class Recovery:
    def __init__(self, middle_contig: dict[str, list[MacSyHit]] | None = None, left_contig: dict[str, list[MacSyHit]] | None = None, 
                 right_contig: dict[str, list[MacSyHit]] | None = None, left_direction: str | None = None, 
                 right_direction: str | None = None, left_gene: Gene | None = None, right_gene: Gene | None = None, 
                 mob: set | None = None, mpf: set | None = None, all_gene_num: set[int] | None = None, 
                 trnalist: list[Gene] | None = [], expanded_left: Gene | None = None, 
                 expanded_right: Gene | None = None, DR1: int | None = 0, DR2: int | None = 0, DR3: int | None = 0, 
                 DR4: int | None = 0, DR_seq: str | None = None, left_seq: Seq | None = None, 
                 right_seq: Seq | None = None, ICE_location: list | None = None, gc: str | None = None, 
                 trnaout: str | None = '-', oritseqs: str | None = '-', length: int | None = None, 
                 all_middle: list | None = None):
        
        # left_gene....left_end | middle_contig | right_end....right_gene
        self.middle_contig    = middle_contig or {}
        # acturally we cannot determine which one is left or right of the MGE from the recovery results, 
        # the reason i write the object this way is to make it easier to understand the structure
        self.left_contig      = left_contig or {}
        self.right_contig     = right_contig or {}
        self.left_direction   = left_direction
        self.right_direction  = right_direction
        self.left_gene        = left_gene   # means the edge gene in the left contig
        self.right_gene       = right_gene
        self.mob              = mob or set()
        self.mpf              = mpf or set()
        self.all_gene_num     = all_gene_num or set()
        self.trnalist         = trnalist or []
        self.expanded_left    = expanded_left
        self.expanded_right   = expanded_right
        self.DR1              = DR1 or 0
        self.DR2              = DR2 or 0
        self.DR3              = DR3 or 0
        self.DR4              = DR4 or 0
        self.DR_seq           = DR_seq or '-'
        self.left_seq         = left_seq or Seq('')
        self.right_seq        = right_seq or Seq('')
        self.ICE_location     = ICE_location or []
        self.gc               = gc or ''
        self.trnaout          = trnaout or '-'
        self.oritseqs         = oritseqs or '-'
        self.all_middle       = all_middle or ['-']
        self.length           = length
    
    @classmethod
    def parse_raw_recovery(cls, mge_name: str, mge_result: list[MacSyHit], assembly_info: Assembly, 
                           EXPAND: int = 20):
        '''
        Parse the raw recovery result from MacSyFinder and construct a Recovery object.
        '''
        contig_locus = {}   # eg. {'contig1': [MacSyHit, MacSyHit], 'contig2': [MacSyHit, MacSyHit]}
        mpf, mob = set(), set()
        for hit in mge_result:
            contig_locus.setdefault(hit.loc_contig, []).append(hit)
            tags = get_feature(hit.gene_name)
            mpf.add(mge_name.split('_')[1])
            if 'Relaxase@' in tags:
                mob.add(tags.split('@')[1])

        middle_contig, edge_contigs, all_gene_num = {}, {}, set()
        for contig, locus_list in contig_locus.items():
            target = [gene.locus_num for gene in assembly_info.contigs[contig].genes.values()]
            mge_result_locus = [hit.locus_num for hit in mge_result if hit.loc_contig == contig]
            all_gene_num.update(mge_result_locus)
            if abs(max(target) - max(mge_result_locus)) < EXPAND and abs(min(target) - min(mge_result_locus)) < EXPAND:
                middle_contig[contig] = locus_list
            else:
                edge_contigs[contig] = locus_list
        
        # acturally, some times the edge contigs may also show character like the middle contigs, because
        # the MGE might be located near the edge of a contig, and still be considered part of the middle region, 
        # in that case, this algorithm cannot define them. However, that is very rare.
        if len(edge_contigs) != 2:
            return None
        
        left_contig, left_direction, left_gene, right_contig, right_direction, right_gene = direction_pending(
            edge_contigs, assembly_info, EXPAND)

        return (cls(middle_contig, left_contig, right_contig, left_direction, right_direction, 
                    left_gene, right_gene, mpf, mob, all_gene_num))

                   
    def merge_tRNA(self, assembly_info: Assembly, temp_dir: Path, EXPAND: int = 20):
        '''
        Search for tRNA genes in the flanking regions of the ICE and update the DR information accordingly.
        '''
        left_trna = search_trna(self.left_contig, self.left_gene, self.left_direction, assembly_info, 'left', 
                                EXPAND)
        if left_trna:
            self.trnalist.append(left_trna)
        right_trna = search_trna(self.right_contig, self.right_gene, self.right_direction, assembly_info, 'right', 
                                 EXPAND)
        if right_trna:
            self.trnalist.append(right_trna)

        self.expanded_left = left_trna if left_trna else self.left_gene
        self.expanded_right = right_trna if right_trna else self.right_gene

        self.DR1 = self.left_gene.start if self.left_direction == '+' else self.left_gene.end
        self.DR4 = self.right_gene.end if self.right_direction == '+' else self.right_gene.start

        for record in SeqIO.parse(assembly_info.sample_path, 'fasta'):
            if record.id == next(iter(self.left_contig.keys())):
                self.left_seq = record.seq if self.left_direction == '+' else record.seq.reverse_complement()
            elif record.id == next(iter(self.right_contig.keys())):
                self.right_seq = record.seq if self.right_direction == '+' else record.seq.reverse_complement()
            if self.left_seq and self.right_seq:
                break
        len_left_seq = len(self.left_seq)
        len_right_seq = len(self.right_seq)
        if self.trnalist:
            temp_DR = temp_dir / 'detect_DR.fasta'

            with open(temp_DR , 'w') as handle:
                handle.write('>temp_DR\n')
                handle.write(str(self.left_seq))
                handle.write(str(self.right_seq))
            
            rDRlist = get_DR(temp_DR, temp_dir)   # in case the direct of contig is wrong

            for DR in rDRlist:
                DRs = DR.split(',')
                if self.left_direction == '+' and self.right_direction == '+':
                    if ((self.expanded_left.start <= int(DRs[0]) <= self.expanded_left.end 
                         and int(DRs[3]) - len_left_seq - self.expanded_right.end <= 10000)
                        or (self.expanded_right.start <= int(DRs[3]) - len_left_seq <= self.expanded_right.end 
                            and self.expanded_left.start - int(DRs[0]) <= 10000)):
                        self.DR1 = int(DRs[0]) + 1
                        self.DR2 = int(DRs[1])
                        self.DR3 = int(DRs[2]) - len(self.left_seq) + 1
                        self.DR4 = int(DRs[3]) - len(self.left_seq)
                        self.DR_seq = str(getfa(temp_DR, 'temp_DR', self.DR1, self.DR2))
                        self.left_seq = self.left_seq[self.DR1:]
                        self.right_seq = self.right_seq[:self.DR4 - len_left_seq]
                elif self.left_direction == '-' and self.right_direction == '+':
                    if (((len_left_seq - self.expanded_left.end) <= int(DRs[0]) <= (len_left_seq - self.expanded_left.start) 
                         and int(DRs[3]) - len_left_seq - self.expanded_right.end <= 10000)
                        or (self.expanded_right.start <= int(DRs[3]) - len_left_seq <= self.expanded_right.end 
                            and len_left_seq - self.expanded_left.start - int(DRs[0]) <= 10000)):
                            self.DR1 = len_left_seq - int(DRs[1])
                            self.DR2 = len_left_seq - int(DRs[0]) + 1
                            self.DR3 = int(DRs[2]) - len(self.left_seq) + 1
                            self.DR4 = int(DRs[3]) - len(self.left_seq)
                            self.DR_seq = str(getfa(temp_DR, 'temp_DR', int(DRs[0]) + 1, int(DRs[1])))
                            self.left_seq = self.left_seq[int(DRs[0]) + 1:]
                            self.right_seq = self.right_seq[:self.DR4 - len_left_seq]
                else:   # self.left_direction == '+' and self.right_direction == '-', there are only 3 posibility because of pending_direction
                    if ((self.expanded_left.start <= int(DRs[0]) <= self.expanded_left.end 
                         and int(DRs[3]) - len_left_seq - len_right_seq + self.expanded_right.end <= 10000)
                        or self.expanded_right.start <= (len_right_seq - int(DRs[3]) + len_left_seq) <= self.expanded_right.end 
                        and self.expanded_left.start - int(DRs[0]) <= 10000):
                            self.DR1 = int(DRs[0]) + 1
                            self.DR2 = int(DRs[1])
                            self.DR3 = len_right_seq - int(DRs[3]) + len_left_seq + 1
                            self.DR4 = len_right_seq - int(DRs[2]) + len_left_seq
                            self.DR_seq = str(getfa(temp_DR, 'temp_DR', self.DR1, self.DR2))
                            self.left_seq = self.left_seq[self.DR1:]
                            self.right_seq = self.right_seq[:int(DRs[3]) - len_left_seq]

        else:
            if self.left_direction == '+':
                self.left_seq = self.left_seq[self.DR1:]
            else:
                self.left_seq = self.left_seq[:len(self.left_seq) - self.DR1]
            if self.right_direction == '+':
                self.right_seq = self.right_seq[:self.DR4]
            else:
                self.right_seq = self.right_seq[:len(self.right_seq) - self.DR4]

    def getrfasta(self, assembly_info: Assembly, rICE_count: int, output_dir: Path, temp_dir: Path, rootdir: Path, 
                  threads: int, verbose: bool = False):
        '''
        Generate the fasta file for the recovery ICE, and also fill in the attributes of the Recovery object.
        '''
        self.ICE_name = f'{assembly_info.sample_name}_rICE_{rICE_count}'
        
        final_dir = output_dir / assembly_info.sample_name
        if not final_dir.exists():
            final_dir.mkdir(parents = True)

        if self.left_direction == '+':
            left_ID = f'{assembly_info.sample_name}_{next(iter(self.left_contig.keys()))}_{int(self.DR1)}-{assembly_info.contigs[next(iter(self.left_contig.keys()))].length}'
            left_loc = f'{next(iter(self.left_contig.keys()))}:{int(self.DR1)}..{assembly_info.contigs[next(iter(self.left_contig.keys()))].length}'
        else:
            left_ID = f'{assembly_info.sample_name}_{next(iter(self.left_contig.keys()))}_1-{assembly_info.contigs[next(iter(self.left_contig.keys()))].length - self.DR1}'
            left_loc = f'{next(iter(self.left_contig.keys()))}:1..{assembly_info.contigs[next(iter(self.left_contig.keys()))].length - self.DR1}'
        
        if self.right_direction == '+':
            right_ID = f'{assembly_info.sample_name}_{next(iter(self.right_contig.keys()))}_1-{self.DR4}'
            right_loc = f'{next(iter(self.right_contig.keys()))}:1..{self.DR4}'
        else:
            right_ID = f'{assembly_info.sample_name}_{next(iter(self.right_contig.keys()))}_{assembly_info.contigs[next(iter(self.right_contig.keys()))].length - self.DR4}_{assembly_info.contigs[next(iter(self.right_contig.keys()))].length}'
            right_loc = f'{next(iter(self.right_contig.keys()))}:{assembly_info.contigs[next(iter(self.right_contig.keys()))].length - self.DR4}..{assembly_info.contigs[next(iter(self.right_contig.keys()))].length}'
        
        self.ICE_location = [left_loc]

        outfa = final_dir / f'{self.ICE_name}.fa'

        with open(outfa, 'w') as out_handle:
                SeqIO.write(SeqRecord(self.left_seq, id = left_ID, description=''), out_handle, 'fasta')

        all_seq = self.left_seq
        
        self.all_middle = [contig for contig, hit in self.middle_contig.items()]

        with open(outfa, 'a') as out_handle:
            for record in SeqIO.parse(assembly_info.sample_path, 'fasta'):
                if record.id in self.all_middle:
                    SeqIO.write(record, out_handle, 'fasta')
                    all_seq += record.seq

            SeqIO.write(SeqRecord(self.right_seq, id = right_ID, description=''), out_handle, 'fasta')
            all_seq += self.right_seq

        orifa = temp_dir / f'{self.ICE_name}_fororit.fasta'

        with open(orifa, 'w') as out_handle:
                SeqIO.write(SeqRecord(all_seq, id = 'fororit', description=''), out_handle, 'fasta')

        self.ICE_location += self.all_middle
        self.ICE_location.append(right_loc)
        self.gc = f'{gc_fraction(all_seq) * 100:.2f}'
        self.oritseqs = oritseq(orifa, self.ICE_name, 'fororit', 1, len(all_seq) + 1, rootdir, temp_dir, threads)
        if self.trnalist:
            trnaout = []
            for trna in self.trnalist:
                trnaout.append(f'{trna.product} ({trna.start}..{trna.end}) [{trna.strand}]')
            self.trnaout = (';').join(trnaout)
        self.length = len(all_seq)

        log(f'{self.ICE_name} recovery fasta generated.', verbose = verbose)

