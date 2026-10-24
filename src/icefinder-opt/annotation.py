# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 15:46:48 2026

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

from __future__ import annotations
from typing import TYPE_CHECKING

import csv
from copy import deepcopy
from dataclasses import dataclass, field, replace
from io import StringIO
from pathlib import Path
import re

from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

from .log import ICEfinderError, log
from .commands import run_command

if TYPE_CHECKING:
    from .assembly import Contig


# ─────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────
@dataclass(frozen = True)
class ResistanceHit:
    '''
    Store resistance-search evidence and its location on the original contig.
    '''
    category: str
    contig: str
    location: SimpleLocation | CompoundLocation
    gene: str
    product: str
    accession: str = ''
    identity: float = 0
    coverage: float = 0
    score: float = 0
    query: str = ''
    database_gene: str = ''

    @property
    def quality(self) -> tuple[float, float, float, float, str, str]:
        '''
        Return a deterministic ranking of resistance hits by score, coverage, identity, and name.
        '''
        return (self.score, self.coverage * self.identity, self.coverage,
                self.identity, self.gene, self.accession)


@dataclass
class ResistanceAnnotations:
    '''
    Return complete resistance-search evidence and gene mappings without changing an Assembly object.
    '''
    hits: list[ResistanceHit] = field(default_factory = list)
    amr_genes: dict[str, str] = field(default_factory = dict)
    metal_genes: dict[str, str] = field(default_factory = dict)


# ─────────────────────────────────────────────
# Core methods
# ─────────────────────────────────────────────
def overlap(a: SimpleLocation | CompoundLocation, b: SimpleLocation | CompoundLocation) -> int:
    '''
    Count shared bases between simple or joined feature locations.
    '''
    return sum(max(0, min(int(x.end), int(y.end)) - max(int(x.start), int(y.start)))
               for x in a.parts for y in b.parts)


def matching_cds(record: SeqRecord, hit: ResistanceHit) -> SeqFeature | None:
    '''
    Find the same-strand CDS covering at least 80 percent of a resistance hit.
    '''
    candidates = []
    for index, feature in enumerate(record.features):
        if feature.type != 'CDS' or feature.location is None:
            continue
        if feature.location.strand != hit.location.strand:
            continue
        covered = overlap(feature.location, hit.location)
        # Require most of the hit to lie within a CDS on the same strand.
        if covered >= 0.8 * len(hit.location):
            candidates.append((covered, -index, feature))
    return max(candidates, key = lambda item: item[:2])[2] if candidates else None


def parse_resfinder(output: str) -> list[ResistanceHit]:
    '''
    Parse ABRicate results into resistance hits with zero-based, half-open source coordinates.
    '''
    hits = []
    if not output.strip():
        return hits
    for row in csv.DictReader(StringIO(output), delimiter = '\t'):
        if not row.get('SEQUENCE'):
            continue
        start, end = sorted((int(row['START']), int(row['END'])))
        # Convert ABRicate's inclusive coordinates to Biopython's half-open coordinates.
        # Remove only the numeric database-entry suffix; retain the full identifier in notes.
        hits.append(ResistanceHit(
            'AMR', row['SEQUENCE'], SimpleLocation(start - 1, end, strand = -1 if row['STRAND'] == '-' else 1),
            re.sub(r'_\d+$', '', row['GENE']), row.get('PRODUCT') or row['GENE'], row.get('ACCESSION', ''),
            float(row['%IDENTITY']), float(row['%COVERAGE']), database_gene = row['GENE']))
    return hits


def parse_metal(output: str, query_features: dict[str, tuple[str, SeqFeature]]) -> list[ResistanceHit]:
    '''
    Parse metal BLAST results and retain the highest-quality qualifying hit for each protein.
    '''
    best = {}
    for line in output.splitlines():
        if not line.strip():
            continue
        fields = line.split('\t', 13)
        if len(fields) < 14:
            raise ICEfinderError('Metal BLAST output is missing subject descriptions.')
        query, subject = fields[:2]
        identity, coverage = float(fields[2]), 100 * int(fields[3]) / int(fields[12])
        # Both values are percentages; divide by 10000 to retain the original 0.64 cutoff.
        if identity * coverage / 10000 < 0.64 or query not in query_features:
            continue
        contig, feature = query_features[query]
        title = fields[13]
        match = re.search(r'\bGN=([^\s]+)', title)
        # The bundled database encodes the gene symbol after its BAC accession.
        parts = subject.removeprefix('lcl|').split('|')
        accession = parts[0]
        gene = parts[1] if len(parts) > 1 else (match.group(1) if match else accession)
        product = re.split(r'\s+OS=', title, maxsplit = 1)[0]
        if product.startswith(subject + ' '):
            product = product[len(subject) + 1:]
        hit = ResistanceHit('metal', contig, feature.location, gene, product or gene,
                            accession, identity, coverage, float(fields[11]), query)
        if query not in best or hit.quality > best[query].quality:
            best[query] = hit
    return list(best.values())


def query_feature_map(contigs: dict[str, Contig], records: dict[str, SeqRecord]) -> dict[str, tuple[str, SeqFeature]]:
    '''
    Link internal protein IDs to actual input CDS locations, never Gene.start/end.
    '''
    result = {}
    for contig, info in contigs.items():
        by_locus = {}
        # Original CDS locations also retain joined features and their strand information.
        for feature in records[contig].features:
            if feature.type == 'CDS':
                for locus in feature.qualifiers.get('locus_tag', []):
                    by_locus[locus] = feature
        for query, gene in info.genes.items():
            if gene.origin_id in by_locus:
                result[query] = (contig, by_locus[gene.origin_id])
    return result


def enrich_records(records: dict[str, SeqRecord], hits: list[ResistanceHit]) -> dict[str, SeqRecord]:
    '''
    Update resistance gene and product qualifiers while preserving the original annotations in notes.
    '''
    enriched = deepcopy(records)
    groups = {}
    for hit in hits:
        record = enriched[hit.contig]
        feature = matching_cds(record, hit)
        if feature is None:
            # A local alignment supports a feature annotation, not a complete predicted CDS.
            feature = next((f for f in record.features
                            if f.type == 'misc_feature' and f.location == hit.location
                            and f.qualifiers.get('inference') == ['similar to sequence:ICEfinder-opt_v1']), None)
        if feature is None:
            feature = SeqFeature(hit.location, type = 'misc_feature', qualifiers = {
                'inference': ['similar to sequence:ICEfinder-opt_v1']})
            record.features.append(feature)
        key = (hit.contig, id(feature))
        groups.setdefault(key, (feature, []))[1].append(hit)

    for feature, evidence in groups.values():
        notes = feature.qualifiers.setdefault('note', [])
        # Retain the original names and descriptions before replacing either qualifier.
        for key in ('gene', 'product'):
            for original in feature.qualifiers.get(key, []):
                note = f'Original {key}: {original}'
                if note not in notes:
                    notes.append(note)
        selected = []
        # Keep the best description from each category when one CDS has both functions.
        for category in ('AMR', 'metal'):
            candidates = [h for h in evidence if h.category == category]
            if candidates:
                selected.append(max(candidates, key = lambda h: h.quality))
        feature.qualifiers['gene'] = list(dict.fromkeys(h.gene for h in selected))
        feature.qualifiers['product'] = ['; '.join(
            f'{h.product} [{"antimicrobial resistance gene" if h.category == "AMR" else "metal resistance gene"}]'
            for h in selected)]
        for hit in sorted(evidence, key = lambda h: (h.category, h.gene, str(h.location))):
            note = (f'ICEfinder-opt_v1 {hit.category}: {hit.gene}; accession={hit.accession}; '
                    f'database_gene={hit.database_gene or hit.gene}; '
                    f'identity={hit.identity:.2f}%; coverage={hit.coverage:.2f}%; '
                    f'source_location={hit.contig}:{hit.location}')
            if note not in notes:
                notes.append(note)
    return enriched


def scan_resistance(sample_path: Path, faa_file: Path, query_features: dict[str, tuple[str, SeqFeature]],
                    rootdir: Path, threads: int = 1, verbose: bool = False) -> ResistanceAnnotations:
    '''
    Search AMR and metal databases and return a complete result for the Assembly to commit.
    '''
    log('Running cached resfinder search on the complete assembly (once per sample).', verbose)
    amr = parse_resfinder(run_command(['abricate', '-db', 'resfinder', sample_path]).stdout)
    log('Running cached metal protein search (once per sample).', verbose)
    metal = parse_metal(run_command([
        'blastp', '-query', str(faa_file), '-db', str(rootdir / 'data' / 'metal'),
        '-evalue', '1e-4', '-num_threads', str(threads), '-max_hsps', '1',
        '-max_target_seqs', '1', '-outfmt', '6 std slen stitle']).stdout, query_features)
    result = ResistanceAnnotations(hits = amr + metal)
    for query, (contig, feature) in query_features.items():
        for category, target in [('AMR', result.amr_genes), ('metal', result.metal_genes)]:
            matching = [h for h in result.hits if h.category == category and h.contig == contig
                        and h.location.strand == feature.location.strand
                        and overlap(h.location, feature.location) >= 0.8 * len(h.location)]
            if matching:
                target[query] = ';'.join(sorted({h.gene for h in matching}))
    return result


def resistance_names(hits: list[ResistanceHit], regions: list[tuple[str, SimpleLocation]]) -> tuple[str, str]:
    '''
    Return unique AMR and metal gene names for hits fully contained in the source regions.
    '''
    names = {'AMR': set(), 'metal': set()}
    for hit in hits:
        if any(hit.contig == contig and overlap(hit.location, location) == len(hit.location)
               for contig, location in regions):
            names[hit.category].add(hit.gene)
    return tuple(';'.join(sorted(names[key])) or '-' for key in ('AMR', 'metal'))


def add_site(record: SeqRecord, label: str, location: SimpleLocation) -> None:
    '''
    Add an att or oriT feature unless an equivalent annotation already exists.
    '''
    for feature in record.features:
        text = ' '.join(str(v) for key in ('label', 'note', 'gene', 'product')
                        for v in feature.qualifiers.get(key, []))
        if (feature.location is not None and re.search(rf'\b{re.escape(label)}\b', text, re.I)
                and [(int(p.start), int(p.end)) for p in feature.location.parts]
                == [(int(p.start), int(p.end)) for p in location.parts]
                and feature.location.strand in (None, location.strand)):
            return
    record.features.append(SeqFeature(location,
        type = 'repeat_region' if label in ('attL', 'attR') else 'misc_feature',
        qualifiers = {'label': [label], 'note': [f'{label} predicted by ICEfinder-opt_v1']}))


def source_region(record: SeqRecord, region: slice, reverse: bool = False) -> SimpleLocation:
    '''
    Convert an oriented fragment slice to its location on the original contig.
    '''
    start, stop, _ = region.indices(len(record))
    if reverse:
        start, stop = len(record) - stop, len(record) - start
    return SimpleLocation(start, stop, strand = -1 if reverse else 1)


def enrich_fragment(record: SeqRecord, contig: str, region: SimpleLocation,
                    hits: list[ResistanceHit]) -> SeqRecord:
    '''
    Map contained resistance hits to an output fragment and update its GenBank annotations.
    '''
    local_hits = []
    for hit in hits:
        if hit.contig != contig or overlap(hit.location, region) != len(hit.location):
            continue
        parts = []
        # Project every part independently so joined CDS locations remain valid.
        for part in hit.location.parts:
            start, end = int(part.start) - int(region.start), int(part.end) - int(region.start)
            strand = part.strand
            if region.strand == -1:
                # Reverse both the coordinates and the strand on reverse-complemented fragments.
                start, end = len(region) - end, len(region) - start
                strand = -strand if strand else strand
            parts.append(SimpleLocation(start, end, strand = strand))
        location = (parts[0] if len(parts) == 1 else
                    CompoundLocation(parts, operator = hit.location.operator))
        local_hits.append(replace(hit, contig = record.id, location = location))
    return enrich_records({record.id: record}, local_hits)[record.id]


def add_att(record: SeqRecord, label: str, start: int, end: int, region: SimpleLocation) -> None:
    '''
    Map a 1-based source repeat onto an oriented output fragment.
    '''
    if not start or not end:
        return
    start, end = min(start, end) - 1, max(start, end)
    if int(region.start) <= start < end <= int(region.end):
        if region.strand == -1:
            start, end = int(region.end) - end, int(region.end) - start
        else:
            start, end = start - int(region.start), end - int(region.start)
        add_site(record, label, SimpleLocation(start, end, strand = region.strand))


def add_orit(records: list[SeqRecord], hits: list[dict[str, int]]) -> None:
    '''
    Annotate oriT hits within individual fragments, excluding matches across artificial joins.
    '''
    offset = 0
    for record in records:
        for hit in hits:
            # A valid hit must fit within one real fragment of the concatenated query.
            if offset <= hit['start'] < hit['end'] <= offset + len(record):
                add_site(record, 'oriT', SimpleLocation(hit['start'] - offset, hit['end'] - offset,
                                                       strand = hit['strand']))
        offset += len(record)
