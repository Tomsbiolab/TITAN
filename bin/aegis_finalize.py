#!/usr/bin/env python3
"""
aegis_finalize.py
Filters redundant annotations and exports the final annotation set using AEGIS.
Called by the AEGIS_FILTER_REDUNDANCY Nextflow process.
"""

import argparse
import os
from aegis.utils.misc import pickle_save
from aegis.annotation import Annotation
from aegis.genome import Genome

parser = argparse.ArgumentParser()
parser.add_argument('--genome', required=True)
parser.add_argument('--masked_genome', required=True)
parser.add_argument('--merged_gff', required=True)
parser.add_argument('--blast_results', nargs='+', required=True)
parser.add_argument('--priority', required=True)
parser.add_argument('--version', required=True)
args = parser.parse_args()

genome_obj        = Genome('RefGenome', args.genome)
masked_genome_obj = Genome('MaskedGenome', args.masked_genome)

merge = Annotation(
    name='merge',
    annot_file_path=args.merged_gff,
    genome=genome_obj,
    hard_masked_genome=masked_genome_obj
)

merge.stats.calculate_transcript_masking()
merge.update()
merge.overlaps.detect()

print('Adding BLAST results...')
for blast_file in args.blast_results:
    source_name = blast_file.replace('_merged.diamond', '')
    print(f"Loading BLAST hits for {source_name}...")
    merge.add_blast_hits(source_name, blast_file)

print('Saving intermediate overlaps...')
pickle_save(f'merge_overlaps_{args.version}.pkl', merge)
merge.update()
merge.export.gff('.', f'merge_overlaps_{args.version}.gff3')

print('Reducing redundancy...')
priority_sources = args.priority.split(',')
merge.redundancy.filter(source_priority=priority_sources)

merge.id   = f'final_annotation_{args.version}_on_genome'
merge.name = f'final_annotation_{args.version}'

print('Exporting final outputs...')
merge.export.gff('.', f'final_annotation_{args.version}.gff3')
merge.overlaps.export(overlap_threshold=0, export_csv=True, export_self=True, NAs=False, custom_path='.')
merge.export.proteins(only_main=True,  verbose=False, custom_path='.')
merge.export.proteins(only_main=False, verbose=False, custom_path='.')
