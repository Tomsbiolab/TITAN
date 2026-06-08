#!/usr/bin/env python3
"""
aegis_merge.py
Merges ab initio and transcriptome annotations using the AEGIS library.
Called by the AEGIS_MERGE_ANNOTATIONS Nextflow process.
"""

import argparse
import os
from aegis.annotation import Annotation
from aegis.genome import Genome

parser = argparse.ArgumentParser()
parser.add_argument('--genome', required=True)
parser.add_argument('--ab_initio', nargs='+', required=True)
parser.add_argument('--transcriptome', nargs='+', required=True)
parser.add_argument('--version', required=True)
args = parser.parse_args()

genome_obj = Genome('RefGenome', args.genome)

def get_name(filepath):
    return os.path.splitext(os.path.basename(filepath))[0]

print('Processing Ab Initio evidences...')
ab_initio_annotations = [
    Annotation(name=get_name(f), annot_file_path=f, genome=genome_obj, rename_source=get_name(f))
    for f in args.ab_initio
]

print('Processing Transcriptome evidences...')
transcriptome_annotations = [
    Annotation(name=get_name(f), annot_file_path=f, genome=genome_obj, rename_source=get_name(f), rework_all_CDSs=True)
    for f in args.transcriptome
]

print('Merging Transcriptome evidences...')
transcript_evidence = transcriptome_annotations[0].copy()
for annotation in transcriptome_annotations[1:]:
    transcript_evidence.merge(annotation)

print('Merging all evidences...')
merged = ab_initio_annotations[0].copy()
for annotation in ab_initio_annotations[1:]:
    merged.merge(annotation)

merged.merge(transcript_evidence)
merged.make_alternative_transcripts_into_genes()

merged.id   = f'merge_{args.version}_on_genome'
merged.name = f'merge_{args.version}'

print('Exporting results...')
merged.export.gff('.', f'merge_{args.version}.gff3')
merged.export.unique_proteins(custom_path='.')
