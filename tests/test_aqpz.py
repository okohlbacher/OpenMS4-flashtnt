#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Oliver Kohlbacher $
"""Exercise AQPZ identification and report comparison with the older app cache."""
import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import subprocess
import time


def check_results(directory, reference):
    tables = {}
    for name, columns in reference['columns'].items():
        with (directory / name).open(encoding='utf-8', newline='') as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            # Pinned upstream renamed the PrSM row index after the 2025 cache.
            if name == 'prsms.tsv':
                columns = ['PrSMIndex' if c == 'ProteoformIndex' else c for c in columns]
            if not set(columns).issubset(reader.fieldnames or []):
                raise ValueError(f'{name}: missing expected output columns')
            tables[name] = list(reader)
        if not tables[name]:
            raise ValueError(f'{name}: no identifications produced')

    best = max(tables['protein.tsv'], key=lambda row: float(row['Score']))
    expected = reference['best_aqpz']
    for field in ('ProteinAccession', 'DatabaseSequence', 'ProteinSequence'):
        if best[field] != expected[field]:
            raise ValueError(f'Highest-scoring identification has unexpected {field}')
    for field in ('Score', 'ProteoformMass', 'MatchingFragments', 'Coverage(%)'):
        value = float(best[field])
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f'AQPZ has invalid {field}: {value}')
    precursor_mass = float(best['PrecursorMass'])
    # The retained deconvolution has an unknown precursor, represented by -1.
    if not math.isfinite(precursor_mass) or not (precursor_mass == -1 or precursor_mass > 0):
        raise ValueError(f'AQPZ has invalid PrecursorMass: {precursor_mass}')
    if not (1 <= int(best['StartPosition']) <= int(best['EndPosition']) <= len(best['DatabaseSequence'])):
        raise ValueError('AQPZ identification positions are outside its database sequence')
    if float(best['Coverage(%)']) > 100:
        raise ValueError('AQPZ coverage exceeds 100 percent')
    aqpz_tags = [row for row in tables['tags.tsv'] if 'AQPZ' in row['ProteinAccession'].split(';')]
    if not aqpz_tags or any(int(row['Length']) < 4 for row in aqpz_tags):
        raise ValueError('Expected AQPZ sequence tags of at least four residues')

    comparison = {}
    for name, rows in tables.items():
        count = reference['row_counts'][name]
        comparison[name + '_rows'] = {'actual': len(rows), 'reference': count, 'match': len(rows) == count}
    actual_tags = Counter(row['TagSequence'] for row in tables['tags.tsv'])
    old_tags = Counter(reference['tag_sequence_counts'])
    comparison['tag_sequence_counts'] = {'match': actual_tags == old_tags,
        'extra_count': sum((actual_tags - old_tags).values()),
        'missing_count': sum((old_tags - actual_tags).values())}
    for field in ('Score', 'PrecursorMass', 'ProteoformMass', 'MatchingFragments',
                  'Coverage(%)', 'StartPosition', 'EndPosition'):
        actual, old = float(best[field]), float(expected[field])
        comparison[field] = {'actual': actual, 'reference': old,
                            'match': math.isclose(actual, old, rel_tol=1e-6, abs_tol=1e-6)}
    return {'positive_identity_checks_passed': True, 'aqpz_tags': len(aqpz_tags),
            'legacy_reference_matches': all(item['match'] for item in comparison.values()),
            'legacy_reference_comparison': comparison}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--strict-reference', action='store_true',
                        help='Also require numeric parity with the older, unversioned binary output')
    args = parser.parse_args()
    fixture = Path(__file__).resolve().parent / 'data' / 'aqpz'
    reference = json.loads((fixture / 'reference.json').read_text())
    parameters = json.loads((fixture / 'FTnT_parameters.json').read_text())
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    command = [str(args.executable.resolve()), '-in', str(fixture / 'out_deconv.mzML'),
               '-fasta', str(fixture / 'database.fasta'), '-threads', '2']
    for key, filename in (('out_tag', 'tags.tsv'), ('out_pro', 'protein.tsv'), ('out_prsm', 'prsms.tsv')):
        command.extend(['-' + key, str(output / filename)])
    for key, value in parameters.items():
        command.extend(['-' + key, *str(value).splitlines()])
    started = time.perf_counter()
    with (output / 'run.log').open('w', encoding='utf-8') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=180)
    report = check_results(output, reference)
    report.update({'wall_seconds': time.perf_counter() - started, 'command': command})
    (output / 'acceptance.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    if args.strict_reference and not report['legacy_reference_matches']:
        raise ValueError('Historical numeric fields differ; see acceptance.json')


if __name__ == '__main__':
    main()
