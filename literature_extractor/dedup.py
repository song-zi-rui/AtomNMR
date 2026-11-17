"""Identity of experimental spectra and their reported atom assignments."""

from collections import defaultdict
from hashlib import sha256
import json
import re
import unicodedata


def fingerprint(value):
    return sha256(json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def normalized(value):
    if isinstance(value, (int, float)):
        return format(value, '.15g')
    return ' '.join(unicodedata.normalize('NFKC', str(value)).casefold().split()) if value is not None else ''


def source_bundle(checksums):
    has_sdf = any(name.lower().endswith('.sdf') for name in checksums)
    contents = checksums if has_sdf else sorted(set(checksums.values()))
    return fingerprint(contents)


def canonical_targets(row, mapping, overrides):
    tags = row.get('paper_atom_numbers')
    if tags is None:
        tags = []
        for label in row.get('source_labels', []):
            if label in overrides:
                tags.append(overrides[label])
                continue
            match = re.fullmatch(r'[HC](\d+)[a-z′’\']?', label)
            if not match:
                return None
            tags.append(int(match[1]))
    if any(tag not in mapping for tag in tags):
        return None
    return sorted({mapping[tag] for tag in tags})


def spectrum_keys(table):
    """Compare mapped spectra independently of page order and paper atom numbering."""
    if not table.get('smiles') or table.get('issues'):
        return {}
    conditions = table.get('conditions', {})
    mapping = {int(tag): index for tag, index in table['mapping'].items()}
    overrides = {f'H{label}': int(tag) for tag, label in table.get('atom_display', {}).items()}
    solvent = normalized(conditions.get('solvent'))
    sample = normalized(table.get('metadata', {}).get('sample'))
    # Without a recorded solvent, compare within the same article or PDF only.
    source = '' if solvent else table.get('doi') or table.get('sha256')
    keys = {}
    for nucleus, rows in table['peaks'].items():
        if not rows:
            continue
        targets_by_row = [canonical_targets(row, mapping, overrides) for row in rows]
        if any(targets is None for targets in targets_by_row):
            continue
        entries, groups = [], defaultdict(list)
        for row, targets in zip(rows, targets_by_row, strict=True):
            bounds = row.get('range_ppm')
            shifts = [normalized(shift) for shift in sorted(bounds)] if bounds else normalized(row.get('shift_ppm'))
            peak = [shifts, normalized(row.get('multiplicity')), normalized(row.get('reported_nH')),
                    sorted(normalized(value) for value in row.get('j_values_hz', [])),
                    normalized(row.get('j_values_reported'))]
            if row.get('group'):
                groups[row['group']].append((peak, targets))
            else:
                entries.append([peak, row.get('relation'), targets])
        entries.extend([sorted(fingerprint(peak) for peak, _ in group), 'interchangeable',
                        sorted({index for _, targets in group for index in targets})] for group in groups.values())
        context = [solvent, sample, normalized(conditions.get('frequencies_mhz', {}).get(nucleus)),
                   normalized(conditions.get('temperature_k')), normalized(conditions.get('chemical_shift_reference')),
                   source]
        keys[nucleus] = fingerprint([table['smiles'], nucleus, context, sorted(fingerprint(entry) for entry in entries)])
    return keys
