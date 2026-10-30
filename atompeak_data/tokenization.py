"""Offline vocabularies and one-pass token/peak-position construction."""
from __future__ import annotations

import math
import re
from pathlib import Path

from .schema import BRANCHES, NUCLEI, finite, j_values
from .vocab_rules_partial import NUMERICAL_RULES, CATEGORICAL_RULES

ROOT = Path(__file__).resolve().parents[1]
RULES = {'NUMERICAL_RULES': NUMERICAL_RULES, 'CATEGORICAL_RULES': CATEGORICAL_RULES}


def numeric_token(value, family):
    value = finite(value)
    if value is None:
        return None
    rule = RULES['NUMERICAL_RULES'][family]
    value = round(value, rule['precision'])
    for low, high, step, precision, template in rule['intervals']:
        if low <= value < high:
            if step is None:
                return template
            if step < 1:
                multiplier = 10 ** rule['precision']
                start = math.floor(round(value * multiplier) / round(step * multiplier)) * round(step * multiplier) / multiplier
            else:
                start = math.floor(value) if step == 1 else low + math.floor(round((value-low)/step, 10)) * step
            return template.format(lower=round(start, precision), upper=round(start + step, precision))
    return None


class Vocabulary:
    def __init__(self, path):
        self.tokens = [line.split()[0] for line in Path(path).read_text(encoding='utf-8').splitlines() if line.strip()]
        # The legacy SMILES dictionary spells the backslash token as two
        # backslashes. It denotes ONE stereobond character, not two bonds.
        self.tokens = ['\\' if token == '\\\\' else token for token in self.tokens]
        self.indices = {value: i for i, value in enumerate(self.tokens)}
        if len(self.indices) != len(self.tokens) or self.tokens[:4] != ['[CLS]', '[PAD]', '[SEP]', '[UNK]']:
            raise ValueError('Vocabulary must have unique entries and [CLS]/[PAD]/[SEP]/[UNK] at 0/1/2/3')
        self.bos, self.pad, self.eos, self.unk = range(4)
        self.by_initial = {}
        for word in self.tokens[4:]:
            self.by_initial.setdefault(word[0], []).append(word)
        for values in self.by_initial.values():
            values.sort(key=lambda x: (-len(x), x))

    def __len__(self):
        return len(self.tokens)

    def token_id(self, token):
        return self.indices.get(token, self.unk)

    def encode_smiles(self, smiles):
        # Longest matching symbols from the bundled explicit SMILES vocabulary;
        # exact round-trip is checked, with no network tokenizer or target offset.
        ids, index = [self.bos], 0
        while index < len(smiles):
            word = next((w for w in self.by_initial.get(smiles[index], []) if smiles.startswith(w, index)), None)
            if word is None:
                raise ValueError(f'SMILES character not represented in vocabulary: {smiles[index:index+12]!r}')
            ids.append(self.indices[word]); index += len(word)
        ids.append(self.eos)
        if self.decode_smiles(ids) != smiles:
            raise ValueError('SMILES tokenizer failed exact round-trip')
        return ids

    def decode_smiles(self, ids):
        result = []
        for token in ids:
            token = int(token)
            if token == self.eos:
                break
            if token in (self.bos, self.pad):
                continue
            result.append(self.tokens[token])
        return ''.join(result)


def condition_tokens(record, branch, use_solvent, use_temperature, use_field_strength):
    if not any((use_solvent, use_temperature, use_field_strength)):
        return []
    prefix = 'cnmr' if branch == 'carbon' else 'hnmr'
    values = record.get('conditions', {}).get('by_nucleus', {}).get(NUCLEI[branch], {})
    tokens = [prefix + '_condition_start']
    if use_solvent:
        solvent = str(values.get('solvent') or '').strip().lower()
        aliases = {'cdcl3':'chloroform_d','dmso-d6':'dmso_d6','d2o':'water_d2','cd3od':'methanol_d4',
                   'acetone-d6':'acetone_d6','c6d6':'benzene_d6','cd3cn':'acetonitrile_d3',
                   'pyridine-d5':'pyridine_d5','thf-d8':'tetrahydrofuran_d8'}
        solvent = aliases.get(solvent, solvent)
        solvent = solvent if solvent in RULES['CATEGORICAL_RULES']['nmr_solvent'] else 'other' if solvent else 'missing'
        tokens.append(f'{prefix}_solvent_{solvent}')
    for key, enabled in (('temperature_k', use_temperature), ('field_strength_mhz', use_field_strength)):
        if enabled:
            tokens.append(numeric_token(values.get(key), prefix + '_' + key) or prefix + '_' + key + '_missing')
    tokens.append(prefix + '_condition_end')
    return tokens


def tokenize_spectrum(record, vocabulary, *, use_formula=True, use_h_category=True,
                      use_h_nh=True, use_h_jvalue=True, use_range_shift=False,
                      use_solvent=False, use_temperature=False, use_field_strength=False,
                      max_length=1024):
    tokens, positions = [], {b: [] for b in BRANCHES}
    for branch in BRANCHES:
        peaks = record['spectra'][branch]['peaks']
        if not peaks:
            continue
        conditions = condition_tokens(record, branch, use_solvent, use_temperature, use_field_strength)
        missing_conditions = [token for token in conditions if token not in vocabulary.indices]
        if missing_conditions:
            raise ValueError('Dataset vocabulary lacks condition tokens: ' + ', '.join(missing_conditions)
                + '. Use the condition_dict.txt exported with this dataset; do not silently replace conditions with UNK.')
        tokens.extend(conditions)
        for peak in peaks:
            positions[branch].append(len(tokens) + 1)  # Leading BOS is added below.
            if branch == 'carbon':
                tokens.append(numeric_token(peak['shift_ppm'], 'cnmr_shift'))
                continue
            shifts = [peak['rangeMax'], peak['rangeMin']] if use_range_shift else [peak['shift_ppm']]
            tokens.extend(numeric_token(value, 'hnmr_shift') for value in shifts)
            if use_h_category:
                category = str(peak.get('category') or 'unk').replace(' ', '_')
                tokens.append('hnmr_category_' + category if category in RULES['CATEGORICAL_RULES']['hnmr_category'] else 'hnmr_unk')
            if use_h_nh:
                nh = finite(peak.get('nH'))
                tokens.append(f'hnmr_nH_{int(nh)}' if nh is not None and nh.is_integer() and 1 <= nh <= 50 else 'hnmr_unk')
            if use_h_jvalue:
                tokens.extend(numeric_token(j, 'hnmr_jvalue') for j in j_values(peak.get('j_values')))
        tokens.append('cnmr_end' if branch == 'carbon' else 'hnmr_end')
    if use_formula:
        for element, count in re.findall(r'([A-Z][a-z]?|[+-])(\d*)', record['molecular_formula']):
            tokens.append('formula_charge_' + element if element in '+-' else 'formula_ele_' + element)
            if count and int(count) > 1:
                tokens.append('formula_cnt_' + count)
    if len(tokens) + 2 > max_length:
        raise ValueError(f'Spectrum sequence length {len(tokens)+2} exceeds {max_length}; peaks are never silently truncated')
    ids = [vocabulary.bos] + [vocabulary.token_id(t) for t in tokens] + [vocabulary.eos]
    return ids, positions
