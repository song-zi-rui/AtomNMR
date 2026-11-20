"""Stored computational spectra, addressed through the same molecule IDs."""

import json

IBM = 'IBM_DFT'
DFT8K = 'CASCADE-DFT8K'
MAGNET = 'MagNET-DFT'
CASPRE = 'CASPRE-JEOL'


def initialize(c):
    c.execute("""CREATE TABLE IF NOT EXISTS computed_spectra(
        id TEXT PRIMARY KEY,
        molecule_id TEXT NOT NULL REFERENCES molecules(id),
        nucleus TEXT NOT NULL CHECK(nucleus IN ('1H','13C')),
        source TEXT NOT NULL,
        source_record_id TEXT NOT NULL,
        peaks_json TEXT NOT NULL,
        source_json TEXT NOT NULL,
        UNIQUE(source,source_record_id,nucleus)
    )""")
    c.execute("CREATE INDEX IF NOT EXISTS computed_spectra_molecule ON computed_spectra(molecule_id,nucleus,id)")


def exists(c):
    return c.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='computed_spectra'").fetchone() is not None


def for_molecule(c, mid):
    if not exists(c):
        return []
    result = []
    for row in c.execute('SELECT * FROM computed_spectra WHERE molecule_id=? ORDER BY nucleus,source,id', (mid,)):
        spectrum = dict(row)
        spectrum['peaks'] = json.loads(spectrum.pop('peaks_json'))
        spectrum['provenance'] = json.loads(spectrum.pop('source_json'))
        spectrum.update(
            spectrum_kind='computed', editable=False, manual_state='active',
            conditions={}, solvent=None, frequency_mhz=None, temperature_k=None,
            quality={'decision': 'computed_reference'},
        )
        result.append(spectrum)
    return result


def pairs(spectra):
    """Pair H/C only when the same computation supplied both nuclei."""
    by_record = {}
    for s in spectra:
        if s.get('spectrum_kind') == 'computed':
            by_record.setdefault((s['source'], s['source_record_id']), {})[s['nucleus']] = s['id']
    return [
        {'id': f'computed:{source}:{record}', 'c_id': pair['13C'], 'h_id': pair['1H'],
         'solvent_match': 'unknown', 'basis': 'same_computed_source_record'}
        for (source, record), pair in sorted(by_record.items()) if '13C' in pair and '1H' in pair
    ]


def statistics(c):
    if not exists(c):
        return {'molecules': 0, 'spectra': 0}
    molecules, spectra = c.execute('SELECT count(DISTINCT molecule_id),count(*) FROM computed_spectra').fetchone()
    return {'molecules': molecules, 'spectra': spectra}
