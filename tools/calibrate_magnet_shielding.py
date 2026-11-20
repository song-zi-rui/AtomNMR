"""Fit each MagNET DFT method to assigned Delta22 chloroform sites."""

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from scan_magnet_structures import OUT

REVISION = 'fe222d1894afb2d954fe8d2e67e70731a006c8d2'
DELTA22_SHA256 = '5238ded71fbf9ed08b6b1cbd3813a81e54ee20fdbc480b7bceacdbf1f556e915'
EXP_GIT = '95987533ce02b7d7b47b7a55e5736d490195e249'
MODELS = {
    'pbe0_gas': ('pbe0_d3bj,pcSseg1,gas,"none"', 'pbe0_d3bj,pcSseg1,gas,"none"'),
    'pepper_gas': ('wp04,pcSseg2,gas,"none"', 'wb97xd,pcSseg2,gas,"none"'),
    'b3lyp_pcm': ('b3lyp_d3bj,pcSseg2,pcm,"chloroform"', 'b3lyp_d3bj,pcSseg2,pcm,"chloroform"'),
}


def fit(x, y):
    intercept, slope = np.linalg.lstsq(np.column_stack((np.ones(len(x)), x)), y, rcond=None)[0]
    return float(intercept), float(slope)


def model_for(ids, x, y):
    full = fit(x, y)
    scores = np.empty(len(y))
    for name in sorted(set(ids)):
        test = ids == name
        intercept, slope = fit(x[~test], y[~test])
        scores[test] = intercept + slope * x[test] - y[test]
    return {'intercept': full[0], 'slope': full[1], 'sites': len(y), 'molecules': len(set(ids)),
            'loso_rmse_ppm': float(np.sqrt(np.mean(scores**2))),
            'loso_mae_ppm': float(np.mean(np.abs(scores))),
            'shielding_min_max_ppm': [float(x.min()), float(x.max())]}


def calibrate(delta22, experimental):
    with delta22.open('rb') as f:
        if hashlib.file_digest(f, 'sha256').hexdigest() != DELTA22_SHA256:
            raise ValueError('Calibration HDF5 checksum mismatch')
    raw = experimental.read_bytes()
    if hashlib.sha1(f'blob {len(raw)}\0'.encode() + raw).hexdigest() != EXP_GIT:
        raise ValueError('Calibration spreadsheet checksum mismatch')
    sites = pd.read_excel(experimental).dropna(subset=['chloroform'])
    models, anchor_rows = {}, []
    with h5py.File(delta22) as f:
        names = [s.decode() for s in f['conventional_nmr_method_names'][:]]
        shieldings = {s: f['solutes'][s]['stationary_and_pcm']['conventional_shieldings'][:, 0, :]
                      for s in set(sites['solute'])}
        for model, methods in MODELS.items():
            models[model] = {}
            for nucleus, method in zip(('1H', '13C'), methods):
                atom_type = {'1H': 'H', '13C': 'C'}[nucleus]
                rows = []
                for row in sites.itertuples(index=False):
                    if not row.site.endswith('_' + atom_type):
                        continue
                    indices = [int(i) - 1 for i in str(row.atom_numbers).split(',')]
                    raw_sigma = shieldings[row.solute][names.index(method), indices]
                    if np.any(raw_sigma == np.iinfo(np.int32).min):
                        continue
                    rows.append((row.solute, row.site, float(np.mean(raw_sigma) * 1e-4), float(row.chloroform)))
                ids = np.array([r[0] for r in rows])
                xy = np.array([r[2:] for r in rows])
                models[model][nucleus] = {**model_for(ids, *xy.T), 'calibration_dft_method': method}
                anchor_rows.extend({'model': model, 'nucleus': nucleus, 'solute': a, 'site': b,
                                    'shielding_ppm': x, 'experimental_chloroform_ppm': y} for a, b, x, y in rows)
    return {
        'id': 'magnet-delta22-chloroform-v1',
        'source_url': f'https://huggingface.co/ekwan16/MagNET/tree/{REVISION}/data/delta22',
        'delta22_sha256': DELTA22_SHA256, 'experimental_git_blob': EXP_GIT, 'models': models,
        'geometry_note': 'Delta22 AIMNet2 geometry. Transfer to sigma-shake B3LYP geometry is not validated here.',
        'validation': 'leave-one-solute-out, all assigned Delta22 sites with measured chloroform shift',
        'measurement_solvent': 'chloroform',
        'anchors': anchor_rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('delta22', type=Path)
    parser.add_argument('experimental', type=Path)
    parser.add_argument('--output', type=Path, default=OUT / 'shielding_calibration.json')
    args = parser.parse_args()
    report = calibrate(args.delta22, args.experimental)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report['models'], indent=2), flush=True)


if __name__ == '__main__':
    main()
