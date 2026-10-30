# rules.py
"""
化学位移和J值的token规则配置
"""

NUMERICAL_RULES = {
    # lower_bound, upper_bound, step, display_precision, template
    'cnmr_shift': {
        'intervals': [
            (-float('inf'), -6, None, None, "cnmr_shift_-∞_-6.0"),
            (-6, 0, 0.5, 1, "cnmr_shift_{lower:.1f}_{upper:.1f}"),
            (0, 200, 0.1, 1, "cnmr_shift_{lower:.1f}_{upper:.1f}"),
            (200, 210, 0.5, 1, "cnmr_shift_{lower:.1f}_{upper:.1f}"),
            (210, 220, 1, 0, "cnmr_shift_{lower:.1f}_{upper:.1f}"),
            (220, 250, None, None, "cnmr_shift_220.0_250.0"),
            (250, float('inf'), None, None, "cnmr_shift_250.0_∞")
        ],
        'precision': 1
    },
    'hnmr_shift': {
        'intervals': [
            (-float('inf'), -1, None, None, "hnmr_shift_-∞_-1.00"),
            (-1, 0, 0.1, 2, "hnmr_shift_{lower:.2f}_{upper:.2f}"),
            (0, 12, 0.05, 2, "hnmr_shift_{lower:.2f}_{upper:.2f}"),
            (12, 13, 0.1, 2, "hnmr_shift_{lower:.2f}_{upper:.2f}"),
            (13, 14, 0.2, 2, "hnmr_shift_{lower:.2f}_{upper:.2f}"),
            (14, 16, 0.5, 2, "hnmr_shift_{lower:.2f}_{upper:.2f}"),
            (16, float('inf'), None, None, "hnmr_shift_16.00_∞")
        ],
        'precision': 2
    },
    'hnmr_jvalue': {
        'intervals': [
            (0, 20, 0.5, 1, "hnmr_j_{lower:.1f}_{upper:.1f}"),
            (20, 30, 1, 0, "hnmr_j_{lower:.1f}_{upper:.1f}"),
            (30, 50, 2, 1, "hnmr_j_{lower:.1f}_{upper:.1f}"),
            (50, 60, 5, 1, "hnmr_j_{lower:.1f}_{upper:.1f}"),
            (60, 300, None, None, "hnmr_j_60.0_300.0"),
            (300, float('inf'), None, None, "hnmr_j_300.0_∞")
        ],
        'precision': 1
    }
}


# Experimental-condition tokens are nucleus-scoped because a molecule can have
# 1H and 13C spectra acquired under different conditions.  The compact bins are
# intentionally coarser than chemical-shift bins: these values are context, not
# prediction targets.  Concentration is represented in log10(mol/L).
_CONDITION_RULES = {
    'temperature_k': {
        'intervals': [
            (-float('inf'), 200, None, None, '{prefix}_temperature_k_-∞_200'),
            (200, 350, 5, 0, '{prefix}_temperature_k_{{lower:.0f}}_{{upper:.0f}}'),
            (350, float('inf'), None, None, '{prefix}_temperature_k_350_∞'),
        ],
        'precision': 0,
    },
    'field_strength_mhz': {
        'intervals': [
            (-float('inf'), 40, None, None, '{prefix}_field_mhz_-∞_40'),
            (40, 1200, 10, 0, '{prefix}_field_mhz_{{lower:.0f}}_{{upper:.0f}}'),
            (1200, float('inf'), None, None, '{prefix}_field_mhz_1200_∞'),
        ],
        'precision': 0,
    },
    'ph': {
        'intervals': [
            (-float('inf'), 0, None, None, '{prefix}_ph_-∞_0.0'),
            (0, 14.5, 0.5, 1, '{prefix}_ph_{{lower:.1f}}_{{upper:.1f}}'),
            (14.5, float('inf'), None, None, '{prefix}_ph_14.5_∞'),
        ],
        'precision': 1,
    },
    'concentration_log_molar': {
        'intervals': [
            (-float('inf'), -9, None, None, '{prefix}_concentration_logM_-∞_-9.0'),
            (-9, 2, 0.5, 1, '{prefix}_concentration_logM_{{lower:.1f}}_{{upper:.1f}}'),
            (2, float('inf'), None, None, '{prefix}_concentration_logM_2.0_∞'),
        ],
        'precision': 1,
    },
    'line_width_hz': {
        'intervals': [
            (-float('inf'), 0, None, None, '{prefix}_line_width_hz_-∞_0.0'),
            (0, 20, 0.5, 1, '{prefix}_line_width_hz_{{lower:.1f}}_{{upper:.1f}}'),
            (20, 100, 5, 0, '{prefix}_line_width_hz_{{lower:.0f}}_{{upper:.0f}}'),
            (100, 500, 20, 0, '{prefix}_line_width_hz_{{lower:.0f}}_{{upper:.0f}}'),
            (500, float('inf'), None, None, '{prefix}_line_width_hz_500_∞'),
        ],
        'precision': 1,
    },
}

for _prefix in ('cnmr', 'hnmr'):
    for _name, _rule in _CONDITION_RULES.items():
        NUMERICAL_RULES[f'{_prefix}_{_name}'] = {
            'precision': _rule['precision'],
            'intervals': [
                (lower, upper, step, display_precision, template.format(prefix=_prefix))
                for lower, upper, step, display_precision, template in _rule['intervals']
            ],
        }



CATEGORICAL_RULES = {
    'hnmr_nH': {
        'min': 1,
        'max': 50,
    },
    'hnmr_category' : ["AA'BB'",
    "AA'BB'C",
    'AB',
    'ABX',
    'ABq',
    'app_d',
    'app_dd',
    'app_dq',
    'app_dt',
    'app_q',
    'app_s',
    'app_t',
    'app_td',
    'br',
    'brd',
    'brdd',
    'brq',
    'brs',
    'brt',
    'd',
    'dd',
    'ddd',
    'dddd',
    'ddddd',
    'dddddd',
    'dddddt',
    'ddddq',
    'ddddt',
    'ddddtd',
    'dddp',
    'dddq',
    'dddqd',
    'dddt',
    'dddtd',
    'dddtt',
    'ddh',
    'ddp',
    'ddpd',
    'ddq',
    'ddqd',
    'ddqdd',
    'ddqt',
    'ddt',
    'ddtd',
    'ddtdd',
    'ddtdt',
    'ddtq',
    'ddtt',
    'ddttd',
    'dh',
    'dhd',
    'dhept',
    'dp',
    'dpd',
    'dpdd',
    'dpt',
    'dq',
    'dqd',
    'dqdd',
    'dqddd',
    'dqdt',
    'dqq',
    'dqt',
    'dqtd',
    'dt',
    'dtd',
    'dtdd',
    'dtddd',
    'dtddt',
    'dtdq',
    'dtdt',
    'dtdtd',
    'dtp',
    'dtq',
    'dtqd',
    'dtt',
    'dttd',
    'dttt',
    'h',
    'hd',
    'hdd',
    'hept',
    'heptd',
    'hex',
    'ht',
    'm',
    'p',
    'pd',
    'pdd',
    'pdt',
    'pq',
    'pt',
    'ptd',
    'q',
    'qd',
    'qdd',
    'qddd',
    'qddt',
    'qdq',
    'qdt',
    'qdtd',
    'qp',
    'qq',
    'qqd',
    'qt',
    'qtd',
    'qtdd',
    'qtt',
    's',
    'spt',
    't',
    'td',
    'tdd',
    'tddd',
    'tdddd',
    'tdddt',
    'tddq',
    'tddt',
    'tddtd',
    'tdp',
    'tdq',
    'tdqd',
    'tdt',
    'tdtd',
    'tdtdd',
    'tdtt',
    'th',
    'tp',
    'tpd',
    'tq',
    'tqd',
    'tqdd',
    'tqt',
    'tt',
    'ttd',
    'ttdd',
    'ttdt',
    'ttq',
    'ttt',
    'tttd'],
    'formula_element': {
        'min': 1,
        'max': 118,
    },
    'formula_count': {
        'min': 2,
        'max': 100,
    },

    'nmr_solvent': [
        'chloroform_d',
        'dmso_d6',
        'water_d2',
        'methanol_d4',
        'acetone_d6',
        'benzene_d6',
        'acetonitrile_d3',
        'pyridine_d5',
        'tetrahydrofuran_d8',
        'other',
    ],

}
