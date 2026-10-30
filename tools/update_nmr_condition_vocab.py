#!/usr/bin/env python
"""Deterministically add experimental-condition tokens to BART dictionaries."""

from __future__ import annotations

import argparse
import runpy
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RULES = runpy.run_path(str(ROOT / "nmrpeak" / "utils" / "vocab_rules.py"))
NUMERICAL_RULES = RULES["NUMERICAL_RULES"]
CATEGORICAL_RULES = RULES["CATEGORICAL_RULES"]


MANAGED_PREFIXES = (
    "cnmr_condition_", "hnmr_condition_", "cnmr_solvent_", "hnmr_solvent_",
    "cnmr_temperature_", "hnmr_temperature_", "cnmr_field_", "hnmr_field_",
    "cnmr_ph_", "hnmr_ph_", "cnmr_concentration_", "hnmr_concentration_",
    "cnmr_line_width_", "hnmr_line_width_",
)


def generate_numerical_vocabulary(token_type: str) -> list[str]:
    rule = NUMERICAL_RULES[token_type]
    vocabulary = []
    for lower_bound, upper_bound, step, display_precision, template in rule["intervals"]:
        if step is None:
            vocabulary.append(template)
            continue
        count = int(round((upper_bound - lower_bound) / step))
        for index in range(count):
            lower = lower_bound + index * step
            upper = lower + step
            if display_precision == 0:
                vocabulary.append(template.format(lower=round(lower), upper=round(upper)))
            elif display_precision == 1:
                vocabulary.append(template.format(lower=round(lower, 1), upper=round(upper, 1)))
            else:
                vocabulary.append(template.format(lower=round(lower, 2), upper=round(upper, 2)))
    return vocabulary


def generated_tokens() -> list[str]:
    tokens: list[str] = []
    for prefix in ("cnmr", "hnmr"):
        tokens.extend((f"{prefix}_condition_start", f"{prefix}_condition_end"))
        tokens.append(f"{prefix}_solvent_missing")
        tokens.extend(f"{prefix}_solvent_{name}" for name in CATEGORICAL_RULES["nmr_solvent"])
        for family, missing_name in (
            ("temperature_k", "temperature_k"),
            ("field_strength_mhz", "field_strength_mhz"),
            ("ph", "ph"),
            ("concentration_log_molar", "concentration"),
        ):
            tokens.append(f"{prefix}_{missing_name}_missing")
            tokens.extend(generate_numerical_vocabulary(f"{prefix}_{family}"))
        tokens.append(f"{prefix}_line_width_hz_missing")
        tokens.extend(generate_numerical_vocabulary(f"{prefix}_line_width_hz"))
    return list(dict.fromkeys(tokens))


def update(check: bool = False) -> int:
    condition_path = ROOT / "dict" / "bart" / "condition_dict.txt"
    total_path = ROOT / "dict" / "bart" / "total_dict.txt"
    condition = [line for line in condition_path.read_text(encoding="utf-8").splitlines() if line]
    total = [line for line in total_path.read_text(encoding="utf-8").splitlines() if line]
    base_condition = [token for token in condition if not token.startswith(MANAGED_PREFIXES)]
    new_condition = base_condition + generated_tokens()
    try:
        smiles_start = total.index("#")
    except ValueError as exc:
        raise RuntimeError("Cannot locate the SMILES vocabulary suffix in total_dict.txt") from exc
    new_total = new_condition + total[smiles_start:]
    changed = condition != new_condition or total != new_total
    if check:
        if changed:
            print("Condition dictionaries are stale; run this tool without --check.")
            return 1
        print(f"Condition dictionaries are current ({len(new_condition)} condition tokens).")
        return 0
    condition_path.write_text("\n".join(new_condition) + "\n", encoding="utf-8")
    total_path.write_text("\n".join(new_total) + "\n", encoding="utf-8")
    print(f"Updated condition_dict.txt={len(new_condition)} and total_dict.txt={len(new_total)} tokens")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    raise SystemExit(update(parser.parse_args().check))
