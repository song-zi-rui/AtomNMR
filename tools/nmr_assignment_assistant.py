#!/usr/bin/env python
"""Entry point for the independent 13C assignment-assistant pipeline."""

from pathlib import Path
import os
import sys
import types


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
os.environ["NMRPEAK_LIGHT_IMPORT"] = "1"

# Importing ``nmrpeak.shift_assignment`` would normally execute
# ``nmrpeak/__init__.py`` first.  Older checkouts eagerly import every Uni-Core
# task there, even though this standalone pipeline does not use Uni-Core.  A
# lightweight namespace package keeps build-db/report/predict usable inside
# Apptainer environments where Uni-Core is stored separately or not installed.
if "nmrpeak" not in sys.modules:
    package = types.ModuleType("nmrpeak")
    package.__path__ = [str(ROOT / "nmrpeak")]
    package.__package__ = "nmrpeak"
    sys.modules["nmrpeak"] = package

from nmrpeak.shift_assignment.cli import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
