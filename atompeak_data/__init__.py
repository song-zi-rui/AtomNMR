"""Portable AtomPeak data contract, independent of PyTorch and Uni-Core."""

from .schema import SCHEMA_VERSION, build_record, split_for_connectivity

__all__ = ['SCHEMA_VERSION', 'build_record', 'split_for_connectivity']
