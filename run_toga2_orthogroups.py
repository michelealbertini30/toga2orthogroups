#!/usr/bin/env python3
"""
toga2orthogroups — standalone entry point.

Delegates to src.toga2orthogroups, which contains the full implementation.
Run from the repo root directory so that the src/ package is on the import path.

Usage:
    python toga2orthogroups.py -t DIR -s FILE -b FILE -i FILE -o DIR [options]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.toga2orthogroups import main

if __name__ == "__main__":
    main()
