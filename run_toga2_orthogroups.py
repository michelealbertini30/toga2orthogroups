#!/usr/bin/env python3
"""
toga2orthogroups — standalone entry point.

Delegates to src/toga2orthogroups.py, which contains the full implementation.

The launcher does not have to sit next to the code. The repository root (the
directory that contains src/) is resolved in this order:

    1. $TOGA2ORTHOGROUPS_HOME, if set
    2. the directory of the file this script really is, following symlinks —
       so a symlink dropped in a shared bin/ finds the clone it points into
    3. the directory this script sits in — the plain "run it from the clone"
       case

(2) and (3) exist so that on a shared filesystem the launcher can live in a
common scripts directory while src/ stays in a personal clone.

Usage:
    run_toga2_orthogroups.py -t DIR -s FILE -b FILE -i FILE -o DIR [options]
    run_toga2_orthogroups.py plot -t DIR -q SPECIES -g GENE -o OUTDIR [options]
"""
import os
import sys


def _repo_root() -> str:
    """Locate the directory holding src/."""
    env = os.environ.get("TOGA2ORTHOGROUPS_HOME")
    if env:
        return os.path.abspath(os.path.expanduser(env))

    here = os.path.dirname(os.path.abspath(__file__))
    real = os.path.dirname(os.path.realpath(__file__))
    for candidate in (real, here):
        if os.path.isdir(os.path.join(candidate, "src")):
            return candidate
    return here


ROOT = _repo_root()
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

try:
    from src.toga2orthogroups import main
except ImportError as exc:
    sys.exit(
        "toga2orthogroups: cannot import src/toga2orthogroups.py\n"
        f"  looked in: {ROOT}\n"
        f"  {exc}\n"
        "Point TOGA2ORTHOGROUPS_HOME at the directory containing src/, or\n"
        "symlink this launcher into your clone instead of copying it."
    )

if __name__ == "__main__":
    main()
