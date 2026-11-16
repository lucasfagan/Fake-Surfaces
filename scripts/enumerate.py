#!/usr/bin/env python3
"""Step 1 of the pipeline: enumerate acyclic fake surfaces of complexity n.

Sweeps every sheet gluing on every 1-skeleton (fakesurfaces/fastscan.c) and
keeps the labeled surfaces that are already in canonical form: each surface
appears exactly once that way, so no separate deduplication is needed.
Writes build/raw_<n>.jsonl, the input of scripts/canonicalize.py.

Pure Python canonicalization of every labeled surface is fast enough for
n <= 5 (minutes).  For n = 6 there are about 3 * 10^8 labeled surfaces; the
published table was produced by the same sweep with a C canonical form,
and verify.py --completeness proves the table complete independently.

    python3 scripts/enumerate.py --complexity 4 --jobs 16
"""

import argparse
import json
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fakesurfaces.canonical import canonical, normalize  # noqa: E402
from fakesurfaces.enumerate import labeled_acyclic  # noqa: E402
from fakesurfaces.skeleta import skeleta  # noqa: E402
from fakesurfaces.surface import FakeSurface  # noqa: E402

BUILD = Path(__file__).resolve().parent.parent / "build"


def work(task):
    n, k, fix = task
    out = []
    for disks in labeled_acyclic(skeleta(n)[k - 1], fix):
        form = canonical(FakeSurface(disks))[1]
        if normalize(disks) == form:
            out.append([list(w) for w in form])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complexity", type=int, required=True)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    n = a.complexity
    fixes = [()] if n <= 4 else [(x,) for x in range(6)] if n == 5 else \
        [(x, y) for x in range(6) for y in range(6)]
    tasks = [(n, k, f) for k in range(1, len(skeleta(n)) + 1) for f in fixes]
    BUILD.mkdir(exist_ok=True)
    total = 0
    with open(BUILD / f"raw_{n}.jsonl", "w") as f, Pool(a.jobs) as pool:
        for found in pool.imap_unordered(work, tasks):
            for form in found:
                f.write(json.dumps(form) + "\n")
            total += len(found)
    print(f"complexity {n}: {total} surfaces -> build/raw_{n}.jsonl")


if __name__ == "__main__":
    main()
