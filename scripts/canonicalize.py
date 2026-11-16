#!/usr/bin/env python3
"""Step 2 of the pipeline: canonicalize the raw enumeration.

Reads raw surfaces (one JSON list of disk words per line, or a JSON list of
them), checks each is a valid acyclic cellular fake surface, puts it in
canonical form, removes duplicates, and writes

    build/canonical_<n>.jsonl     one line per surface: [skeleton, disks, |Aut|]

sorted by skeleton number and then by canonical form.

    python3 scripts/canonicalize.py RAW_FILE --complexity 6 --jobs 16
"""

import argparse
import gzip
import json
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fakesurfaces.canonical import canonical  # noqa: E402
from fakesurfaces.surface import FakeSurface  # noqa: E402

BUILD = Path(__file__).resolve().parent.parent / "build"


def read_raw(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        if ".jsonl" not in str(path):
            yield from json.load(f)
        else:
            for line in f:
                if line.strip():
                    yield json.loads(line)


def sort_key(kf):
    """Skeleton, then the canonical order of words (1 < -1 < 2 < ...)."""
    k, form = kf
    return k, tuple(tuple(2 * abs(x) + (x < 0) for x in w) for w in form)


def work(disks):
    s = FakeSurface(disks)
    s.validate()
    if not s.is_acyclic():
        raise ValueError(f"not acyclic: {disks}")
    k, form, aut = canonical(s)
    return k, form, aut


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw")
    ap.add_argument("--complexity", type=int, required=True)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    seen = {}
    n = 0
    with Pool(a.jobs) as pool:
        for k, form, aut in pool.imap_unordered(work, read_raw(a.raw), chunksize=500):
            n += 1
            if (k, form) in seen:
                if seen[(k, form)] != aut:
                    raise SystemExit("inconsistent |Aut|")
                continue
            seen[(k, form)] = aut
            if n % 100000 == 0:
                print(f"  {n} read, {len(seen)} distinct", flush=True)
    BUILD.mkdir(exist_ok=True)
    out = BUILD / f"canonical_{a.complexity}.jsonl"
    with open(out, "w") as f:
        for (k, form) in sorted(seen, key=sort_key):
            f.write(json.dumps([k, [list(w) for w in form], seen[(k, form)]]) + "\n")
    print(f"complexity {a.complexity}: {n} read, {len(seen)} distinct -> {out}")


if __name__ == "__main__":
    main()
