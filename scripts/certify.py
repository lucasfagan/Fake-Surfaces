#!/usr/bin/env python3
"""Step 3 of the pipeline: find a reduction certificate for every surface.

For each surface of complexity n (from build/canonical_<n>.jsonl) search for
a short sequence of moves T, U, T^-1, U^-1 ending at a surface of complexity
n - 1, replay it, and identify the end surface in the table of complexity
n - 1 by its canonical form.  Writes build/certificates_<n>.jsonl, one line
per surface in table order:

    [name, "move; move; ...", target name]      or  [name, null, null]

Resumable: rerunning with the same --complexity skips finished lines.

    python3 scripts/certify.py --complexity 5 --jobs 16
"""

import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fakesurfaces.canonical import canonical  # noqa: E402
from fakesurfaces.moves import apply_move, format_move  # noqa: E402
from fakesurfaces.search import best_first_reduction, find_reduction  # noqa: E402
from fakesurfaces.surface import FakeSurface  # noqa: E402
from fakesurfaces.table import read_canonical  # noqa: E402

BUILD = Path(__file__).resolve().parent.parent / "build"
PLAN = ((3, 1), (5, 1), (5, 2), (7, 2))


def work(arg):
    name, disks, budget, fallback = arg
    h = find_reduction(disks, plan=PLAN, budget=budget)
    if h is None and fallback:
        h = best_first_reduction(disks, headroom=3, budget=fallback)
    if h is None:
        return name, None, None
    s = FakeSurface(disks)
    for m in h:
        s = apply_move(s, m)
    k, form, _ = canonical(s)
    return name, "; ".join(format_move(*m) for m in h), (k, form)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complexity", type=int, required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--budget", type=int, default=200_000)
    ap.add_argument("--fallback", type=int, default=0,
                    help="budget of the best-first search tried when the short search fails")
    ap.add_argument("--fill", action="store_true",
                    help="only retry the surfaces that have no certificate yet")
    a = ap.parse_args()
    n = a.complexity
    if n == 1:
        return complexity_one()
    lower = {(k, form): name for name, k, form, _ in read_canonical(n - 1)}
    rows = [(name, form) for name, k, form, _ in read_canonical(n)]
    out = BUILD / f"certificates_{n}.jsonl"
    if a.fill:
        return fill(out, dict(rows), lower, a)
    done = sum(1 for _ in open(out)) if out.exists() else 0
    print(f"complexity {n}: {len(rows)} surfaces, {done} already done", flush=True)
    jobs = [(name, form, a.budget, a.fallback) for name, form in rows[done:]]
    t0 = time.time()
    fails = 0
    with open(out, "a") as f, Pool(a.jobs) as pool:
        for i, (name, moves, end) in enumerate(pool.imap(work, jobs, chunksize=32), 1):
            target = lower.get(end) if end else None
            if moves and target is None:
                raise SystemExit(f"{name}: end surface {end} is not in the table")
            fails += moves is None
            f.write(json.dumps([name, moves, target]) + "\n")
            if i % 20000 == 0:
                f.flush()
                print(f"  {done + i}/{len(rows)}  no certificate yet: {fails}  "
                      f"{time.time() - t0:.0f}s", flush=True)
    print(f"done: {fails} surfaces without a certificate", flush=True)


def complexity_one():
    """Nothing lies below complexity 1.  The abalone (the spine, all T-bundles
    trivial) is a spine of the 3-ball and goes to a point; the other surface
    gets a certificate taking it to the abalone."""
    rows = list(read_canonical(1))
    abalone = next(r for r in rows if FakeSurface(r[2]).nontrivial_count() == 0)
    certs = {abalone[0]: ["", "point"]}
    for name, k, form, _ in rows:
        if name == abalone[0]:
            continue
        h = find_reduction(form, plan=((4, 1), (6, 2)), goal_complexity=1,
                           goal=lambda r: canonical(r)[1] == abalone[2])
        certs[name] = ["; ".join(format_move(*m) for m in h), abalone[0]]
    (BUILD / "certificates_1.json").write_text(json.dumps(certs))
    print(certs)


def fill(out, forms, lower, a):
    lines = [json.loads(x) for x in open(out)]
    todo = [(name, forms[name], a.budget, a.fallback or 300_000)
            for name, moves, _ in lines if moves is None]
    print(f"retrying {len(todo)} surfaces", flush=True)
    found = {}
    with Pool(a.jobs) as pool:
        for name, moves, end in pool.imap_unordered(work, todo):
            if moves:
                found[name] = [name, moves, lower[end]]
            print(f"  {name}: {moves and lower[end]}", flush=True)
    with open(out, "w") as f:
        for line in lines:
            f.write(json.dumps(found.get(line[0], line)) + "\n")
    print(f"done: {sum(1 for x in lines if x[1] is None) - len(found)} surfaces without a certificate")


if __name__ == "__main__":
    main()
