#!/usr/bin/env python3
"""Step 4 of the pipeline: write the tables in data/.

Combines build/canonical_<n>.jsonl (the surfaces) and
build/certificates_<n>.jsonl (their reductions) into the published tables,
computes the embedded / T-bundle columns, and decides contractibility by
following the certificates down to complexity 1.

    python3 scripts/build_tables.py
"""

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fakesurfaces.surface import FakeSurface  # noqa: E402
from fakesurfaces.table import (BUILD, DATA, POINT, Row, read_canonical,  # noqa: E402
                                write_rows)

MAX = 6


def main():
    rows: dict[int, list[Row]] = {}
    target: dict[str, str] = {}
    nonreducers = []
    for n in range(1, MAX + 1):
        certs = {}
        if n >= 2:
            with open(BUILD / f"certificates_{n}.jsonl") as f:
                for line in f:
                    name, moves, to = json.loads(line)
                    certs[name] = (moves, to)
        else:
            # complexity 1: the abalone 1.1.1 is a spine of the 3-ball, so it
            # 3-deforms to a point; 1.1.2 has a certificate taking it to 1.1.1
            certs = json.loads((BUILD / "certificates_1.json").read_text())
        out = []
        for name, k, form, aut in read_canonical(n):
            s = FakeSurface(form)
            r = Row()
            r.name, r.disks, r.automorphisms = name, list(form), aut
            r.embedded = [s.is_embedded(i) for i in range(len(form))]
            r.trivial = [s.has_trivial_tbundle(i) for i in range(len(form))]
            r.nontrivial = r.trivial.count(False)
            moves, to = certs[name]
            r.certificate, r.reduces_to = moves or "", to or ""
            if to:
                target[name] = to
            else:
                nonreducers.append(name)
            out.append(r)
        rows[n] = out
        print(f"complexity {n}: {len(out)} surfaces", flush=True)

    # contractible iff the chain of certificates ends at a point
    status: dict[str, bool | None] = {}

    def contractible(name):
        chain = []
        while name not in status:
            chain.append(name)
            if name == POINT:
                status[name] = True
                break
            if name not in target:
                status[name] = None  # decided separately
                break
            name = target[name]
        for c in chain:
            status[c] = status[name]
        return status[name]

    noncontractible = {e["surface"] for e in json.loads((DATA / "noncontractible.json").read_text())}
    for name in nonreducers:
        status[name] = False if name in noncontractible else None
    undecided = []
    for n in rows:
        for r in rows[n]:
            c = contractible(r.name)
            if c is None:
                undecided.append(r.name)
            r.contractible = bool(c)
    if undecided:
        raise SystemExit(f"undecided: {undecided[:20]} ({len(undecided)})")

    for n, out in rows.items():
        if n <= 4:
            write_rows(DATA / f"complexity_{n}.csv", n, out)
        elif n == 5:
            write_rows(DATA / "complexity_5.csv.gz", n, out)
        else:
            by_k = defaultdict(list)
            for r in out:
                by_k[r.skeleton].append(r)
            for k, rs in sorted(by_k.items()):
                write_rows(DATA / f"complexity_{n}" / f"skeleton_{k:02d}.csv.gz", n, rs)

    # summary: counts by complexity and number of nontrivial T-bundles
    with open(DATA / "summary.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["complexity", "acyclic", "contractible", "spines"]
                   + [f"nontrivial_{j}" for j in range(MAX + 2)] + ["no_small_disk"])
        for n, out in rows.items():
            ntb = Counter(r.nontrivial for r in out)
            w.writerow([n, len(out), sum(r.contractible for r in out), ntb[0]]
                       + [ntb[j] for j in range(MAX + 2)]
                       + [sum(all(len(d) >= 3 for d in r.disks) for r in out)])
    with open(DATA / "summary_by_skeleton.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["complexity", "skeleton", "acyclic", "contractible", "spines"])
        for n, out in rows.items():
            c = defaultdict(lambda: [0, 0, 0])
            for r in out:
                c[r.skeleton][0] += 1
                c[r.skeleton][1] += r.contractible
                c[r.skeleton][2] += r.nontrivial == 0
            for k in sorted(c):
                w.writerow([n, k, *c[k]])
    print("tables written to", DATA)


if __name__ == "__main__":
    main()
