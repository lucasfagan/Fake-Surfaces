#!/usr/bin/env python3
"""Property tests for the library, on random surfaces from the tables.

* every move produces a valid acyclic fake surface of the right complexity,
  and the disk it creates has the right type (T: embedded triangle, trivial
  T-bundle; U: embedded 2-gon, nontrivial T-bundle);
* T followed by T^-1 on the new triangle, and U followed by U^-1 on the new
  2-gon, give back the surface we started from;
* the canonical form does not change under random relabeling, reversal of
  edges, rotation and reversal of words, and reordering of disks.

    python3 scripts/selftest.py
"""

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fakesurfaces.canonical import canonical  # noqa: E402
from fakesurfaces.moves import t_inverse, t_move, u_inverse, u_move, u_splits  # noqa: E402
from fakesurfaces.surface import FakeSurface  # noqa: E402
from fakesurfaces.table import read_rows  # noqa: E402


def scramble(disks, rng):
    edges = sorted({abs(x) for d in disks for x in d})
    new = edges[:]
    rng.shuffle(new)
    m = {e: new[i] * rng.choice((1, -1)) for i, e in enumerate(edges)}
    out = []
    for d in disks:
        w = [m[abs(x)] * (1 if x > 0 else -1) for x in d]
        if rng.random() < 0.5:
            w = [-x for x in reversed(w)]
        r = rng.randrange(len(w))
        out.append(w[r:] + w[:r])
    rng.shuffle(out)
    return out


def main():
    rng = random.Random(2026)
    sample = []
    for n in (2, 3, 4, 5):
        rows = list(read_rows(n))
        sample += rng.sample(rows, min(150, len(rows)))
    moves = 0
    for row in sample:
        s = FakeSurface(row.disks)
        form = canonical(s)[1]
        assert form == tuple(row.disks)
        assert canonical(FakeSurface(scramble(row.disks, rng)))[1] == form
        for e in s.edges:
            a, b = s.edge_endpoints(e)
            if a == b:
                continue
            r = t_move(s, e)
            k = len(r.disks)
            r.validate()
            assert r.is_acyclic() and r.complexity == s.complexity + 1
            assert r.is_embedded(k - 1) and r.has_trivial_tbundle(k - 1)
            if r.complexity <= 6:
                assert canonical(t_inverse(r, k))[1] == form
            moves += 1
        for v in range(s.complexity):
            for a, b in u_splits(s, v):
                r = u_move(s, a, b)
                k = len(r.disks)
                r.validate()
                assert r.is_acyclic() and r.complexity == s.complexity + 1
                assert r.is_embedded(k - 1) and not r.has_trivial_tbundle(k - 1)
                assert canonical(u_inverse(r, k))[1] == form
                moves += 1
    print(f"ok: {len(sample)} surfaces, {moves} moves and their inverses")


if __name__ == "__main__":
    main()
