#!/usr/bin/env python3
"""Check the classification and its certificates from scratch.

    python3 verify.py                  # everything except completeness
    python3 verify.py --max 5          # complexities 1..5 only (a few minutes)
    python3 verify.py --completeness   # also re-count every skeleton (slow at n = 6)

For every row of every table this checks that
  1. the attaching maps form a cellular fake surface (each edge on three
     sheets, each vertex link K4) with the stated complexity, a connected
     1-skeleton, and trivial integral homology;
  2. the row is in canonical form (fakesurfaces/canonical.py), on the stated
     1-skeleton, with the stated number of automorphisms; hence no two rows
     are homeomorphic;
  3. the embedded / trivial T-bundle columns and the count of nontrivial
     T-bundles are right;
  4. its certificate is a valid sequence of moves T, U, T^-1, U^-1 ending at
     the stated surface, which has lower complexity.
It then follows the certificates: every row marked contractible must reach
the abalone (a spine of the 3-ball, so a point), and every row marked not
contractible must reach a surface with an explicit nontrivial representation
of its fundamental group (data/noncontractible.json).

With --completeness it also proves that nothing is missing: on each
1-skeleton, the number of labeled acyclic surfaces (counted by an exhaustive
sweep over all sheet gluings, fakesurfaces/fastscan.c) must equal
sum over rows of |G| / |Aut F|, where |G| is the number of labelings.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections import defaultdict
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from fakesurfaces.canonical import canonical  # noqa: E402
from fakesurfaces.moves import MoveError, apply_move, parse_move  # noqa: E402
from fakesurfaces.skeleta import generate_skeleta, skeleta  # noqa: E402
from fakesurfaces.surface import FakeSurface  # noqa: E402
from fakesurfaces.table import DATA, POINT, read_rows, table_files  # noqa: E402

ABALONE = "1.1.1"


def digest(form) -> bytes:
    return hashlib.blake2b(repr(tuple(tuple(w) for w in form)).encode(), digest_size=12).digest()


class Failure(Exception):
    pass


def check(row):
    """Checks 1-4 for one row; returns (name, digest, target, end digest)."""
    def need(cond, what):
        if not cond:
            raise Failure(f"{row.name}: {what}")

    n, k = row.complexity, row.skeleton
    s = FakeSurface(row.disks)
    need(s.is_valid(), "not a fake surface")
    need(s.complexity == n and len(row.disks) == n + 1, "wrong complexity")
    need(s.is_acyclic(), "not acyclic")
    k2, form, aut = canonical(s)
    need(form == tuple(row.disks), "not in canonical form")
    need(k2 == k, f"lies on Gamma^{n}_{k2}, not Gamma^{n}_{k}")
    need(aut == row.automorphisms, "wrong number of automorphisms")
    emb = [s.is_embedded(i) for i in range(n + 1)]
    triv = [s.has_trivial_tbundle(i) for i in range(n + 1)]
    need(emb == row.embedded and triv == row.trivial, "wrong embedded/trivial columns")
    need(row.nontrivial == triv.count(False), "wrong nontrivial_bundles")
    end = None
    if row.reduces_to == POINT:
        # the abalone is a special spine of the 3-ball (Ikeda 1971)
        need(row.name == ABALONE and row.nontrivial == 0 and not row.certificate,
             "only the abalone is taken to a point directly")
    elif row.reduces_to:
        t = FakeSurface(row.disks)
        for text in row.certificate.split(";"):
            try:
                t = apply_move(t, parse_move(text))
            except MoveError as e:
                raise Failure(f"{row.name}: invalid move {text.strip()!r}: {e}")
        need(t.complexity < n or (n == 1 and row.reduces_to == ABALONE),
             "certificate does not reduce complexity")
        end = digest(canonical(t)[1])
    return row.name, digest(row.disks), row.reduces_to, end


def check_pi1(entry) -> None:
    """The images of the generators are permutations satisfying every
    relator, and not all trivial: so pi_1 is nontrivial."""
    rows = {r.name: r for r in read_rows(int(entry["surface"].split(".")[0]))}
    s = FakeSurface(rows[entry["surface"]].disks)
    gens, rels = s.presentation()
    img = {int(g): tuple(p) for g, p in entry["images"].items()}
    if sorted(img) != gens:
        raise Failure("generators do not match the presentation")
    deg = len(next(iter(img.values())))
    ident = tuple(range(deg))

    def mul(p, q):  # first p, then q
        return tuple(q[p[i]] for i in range(deg))

    def inv(p):
        out = [0] * deg
        for i, x in enumerate(p):
            out[x] = i
        return tuple(out)

    for w in rels:
        acc = ident
        for x in w:
            acc = mul(acc, img[x] if x > 0 else inv(img[-x]))
        if acc != ident:
            raise Failure(f"{entry['surface']}: relator {w} is not satisfied")
    if all(p == ident for p in img.values()):
        raise Failure("trivial representation")


def completeness(n: int, jobs: int):
    from fakesurfaces.enumerate import count_acyclic, labelings
    expect = defaultdict(Fraction)
    for row in iter_rows(n):
        expect[row.skeleton] += Fraction(labelings(skeleta(n)[row.skeleton - 1]), row.automorphisms)
    fixes = [()] if n <= 4 else [(a,) for a in range(6)] if n == 5 else \
        [(a, b) for a in range(6) for b in range(6)]
    tasks = [(k, m, f) for k, m in enumerate(skeleta(n), 1) for f in fixes]
    found = defaultdict(int)
    with Pool(jobs) as pool:
        for k, c in pool.imap_unordered(_count, tasks):
            found[k] += c
    for k in range(1, len(skeleta(n)) + 1):
        if found[k] != expect[k]:
            raise Failure(f"complexity {n}, Gamma_{k}: sweep finds {found[k]} labeled "
                          f"surfaces, the table accounts for {expect[k]}")
    print(f"  complexity {n}: complete ({sum(found.values())} labeled surfaces, "
          f"all accounted for on every 1-skeleton)", flush=True)


def _count(task):
    from fakesurfaces.enumerate import count_acyclic
    k, m, fix = task
    return k, count_acyclic(m, fix)


def iter_rows(n):
    yield from read_rows(n)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--max", type=int, default=6)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--completeness", action="store_true")
    a = ap.parse_args()

    # the 1-skeleta table
    for n in range(1, min(a.max, 5) + 1):
        if generate_skeleta(n) != skeleta(n):
            raise Failure(f"data/one_skeleta.csv is wrong at complexity {n}")
    print("1-skeleta: ok (regenerated for n <= 5)", flush=True)

    digests: dict[str, bytes] = {}
    targets: dict[str, tuple[str, bytes]] = {}
    contractible: dict[str, bool] = {}
    with Pool(a.jobs) as pool:
        for n in range(1, a.max + 1):
            t0 = time.time()
            names = set()
            forms = set()
            count = defaultdict(int)
            for name, dg, to, end in pool.imap(check, _tagged(iter_rows(n), contractible),
                                               chunksize=64):
                if name in names or dg in forms:
                    raise Failure(f"{name}: duplicate")
                names.add(name)
                forms.add(dg)
                digests[name] = dg
                if to:
                    targets[name] = (to, end)
                _, k, i = name.split(".")
                count[int(k)] += 1
                if int(i) != count[int(k)]:
                    raise Failure(f"{name}: surfaces are not numbered 1, 2, ... on Gamma_{k}")
            for name in names:
                if name in targets:
                    to, end = targets[name]
                    if to != POINT and (to not in digests or digests[to] != end):
                        raise Failure(f"{name}: certificate does not end at {to}")
            print(f"complexity {n}: {len(names)} surfaces ok ({time.time() - t0:.0f}s)", flush=True)

    # contractibility
    nc = [e for e in json.loads((DATA / "noncontractible.json").read_text())
          if int(e["surface"].split(".")[0]) <= a.max]
    for entry in nc:
        check_pi1(entry)
    roots = {e["surface"] for e in nc}

    def end_of(name):
        seen = []
        while name in targets and name not in roots:
            seen.append(name)
            name = targets[name][0]
        return name

    bad = 0
    for name, flag in contractible.items():
        e = end_of(name)
        if flag and e != POINT:
            raise Failure(f"{name} is marked contractible but its certificates end at {e}")
        if not flag and e not in roots:
            raise Failure(f"{name} is marked not contractible but its certificates end at {e}")
        bad += not flag
    print(f"contractibility: ok ({len(contractible) - bad} contractible, every one "
          f"3-deforms to a point; {bad} acyclic but not simply connected)", flush=True)

    if a.completeness:
        for n in range(1, a.max + 1):
            completeness(n, a.jobs)
    print("ALL CHECKS PASSED")


def _tagged(rows, contractible):
    for r in rows:
        contractible[r.name] = r.contractible
        yield r


if __name__ == "__main__":
    try:
        main()
    except Failure as e:
        sys.exit(f"FAILED: {e}")
