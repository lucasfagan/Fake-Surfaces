"""Enumeration of the cellular fake surfaces on a given 1-skeleton.

A fake surface structure on a 4-regular graph Gamma is exactly a choice,
for every edge e with ends x and y, of a bijection sigma_e between the three
corners at x that contain x and the three corners at y that contain y (a
corner is a pair of edge ends at a vertex; the vertex links are then
automatically K4).  Sheet k of e turns corner c_x[k] at one end and
sigma_e(c_x[k]) at the other.  The 2-cells are the orbits of walking along a
sheet, turning its corner, and continuing into the sheet of the next edge
that shares the corner.  A structure is cellular with t + 1 disks iff there
are exactly t + 1 orbits, and is then acyclic iff the paper's determinant
test passes.

So the cellular fake surfaces on Gamma, with Gamma's edges labeled and
oriented, are in bijection with the 6^(2t) assignments sigma that pass these
tests.  The sweep over all assignments runs in C (fastscan.c).  Edges are
labeled exactly as in the paper's convention (see skeleta.py): loops at a
vertex first, then edges in upper-triangle order, oriented from the earlier
vertex to the later one.
"""

from __future__ import annotations

import subprocess
import tempfile
from math import factorial
from pathlib import Path

from .skeleta import edge_slots, isomorphisms

HERE = Path(__file__).resolve().parent
SRC = HERE / "fastscan.c"
BIN = HERE.parent / "build" / "fastscan"


class SkeletonContext:
    """The edge ends ("slots") and corners of one labeled skeleton."""

    def __init__(self, mat):
        t = len(mat)
        self.t = t
        slots_at = [[] for _ in range(t)]
        nslots = 0

        def new_slot(v):
            nonlocal nslots
            slots_at[v].append(nslots)
            nslots += 1
            return nslots - 1

        # one (x, y) pair per edge, in the paper's edge order
        self.edges = [(new_slot(i), new_slot(j)) for i, j in edge_slots(mat)]
        self.vertex_of_slot = {s: v for v, ss in enumerate(slots_at) for s in ss}
        self.corner_id = {}
        for ss in slots_at:
            for a in range(4):
                for b in range(a + 1, 4):
                    self.corner_id[frozenset((ss[a], ss[b]))] = len(self.corner_id)
        # side_corners[m][side] = the 3 corners containing that end of edge m
        self.side_corners = []
        for x, y in self.edges:
            side = []
            for s in (x, y):
                v = self.vertex_of_slot[s]
                side.append([self.corner_id[frozenset((s, o))] for o in slots_at[v] if o != s])
            self.side_corners.append(side)

    def spec(self) -> str:
        """The skeleton data in the format read by fastscan."""
        E = len(self.edges)
        membership = [[] for _ in self.corner_id]
        for m in range(E):
            for s in (0, 1):
                for c in self.side_corners[m][s]:
                    membership[c].append((m, s))
        mem_index = {(c, m, s): 2 * c + j
                     for c, mem in enumerate(membership) for j, (m, s) in enumerate(mem)}
        other = {}
        for c, ((m1, s1), (m2, s2)) in enumerate(membership):
            other[(c, m1, s1)] = (m2, s2)
            other[(c, m2, s2)] = (m1, s1)
        blocks = [[], [], [], []]
        for m in range(E):
            rows = [[], [], [], []]
            for s in (0, 1):
                for k in range(3):
                    c = self.side_corners[m][s][k]
                    m2, s2 = other[(c, m, s)]
                    rows[0].append(mem_index[(c, m, s)])
                    rows[1].append(m2)
                    rows[2].append(s2)
                    rows[3].append(mem_index[(c, m2, s2)])
            for b in range(4):
                blocks[b].append(rows[b])
        everts = [[self.vertex_of_slot[x], self.vertex_of_slot[y]] for x, y in self.edges]
        lines = [f"{self.t} {E}"]
        for block in blocks + [everts]:
            lines += [" ".join(map(str, r)) for r in block]
        return "\n".join(lines) + "\n"


def ensure_built() -> Path:
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        BIN.parent.mkdir(exist_ok=True)
        subprocess.run(["cc", "-O2", "-o", str(BIN), str(SRC)], check=True)
    return BIN


def _run(mat, fix=(), extra=()):
    ensure_built()
    with tempfile.NamedTemporaryFile("w", suffix=".spec", delete=False) as f:
        f.write(SkeletonContext(mat).spec())
    cmd = [str(BIN), f.name, *extra]
    if fix:
        cmd += ["--fix", *map(str, fix)]
    return cmd, Path(f.name)


def count_acyclic(mat, fix=()) -> int:
    """Number of acyclic cellular fake surfaces on the labeled skeleton
    (labeled surfaces, not up to homeomorphism)."""
    cmd, spec = _run(mat, fix, ("--count-only", "--acyclic"))
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    finally:
        spec.unlink(missing_ok=True)
    return int(out.stderr.split()[3])


def labeled_acyclic(mat, fix=()):
    """Yield every acyclic labeled surface on the skeleton (disk lists)."""
    cmd, spec = _run(mat, fix, ("--acyclic",))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
                            bufsize=1 << 20)
    try:
        for line in proc.stdout:
            yield [[int(x) for x in d.split()] for d in line.split(";")]
    finally:
        proc.stdout.close()
        proc.wait()
        spec.unlink(missing_ok=True)


def labelings(mat) -> int:
    """|G|: the number of ways to label a surface on this skeleton in the
    paper's convention (graph automorphisms, orders of parallel edges and of
    loops, directions of loops).  A surface F has |G| / |Aut F| distinct
    labeled versions on the labeled skeleton."""
    n = len(mat)
    out = len(isomorphisms(mat, tuple(map(tuple, mat))))
    for i in range(n):
        loops = mat[i][i] // 2
        out *= factorial(loops) * 2 ** loops
        for j in range(i + 1, n):
            out *= factorial(mat[i][j])
    return out
