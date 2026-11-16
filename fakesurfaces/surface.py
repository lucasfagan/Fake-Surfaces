"""Cellular fake surfaces, stored as their 2-cell attaching maps.

A fake surface is a list of disks; each disk is a cyclic word of signed edge
labels (a negative entry runs along the edge against its orientation).  This
is the paper's notation (Fagan-Qiu-Wang, arXiv:2406.09439, Section 3.1):
the 1-skeleton, its vertices and the vertex links are all recovered from the
attaching maps.

* Every edge occurs exactly three times (three sheets meet along an edge).
* Every edge has two ends, (e, 0) = tail and (e, 1) = head; +e runs tail to
  head and -e runs head to tail.
* A corner of a disk (two cyclically consecutive entries x, y) identifies the
  end where x arrives with the end where y departs: the two ends meet at a
  vertex.  Vertices are the classes of edge ends under these identifications.
* The word list is a cellular fake surface iff every vertex has exactly four
  edge ends and its six corners form the complete graph K4 on them (the link
  of a vertex of a fake surface is the 1-skeleton of a tetrahedron).
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from typing import Iterable, Sequence

End = tuple[int, int]         # (edge, 0 = tail | 1 = head)
Occurrence = tuple[int, int]  # (disk index, position in the disk word)


def arrival_end(x: int) -> End:
    """The edge end reached after running along the signed entry x."""
    return (x, 1) if x > 0 else (-x, 0)


def departure_end(x: int) -> End:
    """The edge end the signed entry x starts from."""
    return (x, 0) if x > 0 else (-x, 1)


class InvalidSurface(ValueError):
    pass


class FakeSurface:
    """A cellular fake surface given by its attaching maps."""

    __slots__ = ("disks", "_c")

    def __init__(self, disks: Iterable[Sequence[int]]):
        self.disks: tuple[tuple[int, ...], ...] = tuple(
            tuple(int(x) for x in d) for d in disks)
        self._c: dict = {}
        if any(len(d) == 0 for d in self.disks) or any(
                x == 0 for d in self.disks for x in d):
            raise InvalidSurface("empty disk or edge label 0")

    def __repr__(self) -> str:
        return f"FakeSurface({[list(d) for d in self.disks]})"

    # ----------------------------------------------------------- basics
    @property
    def edges(self) -> tuple[int, ...]:
        if "edges" not in self._c:
            self._c["edges"] = tuple(sorted({abs(x) for d in self.disks for x in d}))
        return self._c["edges"]

    @property
    def occurrences(self) -> dict[int, list[Occurrence]]:
        """edge -> the (disk, position) pairs where it occurs."""
        if "occ" not in self._c:
            occ: dict[int, list[Occurrence]] = {}
            for di, d in enumerate(self.disks):
                for p, x in enumerate(d):
                    occ.setdefault(abs(x), []).append((di, p))
            self._c["occ"] = occ
        return self._c["occ"]

    @property
    def corners(self) -> list[tuple[End, End, Occurrence]]:
        """(arrival end, departure end, occurrence of the second entry) for
        every corner of every disk."""
        if "corners" not in self._c:
            cs = []
            for di, d in enumerate(self.disks):
                k = len(d)
                for i in range(k):
                    j = (i + 1) % k
                    cs.append((arrival_end(d[i]), departure_end(d[j]), (di, j)))
            self._c["corners"] = cs
        return self._c["corners"]

    @property
    def vertex_of_end(self) -> dict[End, int]:
        """Edge end -> vertex number (0, 1, ...), numbered by smallest end."""
        if "voe" not in self._c:
            parent = {}
            for e in self.edges:
                parent[(e, 0)] = (e, 0)
                parent[(e, 1)] = (e, 1)

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x

            for a, b, _ in self.corners:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
            roots = sorted({find(x) for x in parent})
            rid = {r: i for i, r in enumerate(roots)}
            self._c["voe"] = {x: rid[find(x)] for x in parent}
        return self._c["voe"]

    @property
    def vertices(self) -> list[list[End]]:
        """vertex -> its (sorted) edge ends."""
        if "verts" not in self._c:
            vs: list[list[End]] = [[] for _ in range(len(set(self.vertex_of_end.values())))]
            for end, v in self.vertex_of_end.items():
                vs[v].append(end)
            self._c["verts"] = [sorted(v) for v in vs]
        return self._c["verts"]

    @property
    def complexity(self) -> int:
        return len(self.vertices)

    @property
    def euler_characteristic(self) -> int:
        return self.complexity - len(self.edges) + len(self.disks)

    def edge_endpoints(self, e: int) -> tuple[int, int]:
        return self.vertex_of_end[(e, 0)], self.vertex_of_end[(e, 1)]

    # ------------------------------------------------------- validation
    def validate(self) -> None:
        """Raise InvalidSurface unless every edge has three sheets and every
        vertex link is K4."""
        for e, occ in self.occurrences.items():
            if len(occ) != 3:
                raise InvalidSurface(f"edge {e} occurs {len(occ)} times")
        pairs: dict[int, list[frozenset]] = {}
        for a, b, _ in self.corners:
            if a == b:
                raise InvalidSurface(f"disk backtracks at {a}")
            pairs.setdefault(self.vertex_of_end[a], []).append(frozenset((a, b)))
        for v, ends in enumerate(self.vertices):
            want = {frozenset(p) for p in combinations(ends, 2)}
            got = pairs.get(v, [])
            if len(ends) != 4 or len(got) != 6 or set(got) != want:
                raise InvalidSurface(f"link of vertex {v} is not K4")

    def is_valid(self) -> bool:
        try:
            self.validate()
            return True
        except InvalidSurface:
            return False

    # --------------------------------------------- skeleton and homology
    def spanning_tree(self) -> list[int]:
        """Edges of a maximal forest of the 1-skeleton (greedy, by label)."""
        parent = list(range(self.complexity))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        tree = []
        for e in self.edges:
            a, b = (find(v) for v in self.edge_endpoints(e))
            if a != b:
                parent[a] = b
                tree.append(e)
        return tree

    def is_connected(self) -> bool:
        return len(self.spanning_tree()) == self.complexity - 1

    def presentation(self) -> tuple[list[int], list[tuple[int, ...]]]:
        """Collapse a spanning tree: generators are the remaining edges,
        relators are the disks with tree edges deleted."""
        tree = set(self.spanning_tree())
        gens = [e for e in self.edges if e not in tree]
        return gens, [tuple(x for x in d if abs(x) not in tree) for d in self.disks]

    def is_acyclic(self) -> bool:
        """H_*(F) = H_*(point): chi = 1, connected 1-skeleton, and the
        square boundary matrix of the presentation has determinant +-1
        (the paper's test, Section 3.2.2)."""
        if self.euler_characteristic != 1 or not self.is_connected():
            return False
        gens, rels = self.presentation()
        idx = {g: i for i, g in enumerate(gens)}
        mat = [[0] * len(gens) for _ in rels]
        for r, w in zip(mat, rels):
            for x in w:
                r[idx[abs(x)]] += 1 if x > 0 else -1
        return len(mat) == len(gens) and abs(det(mat)) == 1

    # ---------------------------------------------- disks and T-bundles
    def boundary_vertices(self, i: int) -> list[int]:
        d = self.disks[i]
        return [self.vertex_of_end[arrival_end(x)] for x in d]

    def is_embedded(self, i: int) -> bool:
        """A 2-cell is an embedded disk iff its boundary passes through each
        vertex at most once (paper, Section 2.6)."""
        vs = self.boundary_vertices(i)
        return len(vs) == len(set(vs))

    def corner_partner(self, occ: Occurrence, end: End) -> End:
        """The other end of the corner that occurrence `occ` turns at `end`."""
        di, p = occ
        d = self.disks[di]
        if arrival_end(d[p]) == end:
            return departure_end(d[(p + 1) % len(d)])
        if departure_end(d[p]) == end:
            return arrival_end(d[(p - 1) % len(d)])
        raise ValueError(f"{occ} does not touch {end}")

    def prong_transport(self, i: int) -> tuple[list, list]:
        """Follow the two other sheets ("prongs") of the T-bundle around the
        boundary of disk i.  Returns (ring A, ring B): ringX[k] is the
        occurrence over boundary edge k reached from prong X at edge 0, for
        k = 0..len, the last entry being where the prong returns."""
        d = self.disks[i]
        k = len(d)

        def others(j):
            return [o for o in self.occurrences[abs(d[j])] if o != (i, j)]

        a, b = others(0)
        ra, rb = [a], [b]
        for j in range(k):
            nxt = (j + 1) % k
            pin = [self.corner_partner(ra[-1], arrival_end(d[j])),
                   self.corner_partner(rb[-1], arrival_end(d[j]))]
            cand = others(nxt)
            pout = [self.corner_partner(o, departure_end(d[nxt])) for o in cand]
            if set(pin) != set(pout) or len(set(pin)) != 2:
                raise InvalidSurface(f"inconsistent link around disk {i}")
            ra.append(cand[pout.index(pin[0])])
            rb.append(cand[pout.index(pin[1])])
        return ra, rb

    def has_trivial_tbundle(self, i: int) -> bool:
        """The T-bundle over the boundary of disk i is trivial (an annulus
        with a disk on its core) iff the two prongs come back to themselves,
        and nontrivial (a Moebius band with a disk on its core) iff they
        come back swapped."""
        ra, rb = self.prong_transport(i)
        if ra[-1] == ra[0]:
            return True
        assert ra[-1] == rb[0]
        return False

    def nontrivial_count(self) -> int:
        return sum(not self.has_trivial_tbundle(i) for i in range(len(self.disks)))


def det(mat: list[list[int]]) -> int:
    """Exact integer determinant (Gaussian elimination over Q)."""
    n = len(mat)
    m = [[Fraction(x) for x in row] for row in mat]
    out = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if m[r][c] != 0), None)
        if p is None:
            return 0
        if p != c:
            m[c], m[p] = m[p], m[c]
            out = -out
        out *= m[c][c]
        for r in range(c + 1, n):
            if m[r][c]:
                f = m[r][c] / m[c][c]
                for cc in range(c, n):
                    m[r][cc] -= f * m[c][cc]
    return int(out)
