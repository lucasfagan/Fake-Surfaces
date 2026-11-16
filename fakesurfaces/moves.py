"""Matveev's moves T, U and their inverses on cellular fake surfaces.

Two fake surfaces 3-deform to each other iff they are related by a sequence
of these moves (Matveev, Algorithmic Topology and Classification of
3-Manifolds, Theorem 1.3.8).  Each move is an elementary expansion of a
3-ball followed by an elementary collapse, so every sequence of moves is a
3-deformation.

* T e     the T move on a non-loop edge e (dual to the 2-3 Pachner move):
          the two endpoints of e become three vertices, joined pairwise by
          three new edges that bound a new triangular disk (embedded, trivial
          T-bundle).  Complexity +1.
* T- k    the inverse move: collapse the k-th disk, which must be an
          embedded triangle with trivial T-bundle.  Complexity -1.  This is
          the trivial embedded disk move of arXiv:2412.12293 on a 3-vertex
          disk.
* U a b   the U move at the vertex holding the edge ends a and b: that
          vertex splits into two, one keeping a, b and the other the two
          remaining ends, joined by two new parallel edges that bound a new
          2-gon (embedded, nontrivial T-bundle).  Complexity +1.
* U- k    the inverse move: collapse the k-th disk, which must be an
          embedded 2-gon on two distinct edges with nontrivial T-bundle.
          Complexity -1.  This is the nontrivial embedded disk move of
          arXiv:2412.12293 on a 2-vertex disk (its Figure 12).

Conventions, so that a sequence of moves can be replayed exactly:
disks are numbered from 1 in the order they are listed; an edge end is
written `5t` (tail of edge 5) or `5h` (head of edge 5); new edges get the
labels max+1, max+2, ... in the order described in each function; a disk
created by a move is appended at the end of the list; deleted disks and
entries are removed and everything else keeps its order.
"""

from __future__ import annotations

from .surface import End, FakeSurface, InvalidSurface, Occurrence, arrival_end, departure_end


class MoveError(ValueError):
    """The move does not apply at the given site."""


class _Rewrite:
    """Edits to the attaching maps, applied in one pass."""

    def __init__(self, s: FakeSurface):
        self.s = s
        self.deleted: set[Occurrence] = set()
        self.inserted: dict[Occurrence, list[int]] = {}
        self.deleted_disks: set[int] = set()
        self.new_disks: list[tuple[int, ...]] = []

    def delete(self, occ: Occurrence) -> None:
        if occ in self.deleted:
            raise MoveError(f"entry {occ} deleted twice")
        self.deleted.add(occ)

    def insert_before(self, occ: Occurrence, syms: list[int]) -> None:
        """Insert symbols into the corner just before position occ."""
        self.inserted.setdefault(occ, []).extend(syms)

    def result(self) -> FakeSurface:
        out = []
        for di, d in enumerate(self.s.disks):
            if di in self.deleted_disks:
                continue
            w: list[int] = []
            for p, x in enumerate(d):
                w.extend(self.inserted.get((di, p), ()))
                if (di, p) not in self.deleted:
                    w.append(x)
            if not w:
                raise MoveError(f"the move would empty disk {di + 1}")
            out.append(tuple(w))
        out.extend(self.new_disks)
        return FakeSurface(out)


def _fresh(s: FakeSurface, n: int) -> list[int]:
    m = max(s.edges)
    return [m + i + 1 for i in range(n)]


def _corner(s: FakeSurface, x: End, y: End) -> tuple[End, Occurrence]:
    """The unique corner turning between ends x and y: (its arrival end,
    the occurrence just after it)."""
    for a, b, occ in s.corners:
        if {a, b} == {x, y}:
            return a, occ
    raise MoveError(f"no corner between {x} and {y}")


# ------------------------------------------------------------------ T
def t_move(s: FakeSurface, edge: int) -> FakeSurface:
    """T move on the non-loop edge `edge`.  New edges g1, g2, g3 (in that
    order); the new triangle is (g1, g2, g3)."""
    if edge not in s.occurrences:
        raise MoveError(f"no edge {edge}")
    if s.vertex_of_end[(edge, 0)] == s.vertex_of_end[(edge, 1)]:
        raise MoveError(f"edge {edge} is a loop")
    occs = s.occurrences[edge]
    tail = [s.corner_partner(o, (edge, 0)) for o in occs]
    head = [s.corner_partner(o, (edge, 1)) for o in occs]
    g = _fresh(s, 3)
    rw = _Rewrite(s)
    for o in occs:
        rw.delete(o)
    rw.new_disks.append(tuple(g))
    # sheet i of the old edge becomes the corner region w_i; the corner
    # {slot_j, slot_k} at either end is rerouted along g_i, which runs
    # from w_j to w_k for (i, j, k) cyclic.
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        for slots in (tail, head):
            arr, occ = _corner(s, slots[j], slots[k])
            rw.insert_before(occ, [g[i] if arr == slots[j] else -g[i]])
    return rw.result()


def t_inverse(s: FakeSurface, k: int) -> FakeSurface:
    """T^-1 on disk k (1-based).  One new edge."""
    i = k - 1
    if not 0 <= i < len(s.disks):
        raise MoveError(f"no disk {k}")
    d = s.disks[i]
    if len(d) != 3 or not s.is_embedded(i) or not s.has_trivial_tbundle(i):
        raise MoveError(f"disk {k} is not an embedded triangle with trivial T-bundle")
    ra, rb = s.prong_transport(i)
    new = _fresh(s, 1)[0]
    rw = _Rewrite(s)
    rw.deleted_disks.add(i)
    for j in range(3):
        rw.delete(ra[j])
        rw.delete(rb[j])
    # at the corner vertex after boundary edge j the two rings turn
    # corners into slots a_j and b_j; the sheet between them gains the new
    # edge, which runs from the ring-A side to the ring-B side.
    for j in range(3):
        a = s.corner_partner(ra[j], arrival_end(d[j]))
        b = s.corner_partner(rb[j], arrival_end(d[j]))
        arr, occ = _corner(s, a, b)
        rw.insert_before(occ, [new if arr == a else -new])
    return rw.result()


# ------------------------------------------------------------------ U
def parse_end(tok: str) -> End:
    return (int(tok[:-1]), 0 if tok[-1] == "t" else 1)


def format_end(end: End) -> str:
    return f"{end[0]}{'t' if end[1] == 0 else 'h'}"


def u_move(s: FakeSurface, a: End, b: End) -> FakeSurface:
    """U move splitting the ends a, b off their vertex.  New edges h1, h2,
    both running from the new vertex holding a, b to the other one; the new
    2-gon is (h1, -h2)."""
    v = s.vertex_of_end.get(a)
    if v is None or a == b or s.vertex_of_end.get(b) != v:
        raise MoveError(f"{format_end(a)}, {format_end(b)} are not two ends at one vertex")
    c, d = [x for x in s.vertices[v] if x not in (a, b)]
    h1, h2 = _fresh(s, 2)
    rw = _Rewrite(s)
    rw.new_disks.append((h1, -h2))
    for (x, y), h in (((a, c), h1), ((b, d), h1), ((a, d), h2), ((b, c), h2)):
        arr, occ = _corner(s, x, y)
        rw.insert_before(occ, [h if arr in (a, b) else -h])
    return rw.result()


def u_splits(s: FakeSurface, v: int) -> list[tuple[End, End]]:
    """The three ways to split vertex v, each named by the pair holding the
    smallest end."""
    ends = s.vertices[v]
    return [(ends[0], x) for x in ends[1:]]


def u_inverse(s: FakeSurface, k: int) -> FakeSurface:
    """U^-1 on disk k (1-based)."""
    i = k - 1
    if not 0 <= i < len(s.disks):
        raise MoveError(f"no disk {k}")
    d = s.disks[i]
    if (len(d) != 2 or abs(d[0]) == abs(d[1]) or not s.is_embedded(i)
            or s.has_trivial_tbundle(i)):
        raise MoveError(f"disk {k} is not an embedded 2-gon with nontrivial T-bundle")
    rw = _Rewrite(s)
    rw.deleted_disks.add(i)
    for p in (0, 1):
        for occ in s.occurrences[abs(d[p])]:
            if occ != (i, p):
                rw.delete(occ)
    return rw.result()


# --------------------------------------------------------- notation
def format_move(kind: str, arg) -> str:
    if kind == "U":
        return f"U {format_end(arg[0])} {format_end(arg[1])}"
    return f"{kind} {arg}"


def parse_move(text: str) -> tuple[str, object]:
    parts = text.split()
    if parts[0] == "U":
        return "U", (parse_end(parts[1]), parse_end(parts[2]))
    if parts[0] in ("T", "T-", "U-"):
        return parts[0], int(parts[1])
    raise MoveError(f"unknown move {text!r}")


def apply_move(s: FakeSurface, move: tuple[str, object]) -> FakeSurface:
    kind, arg = move
    try:
        if kind == "T":
            out = t_move(s, arg)
        elif kind == "T-":
            out = t_inverse(s, arg)
        elif kind == "U":
            out = u_move(s, *arg)
        elif kind == "U-":
            out = u_inverse(s, arg)
        else:
            raise MoveError(f"unknown move {kind}")
        out.validate()
    except InvalidSurface as e:
        raise MoveError(str(e)) from e
    return out
