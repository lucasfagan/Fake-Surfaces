"""Canonical form of a fake surface, written in the paper's convention.

Two cellular fake surfaces are PL homeomorphic iff their attaching maps are
related by relabeling edges (a graph automorphism), reversing edges,
reversing or rotating disk words, and reordering disks (arXiv:2406.09439,
Section 3.2.3).  The canonical form picks one representative:

1. the vertices are ordered so that the adjacency matrix is the paper's
   D(A)-maximal matrix of the surface's 1-skeleton Gamma^n_k;
2. edges are numbered by their entry (i, j), i <= j, of the upper triangle,
   and every non-loop edge runs from vertex i to vertex j, i < j;
3. among everything this leaves open (graph automorphisms, the order of
   parallel edges, the direction of loops), and over all rotations and
   reversals of each disk word and all orders of the disks, take the
   smallest result.  Words are compared entry by entry with
   1 < -1 < 2 < -2 < ...; disks are listed longest first.

Two surfaces are homeomorphic iff their canonical forms are equal.  The
number of choices in step 3 that produce the canonical form itself is the
order of the surface's automorphism group (as a cell complex), |Aut F|.
"""

from __future__ import annotations

from itertools import permutations, product

from .skeleta import edge_slots, identify, matrix_of_surface, skeleta
from .surface import FakeSurface

Form = tuple[tuple[int, ...], ...]


def _key(label: int, neg: bool) -> int:
    return 2 * label + neg


def _unkey(k: int) -> int:
    return -(k >> 1) if k & 1 else k >> 1


def _best_rotation(w: tuple[int, ...]) -> tuple[int, ...]:
    inv = tuple(k ^ 1 for k in reversed(w))
    n = len(w)
    return min(min(w[i:] + w[:i] for i in range(n)),
               min(inv[i:] + inv[:i] for i in range(n)))


def _labelings(s: FakeSurface, target, isos):
    """Yield every edge -> signed-label map compatible with the convention."""
    n = len(target)
    slots = edge_slots(target)
    first = {}
    for lab, ij in enumerate(slots, start=1):
        first.setdefault(ij, lab)
    for perm in isos:
        pos = {v: i for i, v in enumerate(perm)}
        groups: dict[tuple[int, int], list[tuple[int, int]]] = {}
        for e in s.edges:
            a, b = s.edge_endpoints(e)
            i, j = pos[a], pos[b]
            # sign +1 if e already runs from the earlier to the later vertex
            sign = 1 if i <= j else -1
            groups.setdefault((min(i, j), max(i, j)), []).append((e, sign))
        keys = sorted(groups)
        choices = []
        for ij in keys:
            g = groups[ij]
            if ij[0] == ij[1]:  # loops: any order, any direction
                choices.append([tuple(zip(p, sg)) for p in permutations([e for e, _ in g])
                                for sg in product((1, -1), repeat=len(g))])
            else:
                choices.append([tuple(p) for p in permutations(g)])
        for pick in product(*choices):
            lab = {}
            for ij, assignment in zip(keys, pick):
                for off, (e, sign) in enumerate(assignment):
                    lab[e] = (first[ij] + off) * sign
            yield lab


def canonical(s: FakeSurface) -> tuple[int, Form, int]:
    """(skeleton number k, canonical disk words, |Aut F|) for a surface on
    a connected 1-skeleton."""
    k, isos = identify(matrix_of_surface(s))
    target = skeleta(s.complexity)[k - 1]
    best, count = None, 0
    for lab in _labelings(s, target, isos):
        words = []
        for d in s.disks:
            w = []
            for x in d:
                y = lab[abs(x)] if x > 0 else -lab[abs(x)]
                w.append(_key(abs(y), y < 0))
            words.append(_best_rotation(tuple(w)))
        form = tuple(sorted(words, key=lambda w: (-len(w), w)))
        if best is None or form < best:
            best, count = form, 1
        elif form == best:
            count += 1
    return k, tuple(tuple(_unkey(x) for x in w) for w in best), count


def normalize(disks) -> Form:
    """The disks with each word rotated/reversed to its smallest version and
    the disks sorted, without relabeling any edge."""
    words = [_best_rotation(tuple(_key(abs(x), x < 0) for x in d)) for d in disks]
    words.sort(key=lambda w: (-len(w), w))
    return tuple(tuple(_unkey(x) for x in w) for w in words)


def is_canonical(disks) -> bool:
    return normalize(disks) == canonical(FakeSurface(disks))[1]
