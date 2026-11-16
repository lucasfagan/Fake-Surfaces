"""Search for a complexity-reducing sequence of moves.

Given a surface of complexity c, find the shortest sequence of T, U, T^-1,
U^-1 moves ending at a surface of complexity < c, staying at complexity
<= c + headroom.  Iterative deepening: lengths 1, 3, 5, ... (a reducing
sequence has odd length when it ends one step below c).  The search is only
a producer; the result is checked by replaying it (see verify.py).
"""

from __future__ import annotations

from .moves import MoveError, t_inverse, t_move, u_inverse, u_move, u_splits
from .surface import FakeSurface


class Budget(Exception):
    pass


def down_moves(s: FakeSurface, skip: int = 0):
    """All T^-1 and U^-1 moves, except on disk number `skip`."""
    for i, d in enumerate(s.disks):
        if i + 1 == skip:
            continue
        if len(d) == 3:
            try:
                yield ("T-", i + 1), t_inverse(s, i + 1)
            except MoveError:
                pass
        elif len(d) == 2:
            try:
                yield ("U-", i + 1), u_inverse(s, i + 1)
            except MoveError:
                pass


def up_moves(s: FakeSurface):
    """All T and U moves."""
    for e in s.edges:
        a, b = s.edge_endpoints(e)
        if a != b:
            yield ("T", e), t_move(s, e)
    for v in range(s.complexity):
        for a, b in u_splits(s, v):
            yield ("U", (a, b)), u_move(s, a, b)


def find_reduction(disks, plan=((3, 1), (5, 1), (5, 2), (7, 2)), budget=2_000_000,
                   goal=None, goal_complexity=None):
    """A shortest found sequence of moves from `disks` to a surface of lower
    complexity, or, if `goal` is given, to a surface of complexity
    `goal_complexity` satisfying `goal`.  `plan` lists (max length,
    headroom) stages.  Returns None if nothing is found within the budget
    (a number of expanded surfaces)."""
    s0 = FakeSurface(disks)
    c0 = s0.complexity
    if goal is None:
        def reached(r, x):
            return x < c0

        def need(x):  # moves still needed from complexity x
            return x - c0 + 1
        step = 2
    else:
        def reached(r, x):
            return x == goal_complexity and goal(r)

        def need(x):
            return max(1, abs(x - goal_complexity))
        step = 1
    left = [budget]

    def dfs(s, x, depth, cap, undo):
        """s has complexity x and `depth` moves remain; `undo` is the disk
        made by the previous up move (collapsing it would just go back)."""
        left[0] -= 1
        if left[0] < 0:
            raise Budget
        for m, r in down_moves(s, undo):
            if reached(r, x - 1):
                return [m]
            if depth - 1 >= need(x - 1):
                h = dfs(r, x - 1, depth - 1, cap, 0)
                if h:
                    return [m] + h
        if x < cap and depth - 1 >= need(x + 1):
            for m, r in up_moves(s):
                h = dfs(r, x + 1, depth - 1, cap, len(r.disks))
                if h:
                    return [m] + h
        return None

    try:
        for maxlen, head in plan:
            for length in range(1, maxlen + 1, step):
                h = dfs(s0, c0, length, c0 + head, 0)
                if h:
                    return h
    except Budget:
        return None
    return None


# --------------------------------------------------------------------------
# Best-first search, for the few surfaces the short search above misses.
# States are deduplicated up to homeomorphism with `search_key`, a canonical
# form that works at any complexity (unlike canonical.py, which needs the
# list of 1-skeleta); it is used only to prune the search.

def _variants(word):
    rev = tuple(-x for x in reversed(word))
    out = []
    for w in (word, rev):
        for i in range(len(w)):
            v = w[i:] + w[:i]
            if v not in out:
                out.append(v)
    return out


def _relabel(word, labels, next_id):
    labels = dict(labels)
    out = []
    for x in word:
        e = abs(x)
        if e not in labels:
            labels[e] = next_id if x > 0 else -next_id
            next_id += 1
        out.append(labels[e] if x > 0 else -labels[e])
    return tuple(out), labels, next_id


def search_key(s: FakeSurface):
    """Smallest relabeling of the disks by first encounter, over all disk
    orders, rotations and reversals (branching only on ties).  Equal for
    homeomorphic surfaces and different otherwise."""
    disks = s.disks
    n = len(disks)
    variants = [_variants(d) for d in disks]
    cls: dict = {}
    class_of = [cls.setdefault(frozenset(v), len(cls)) for v in variants]
    # an isomorphism-invariant label for each edge, to break ties early
    inv = {}
    for e in s.edges:
        a, b = s.edge_endpoints(e)
        inv[e] = (tuple(sorted(len(disks[di]) for di, _ in s.occurrences[e])), a == b)
    nbrs: dict[int, list[int]] = {e: [] for e in s.edges}
    for x, y, _ in s.corners:
        nbrs[x[0]].append(y[0])
        nbrs[y[0]].append(x[0])
    for _ in range(2):
        raw = {e: (inv[e], tuple(sorted(inv[o] for o in nbrs[e]))) for e in s.edges}
        rank = {v: i for i, v in enumerate(sorted(set(raw.values())))}
        inv = {e: rank[raw[e]] for e in s.edges}
    best = [None]

    def rec(used, labels, nid, form):
        if len(used) == n:
            if best[0] is None or form < best[0]:
                best[0] = form
            return
        if best[0] is not None and form > best[0][:len(form)]:
            return
        top, wins = None, []
        for di in range(n):
            if di in used:
                continue
            for v in variants[di]:
                lw, lab2, nid2 = _relabel(v, labels, nid)
                key = (nid2 - nid, len(lw), lw, tuple(inv[abs(x)] for x in v))
                if top is None or key < top:
                    top, wins = key, [(di, key, lab2, nid2)]
                elif key == top:
                    wins.append((di, key, lab2, nid2))
        seen = set()
        for di, key, lab2, nid2 in wins:
            sig = (class_of[di], tuple(sorted(lab2.items())))
            if sig not in seen:
                seen.add(sig)
                rec(used | {di}, lab2, nid2, form + (key,))

    rec(frozenset(), {}, 1, ())
    return tuple(w for _, _, w, _ in best[0])


def best_first_reduction(disks, headroom=3, budget=300_000):
    """Best-first search for a sequence of moves to lower complexity,
    preferring low complexity and few nontrivial T-bundles."""
    import heapq
    s0 = FakeSurface(disks)
    c0 = s0.complexity
    states = [(s0.disks, None, None)]  # (disks, parent index, move)
    seen = {search_key(s0)}
    heap = []
    count = 0

    def path(i):
        out = []
        while states[i][1] is not None:
            out.append(states[i][2])
            i = states[i][1]
        return out[::-1]

    def expand(i):
        nonlocal count
        s = FakeSurface(states[i][0])
        x = s.complexity
        for m, r in list(down_moves(s)) + (list(up_moves(s)) if x < c0 + headroom else []):
            count += 1
            heapq.heappush(heap, ((r.complexity, r.nontrivial_count()), count, i, m, r.disks))

    expand(0)
    while heap and len(states) < budget:
        (cx, _), _, parent, m, d = heapq.heappop(heap)
        r = FakeSurface(d)
        k = search_key(r)
        if k in seen:
            continue
        seen.add(k)
        states.append((d, parent, m))
        if cx < c0:
            return path(len(states) - 1)
        expand(len(states) - 1)
    return None
