# Classification of Acyclic Cellular Fake Surfaces

This repository contains

* the **complete classification of acyclic cellular fake surfaces of complexity 1 to 6**
  (3,083,732 surfaces), and
* a **certificate for each of them**: an explicit sequence of Matveev's moves T, U,
  T⁻¹, U⁻¹ taking it to a surface of lower complexity in the tables. The one
  exception is the Poincaré spine, which cannot have one. Following the certificates
  proves that **every contractible fake surface of complexity at most 6 3-deforms to
  a point**.

It accompanies L. Fagan, Y. Qiu and Z. Wang,
[*Classification of cellular fake surfaces*](https://arxiv.org/abs/2406.09439) and
[*Stable Andrews–Curtis conjecture via fake surfaces and Zeeman conjecture*](https://arxiv.org/abs/2412.12293).
The notation is the one used in those papers. The data published in 2024 missed 62
surfaces. **[ERRATA.md](ERRATA.md)** lists them, explains why they were missed, and
describes how the format has changed. Everything can be rechecked with one command
(see [Verifying](#verifying)).

## Results

| complexity | 1-skeleta | acyclic | contractible | spines | no small disks |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 2 | 2 | 1 | 0 |
| 2 | 2 | 17 | 17 | 3 | 0 |
| 3 | 4 | 239 | 239 | 20 | 2 |
| 4 | 10 | 4,676 | 4,676 | 128 | 31 |
| 5 | 28 | 111,460 | 111,459 | 1,298 | 517 |
| 6 | 97 | 2,967,338 | 2,967,333 | 13,663 | 7,762 |

(A *spine* is a surface all of whose disks have trivial T-bundle, i.e. a special spine
of a 3-manifold. A *small disk* is a disk of boundary length 1 or 2.) The breakdown
by number of nontrivial T-bundles, which extends Table 1 of the paper, is in
[`data/summary.csv`](data/summary.csv), and the counts per 1-skeleton are in
[`data/summary_by_skeleton.csv`](data/summary_by_skeleton.csv).

**Theorem 1 (classification).** Up to PL homeomorphism there are exactly 2, 17, 239,
4,676, 111,460 and 2,967,338 acyclic cellular fake surfaces of complexity 1, …, 6.
They are listed in [`data/`](data).

**Theorem 2 (reduction).** Every acyclic cellular fake surface of complexity n with
2 ≤ n ≤ 6, except the dodecahedral spine of the Poincaré homology sphere (`5.28.136`,
complexity 5), 3-deforms to an acyclic cellular fake surface of complexity n − 1 by
at most 19 moves T, U, T⁻¹, U⁻¹.

*Proof.* Theorem 1 holds because the tables pass the completeness check and
`verify.py` (see [Verifying](#verifying)). For Theorem 2, each row's certificate is an
explicit sequence of moves, and `verify.py` replays it and checks that it ends at the
row of complexity n − 1 named in `reduces_to`. Each move T, U, T⁻¹, U⁻¹ is an
elementary expansion of a 3-ball followed by an elementary collapse, so every
certificate is a 3-deformation. ∎

Everything else follows from these two theorems. Starting from any surface and
following the certificates lowers the complexity step by step until the chain stops at
one of two places:

* the **abalone** `1.1.1`, a special spine of the 3-ball, which therefore 3-deforms
  to a point (expand to the ball, collapse to a point). The only other surface of
  complexity 1, `1.1.2`, has a certificate taking it to the abalone; or
* the **Poincaré spine** `5.28.136`, whose fundamental group maps onto the
  alternating group A₅. The explicit representation is in
  [`data/noncontractible.json`](data/noncontractible.json) and is checked by
  `verify.py`.

3-deformations preserve homotopy type. So a surface is contractible exactly when its
chain ends at the abalone.

**Corollary 1.** Every contractible cellular fake surface of complexity at most 6
3-deforms to a point. Consequently the stable Andrews–Curtis conjecture holds for
every presentation of the trivial group read off from one of them (collapse a maximal
tree of the 1-skeleton), as in Corollary 2 of arXiv:2412.12293.

**Corollary 2.** Exactly 6 acyclic surfaces of complexity at most 6 are not
contractible: the Poincaré spine `5.28.136` and the five surfaces of complexity 6
whose chains end at it (`6.78.6891`, `6.88.4394`, `6.95.4721`, `6.95.5203`,
`6.97.1238`). Three of these five (`6.88.4394`, `6.95.4721`, `6.97.1238`) are spines:
special spines of complexity 6 of the Poincaré homology sphere minus a ball (by
geometrization, the only 3-manifold they can be spines of). In particular every
acyclic fake surface of complexity at most 4 is contractible (Corollary 2 of the
paper), and the Poincaré spine admits no complexity-reducing 3-deformation (one would
end at a contractible surface of complexity 4).

**Corollary 3.** Every acyclic cellular fake surface of complexity at most 6 has an
embedded disk (checked on the tables of Theorem 1; this extends Theorem 3 of the
paper).

This is the induction of arXiv:2412.12293, carried out for every surface explicitly.
That paper reduced a surface using its Lemma 7 (embedded disk moves) and treated small
disks and disconnected 1-skeleta separately. Here every step is a primitive move,
the moves T, U, T⁻¹, U⁻¹ never disconnect the 1-skeleton, and every surface on the
way down is a row of the tables.

## Reading the data

### Notation

As in the paper (Section 3.1), a fake surface is a pair (Γⁿₖ, M): a 1-skeleton and
the attaching maps of its disks.

* **1-skeleta.** The 1-skeleton of a surface of complexity n is a connected 4-regular
  multigraph on n vertices. Its adjacency matrix A has aᵢⱼ edges between vertices i
  and j and aᵢᵢ = twice the number of loops at i. The vertices are ordered so that
  D(A), the number whose digits are the entries of A read row by row, is as large as
  possible. Γⁿ₁, Γⁿ₂, … are the 1-skeleta of complexity n in decreasing order of D(A).
  All of them are listed in [`data/one_skeleta.csv`](data/one_skeleta.csv), with each
  matrix written row by row (e.g. Γ³₃ = `211 103 130`). This is the same numbering as
  in the paper and [`docs/ordered_one_skeletons.pdf`](docs/ordered_one_skeletons.pdf).
* **Edges** are numbered 1, …, 2n in the order of their entry (i, j), i ≤ j, in the
  upper triangle of A, read left to right and top to bottom. **Every non-loop edge
  runs from the earlier vertex i to the later vertex j.** (The paper says the
  orientations of non-loop edges are implied by the vertex order; the 2024 data did
  not actually follow a convention. Now it does.) Loops and the order of parallel
  edges are fixed by the canonical form below.
* **Attaching maps.** Each disk is a cyclic word in the edges; a negative entry
  traverses the edge backwards.
* **Names.** `n.k.i` is the i-th surface on Γⁿₖ. For example `2.1.9` is the 9th
  surface of complexity 2 on Γ²₁.

### Files

| file | contents |
|---|---|
| `data/complexity_1.csv` … `complexity_4.csv` | complexity 1–4 |
| `data/complexity_5.csv.gz` | complexity 5 (gzip) |
| `data/complexity_6/skeleton_01.csv.gz` … `skeleton_97.csv.gz` | complexity 6, one file per 1-skeleton (gzip) |
| `data/one_skeleta.csv` | the 1-skeleta |
| `data/summary.csv`, `data/summary_by_skeleton.csv` | counts |
| `data/noncontractible.json` | π₁ of the Poincaré spine maps onto A₅ |
| `data/legacy_2024.csv` | each row of the 2024 file and its new name |

A table of complexity n has one row per surface and these columns:

| column | meaning |
|---|---|
| `surface` | the name `n.k.i` |
| `disk_j` | attaching map of disk j (j = 1, …, n + 1), entries separated by spaces |
| `embedded_j` | `Y` if disk j is an embedded disk |
| `trivial_j` | `Y` if disk j has trivial T-bundle |
| `nontrivial_bundles` | number of disks with nontrivial T-bundle; 0 means a spine |
| `automorphisms` | order of the automorphism group of the surface |
| `contractible` | `Y` or `N` |
| `reduces_to` | the surface the certificate ends at (`point` for the abalone) |
| `certificate` | the moves, separated by `;` |

The first four kinds of column are the 2024 format. Disks are listed longest first.

### Example: Bing's house

Bing's house with two rooms is `2.1.9`. Its row in `data/complexity_2.csv` is

```
surface,disk_1,embedded_1,trivial_1,disk_2,embedded_2,trivial_2,disk_3,embedded_3,trivial_3,nontrivial_bundles,automorphisms,contractible,reduces_to,certificate
2.1.9,1 2 -3 1 3 4 -3 2 4 -2,N,Y,1,Y,Y,4,Y,Y,0,4,Y,1.1.1,T 2; U 1t 6t; U- 2; U- 2; U- 2
```

Its 1-skeleton is Γ²₁ = `22 22`: a loop at each vertex (edges 1 and 4) and two edges
2 and 3 from vertex 1 to vertex 2. The disks are the long disk, which passes through
each vertex more than once and so is not embedded, and the two walls `1` and `4`.
All three T-bundles are trivial, so it is a spine (of the 3-ball). Its certificate
goes up to complexity 3 and back down to the abalone. The paper's "mutant" of Bing's
house is `2.1.13`.

### Canonical form

Every row is in canonical form ([`fakesurfaces/canonical.py`](fakesurfaces/canonical.py)):
among all the ways of writing the surface in the convention above (automorphisms of
the 1-skeleton, the order of parallel edges, the directions of loops), all rotations
and reversals of the words, and all orders of the disks, it is the smallest. Words
are compared entry by entry with 1 < −1 < 2 < −2 < …, and disks are listed longest
first. Two surfaces are homeomorphic iff their canonical forms are equal, so to look
up a surface, put it in canonical form and search for the row:

```python
from fakesurfaces.surface import FakeSurface
from fakesurfaces.canonical import canonical

# Bing's house as written in the paper
k, form, aut = canonical(FakeSurface([[3, -1, -3, 2, -1, -2, -4, 2, -3, -4], [4], [-1]]))
# k = 1, form = ((1, 2, -3, 1, 3, 4, -3, 2, 4, -2), (1,), (4,)), aut = 4: row 2.1.9
```

## Certificates

A certificate is a sequence of moves, applied in order to the surface exactly as
written in its row:

| move | meaning | complexity |
|---|---|---|
| `T e` | the T move on the non-loop edge e (dual to the 2–3 Pachner move): the two endpoints of e become three vertices, joined by three new edges that bound a new triangle | +1 |
| `T- k` | T⁻¹: collapse disk k, an embedded triangle with trivial T-bundle (the trivial embedded disk move of arXiv:2412.12293 on a 3-vertex disk) | −1 |
| `U a b` | the U move at the vertex holding the edge ends a and b: the vertex splits in two, one keeping a and b, joined by two new edges that bound a new 2-gon. An edge end is written `5t` (tail of edge 5) or `5h` (head) | +1 |
| `U- k` | U⁻¹: collapse disk k, an embedded 2-gon on two distinct edges with nontrivial T-bundle (the nontrivial embedded disk move on a 2-vertex disk, Figure 12 of arXiv:2412.12293) | −1 |

Disks are numbered from 1 in the order listed. A move deletes what it destroys and
keeps everything else in order. New edges are labeled max + 1, max + 2, … and a new
disk goes at the end of the list. [`fakesurfaces/moves.py`](fakesurfaces/moves.py)
gives the exact rewriting rules. After the last move the surface is put in canonical
form, and it must equal the row named in `reduces_to`.

For example, `3.1.3` has certificate `U- 3` and reduces to `2.1.1`. `2.1.1` has
certificate `T 2; U- 3; U- 2` and reduces to the abalone `1.1.1`.

Lengths of the certificates:

| moves | 1 | 3 | 5 | 7 | 9–19 |
|---|---:|---:|---:|---:|---:|
| complexity 2 | 4 | 6 | 6 | 1 | |
| complexity 3 | 87 | 116 | 31 | 5 | |
| complexity 4 | 2,149 | 2,146 | 364 | 17 | |
| complexity 5 | 58,454 | 47,642 | 5,228 | 135 | |
| complexity 6 | 1,726,719 | 1,157,884 | 80,977 | 1,721 | 37 |

## Verifying

Requirements: Python 3.10 or later. Nothing else is needed, except a C compiler for
`--completeness`.

```sh
python3 verify.py --max 5 --jobs 16                 # complexities 1-5: under a minute
python3 verify.py --jobs 16                         # everything: about 10 minutes
python3 verify.py --jobs 16 --completeness          # plus completeness: about 1 hour more
```

For every row, `verify.py` checks that

1. the attaching maps form a cellular fake surface (every edge on three sheets, every
   vertex link K₄) with the stated complexity, a connected 1-skeleton and the
   homology of a point (the paper's determinant test);
2. the row is in canonical form on the stated 1-skeleton, with the stated number of
   automorphisms, and the rows are numbered 1, 2, … on each skeleton; so no two rows
   are homeomorphic;
3. the embedded and T-bundle columns are correct;
4. the certificate is a valid sequence of moves that ends, after canonicalization, at
   the row named in `reduces_to`, which has lower complexity. The two exceptions are
   at complexity 1: `1.1.2` goes to the abalone, and the abalone goes to `point`.

It then follows the certificates to check the `contractible` column, and checks the A₅
representation of the Poincaré spine.

**Completeness.** That nothing is missing is checked by counting. On a 1-skeleton Γ
with its edges labeled and oriented, a cellular fake surface is the same thing as a
choice, for each edge, of which of the three sheets at one end continues into which
sheet at the other end. That is 6^(2n) choices in all, and a C program
([`fakesurfaces/fastscan.c`](fakesurfaces/fastscan.c)) runs through all of them and
counts the acyclic ones. Each surface F on Γ appears |G_Γ| / |Aut F| times in that
count, where |G_Γ| is the number of ways to label Γ in the convention above. So the
table is complete iff, for every Γ,

    number of acyclic gluings on Γ  =  Σ over rows F on Γ of |G_Γ| / |Aut F|.

This holds for every 1-skeleton up to complexity 6. At complexity 5 there are
6,262,870 labeled acyclic surfaces, and at complexity 6 214,873,560.

## Reproducing the tables

```sh
python3 scripts/enumerate.py    --complexity 5 --jobs 16     # build/raw_5.jsonl
python3 scripts/canonicalize.py build/raw_5.jsonl --complexity 5 --jobs 16
python3 scripts/certify.py      --complexity 5 --jobs 16 --fallback 300000
python3 scripts/build_tables.py                              # writes data/
```

1. **Enumerate** ([`fakesurfaces/enumerate.py`](fakesurfaces/enumerate.py)): sweep
   all sheet gluings on every 1-skeleton and keep the acyclic surfaces that are
   already in canonical form, which yields each surface exactly once. This takes
   minutes up to complexity 5. Complexity 6 has about 3 × 10⁸ labeled surfaces; its
   table was produced by the same sweep with a C canonical form, and the completeness
   check above proves it correct independently of how it was produced.
2. **Canonicalize and deduplicate.**
3. **Certify** ([`fakesurfaces/search.py`](fakesurfaces/search.py)): an
   iterative-deepening search for a shortest reducing sequence (up to 7 moves, at
   most 2 above the starting complexity) certifies all but 37 surfaces of complexity
   6 in about 8 minutes on 16 cores. A best-first search going up to 3 above the
   starting complexity certifies the rest. The searches only produce certificates;
   `verify.py` is what checks them.
4. **Write the tables.** `scripts/legacy_map.py` matches the 2024 rows to their new
   names, and `scripts/selftest.py` runs property tests of the moves and the
   canonical form.

Complexity 1 is special because nothing lies below it. `scripts/certify.py
--complexity 1` marks the abalone as reducing to a point and searches for a
certificate taking `1.1.2` to the abalone.

The 1-skeleta are generated by `fakesurfaces.skeleta.generate_skeleta`, and
`verify.py` regenerates them up to complexity 5.

## Repository layout

```
README.md, ERRATA.md
verify.py                  checks everything
data/                      the classification and certificates
fakesurfaces/              a small library
  surface.py               fake surfaces from attaching maps; embedded disks; T-bundles; acyclicity
  moves.py                 T, U, T^-1, U^-1
  skeleta.py               1-skeleta in the paper's convention
  canonical.py             canonical form
  enumerate.py, fastscan.c the sweep over sheet gluings
  search.py                the certificate search
  table.py, legacy.py      reading and writing the tables
scripts/                   the pipeline that produced data/
docs/                      the 2024 notes on the convention and the ordered 1-skeleta
legacy/                    the 2024 data and code, unchanged
```

## Citation

```bibtex
@article{FQW-classification,
  author  = {Fagan, Lucas and Qiu, Yang and Wang, Zhenghan},
  title   = {Classification of cellular fake surfaces},
  journal = {arXiv preprint arXiv:2406.09439},
  year    = {2024}
}
@article{FQW-stable-AC,
  author  = {Fagan, Lucas and Qiu, Yang and Wang, Zhenghan},
  title   = {Stable {A}ndrews--{C}urtis conjecture via fake surfaces and {Z}eeman conjecture},
  journal = {arXiv preprint arXiv:2412.12293},
  year    = {2024}
}
```
