# Errata for the 2024 classification

This lists what changed between the data published with
[arXiv:2406.09439](https://arxiv.org/abs/2406.09439) in 2024 (kept unchanged in
[`legacy/fakesurfaces_2024.csv`](legacy/fakesurfaces_2024.csv)) and the current
tables in [`data/`](data).

## 1. Missing surfaces

The 2024 tables are missing **62 surfaces**: 1 of complexity 3, 58 of complexity 4,
and 3 of complexity 5 without small disks. Every surface that was published is
correct: each of the 5,389 rows is a valid acyclic cellular fake surface, no two of
them are homeomorphic, and each one appears in the new tables
([`data/legacy_2024.csv`](data/legacy_2024.csv) gives its new name).

| complexity | 2024 | correct | missing |
|---|---|---|---|
| 1 | 2 | 2 | 0 |
| 2 | 17 | 17 | 0 |
| 3 | 238 | **239** | 1 (`3.3.3`) |
| 4 | 4,618 | **4,676** | 58 (on Γ⁴₃, Γ⁴₅, Γ⁴₆, Γ⁴₇, Γ⁴₈) |
| 5, no small disks | 514 | **517** | 3 (`5.19.18`, `5.28.26`, `5.28.102`) |

The missing surfaces are:

```
3.3.3
4.3.2 4.3.3 4.3.8 4.3.10 4.3.11 4.3.16 4.3.18 4.3.20 4.3.22 4.3.23 4.3.25 4.3.27
4.3.35 4.3.48 4.3.50 4.3.67 4.3.85 4.3.103 4.3.121 4.3.123 4.3.128 4.3.131 4.3.133
4.3.135 4.3.137 4.3.138 4.3.144 4.3.145 4.3.169 4.3.174 4.3.176 4.3.177 4.3.180
4.3.200 4.3.339 4.3.347
4.5.23 4.5.26 4.5.27 4.5.29 4.5.55 4.5.56 4.5.61 4.5.63 4.5.73 4.5.365 4.5.395
4.6.2 4.7.4 4.7.213 4.7.268 4.7.269
4.8.3 4.8.6 4.8.7 4.8.15 4.8.31 4.8.47
5.19.18 5.28.26 5.28.102
```

**Cause.** The 2024 code removed duplicates by checking whether two surfaces
were related by an automorphism of the 1-skeleton (`are_same` in
[`legacy/fakesurfaces_cla_6.py`](legacy/fakesurfaces_cla_6.py)). That check also
accepted some relabelings that are not automorphisms. So it sometimes declared two
non-homeomorphic surfaces equal and kept only one of them. All 62 missing surfaces
were discarded this way.

**How we know the new tables are complete.** They were produced by a different
enumeration (a sweep over all sheet gluings on each 1-skeleton, see the README),
and on every 1-skeleton they pass an exact counting check:

    number of labeled acyclic surfaces on Γ  =  Σ over rows on Γ of  |G_Γ| / |Aut F|,

where the left side is counted by brute force and |G_Γ| is the number of ways
to label Γ. A missing surface would make the right side too small.
`python3 verify.py --completeness` reruns this check.

Per 1-skeleton, complexity 3 now has 97, 65, 63 and 14 surfaces on Γ³₁, …, Γ³₄.
The 2024 data had 97, 65, 62 and 14. Note that Section 3.3.3 of the paper lists the
2024 counts in a different order ("14, 65, 97 and 62 coming from the first, second,
third, and fourth"). At complexity 4 the counts per 1-skeleton are in
[`data/summary_by_skeleton.csv`](data/summary_by_skeleton.csv). Γ⁴₂ is still the
largest (1,171) and Γ⁴₉ the smallest (35).

## 2. Table 1 of arXiv:2406.09439

Corrected counts, by number of disks with a nontrivial T-bundle:

| complexity | 0 (spines) | 1 | 2 | 3 | 4 | 5 | total |
|---|---|---|---|---|---|---|---|
| 1 | 1 | 1 | | | | | 2 |
| 2 | 3 | 6 | 6 | 2 | | | 17 |
| 3 | 20 | **55** | 89 | 62 | 13 | | **239** |
| 4 | 128 | **610** | **1,468** | **1,563** | **751** | **156** | **4,676** |

Bold entries changed. No spine was missing, so the spine counts and the
conclusions drawn from them are unaffected. The current tables for complexity
up to 6 are in [`data/summary.csv`](data/summary.csv) and in the README.

## 3. Statements that still hold

* **Corollary 2** (every acyclic cellular fake surface of complexity at most 4
  is contractible) holds for the missing surfaces too. All of them are
  contractible: their certificates end at a point.
* **Theorem 3** (every contractible fake surface of complexity at most 5 has an
  embedded disk) holds for the missing surfaces too, and is now checked through
  complexity 6 (see the README).
* **Theorem 1 of arXiv:2412.12293** (every contractible fake surface of
  complexity less than 6 with connected 1-skeleton 3-deforms to a point) was
  proved by checking the conditions of its Lemma 7 on the 2024 tables, and
  needed separate arguments for small disks and for disconnected 1-skeleta. The
  certificates here prove it directly for every surface, including the missing
  ones, and extend it to complexity 6.

## 4. Changes of format

* **Disk words** are written as space-separated integers (`1 2 -3`) instead of
  digit strings (`12-3`). The old form was hard to read once edge 10 appeared.
* **Edge orientations.** The paper fixes the vertex order and the edge
  numbering of each 1-skeleton (Section 3.1). The 2024 data followed that
  numbering, but oriented non-loop edges arbitrarily. Now every non-loop edge
  runs from the earlier vertex to the later one, so only loops carry a free
  orientation.
* **One representative per surface.** Each row is now in a canonical form, so
  two rows describe homeomorphic surfaces exactly when they are equal, and
  anyone can recompute a surface's row. The 2024 rows were not normalized, so
  most of them look different now even though they describe the same surfaces.
* **Names.** A surface is named `n.k.i`: the i-th surface on the 1-skeleton
  Γⁿₖ. The skeleton numbering is the paper's and is unchanged.
* **Complexity 5** is now complete. It used to cover only surfaces without small
  disks.
