"""Reading and writing the classification tables.

Layout (see README.md):

    data/complexity_1.csv ... data/complexity_4.csv
    data/complexity_5.csv.gz
    data/complexity_6/skeleton_01.csv.gz ... skeleton_97.csv.gz

Columns of a complexity-n table (n + 1 disks):

    surface                  name n.k.i: the i-th surface on 1-skeleton Gamma^n_k
    disk_j                   attaching map of disk j, entries separated by spaces
    embedded_j, trivial_j    Y/N: disk j is embedded / has trivial T-bundle
    nontrivial_bundles       number of disks with nontrivial T-bundle (0 = spine)
    automorphisms            order of the automorphism group of the surface
    contractible             Y/N
    reduces_to               the surface the certificate ends at
    certificate              the moves, separated by ';' (see moves.py)
"""

from __future__ import annotations

import csv
import gzip
import io
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUILD = ROOT / "build"
DATA = ROOT / "data"

POINT = "point"


def name_of(n: int, k: int, i: int) -> str:
    """Surface i on 1-skeleton Gamma^n_k."""
    return f"{n}.{k}.{i}"


def parse_name(name: str) -> tuple[int, int, int]:
    n, k, i = name.split(".")
    return int(n), int(k), int(i)


def read_canonical(n: int):
    """Yield (name, skeleton, canonical form, |Aut|) from the build step."""
    counter: dict[int, int] = {}
    with open(BUILD / f"canonical_{n}.jsonl") as f:
        for line in f:
            k, form, aut = json.loads(line)
            counter[k] = counter.get(k, 0) + 1
            yield name_of(n, k, counter[k]), k, tuple(tuple(w) for w in form), aut


def header(n: int) -> list[str]:
    cols = ["surface"]
    for j in range(1, n + 2):
        cols += [f"disk_{j}", f"embedded_{j}", f"trivial_{j}"]
    return cols + ["nontrivial_bundles", "automorphisms", "contractible",
                   "reduces_to", "certificate"]


def word_str(w) -> str:
    return " ".join(map(str, w))


def table_files(n: int) -> list[Path]:
    if n <= 4:
        return [DATA / f"complexity_{n}.csv"]
    if n == 5:
        return [DATA / "complexity_5.csv.gz"]
    return sorted((DATA / f"complexity_{n}").glob("skeleton_*.csv.gz"))


def _open(path: Path, mode: str):
    if path.suffix == ".gz":
        if "w" in mode:
            # mtime=0 keeps the compressed files byte-for-byte reproducible
            return io.TextIOWrapper(gzip.GzipFile(path, "wb", mtime=0), newline="")
        return gzip.open(path, "rt", newline="")
    return open(path, mode, newline="")


class Row:
    __slots__ = ("name", "disks", "embedded", "trivial", "nontrivial",
                 "automorphisms", "contractible", "reduces_to", "certificate")

    @property
    def complexity(self) -> int:
        return parse_name(self.name)[0]

    @property
    def skeleton(self) -> int:
        return parse_name(self.name)[1]


def read_rows(n: int, path: Path | None = None):
    """Yield the rows of the complexity-n table (or of one of its files)."""
    for p in [path] if path else table_files(n):
        with _open(p, "r") as f:
            r = csv.reader(f)
            next(r)
            for rec in r:
                row = Row()
                row.name = rec[0]
                m = (len(rec) - 6) // 3
                row.disks = [tuple(int(x) for x in rec[1 + 3 * j].split()) for j in range(m)]
                row.embedded = [rec[2 + 3 * j] == "Y" for j in range(m)]
                row.trivial = [rec[3 + 3 * j] == "Y" for j in range(m)]
                row.nontrivial = int(rec[1 + 3 * m])
                row.automorphisms = int(rec[2 + 3 * m])
                row.contractible = rec[3 + 3 * m] == "Y"
                row.reduces_to = rec[4 + 3 * m]
                row.certificate = rec[5 + 3 * m]
                yield row


def write_rows(path: Path, n: int, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with _open(path, "w") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header(n))
        for row in rows:
            rec = [row.name]
            for d, e, t in zip(row.disks, row.embedded, row.trivial):
                rec += [word_str(d), "Y" if e else "N", "Y" if t else "N"]
            rec += [row.nontrivial, row.automorphisms, "Y" if row.contractible else "N",
                    row.reduces_to, row.certificate]
            w.writerow(rec)
