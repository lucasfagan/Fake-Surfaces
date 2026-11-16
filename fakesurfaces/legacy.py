"""Reader for the 2024 data file (legacy/fakesurfaces_2024.csv).

Its disk words are written as digit strings with a '-' before negative
entries, e.g. "211-2-1" = (2, 1, 1, -2, -1).  Edge 10 is written "10"; this
is unambiguous because 0 is never an edge label.
"""

import csv


def parse_word(s: str) -> tuple[int, ...]:
    out, i, sign = [], 0, 1
    while i < len(s):
        if s[i] == "-":
            sign = -1
            i += 1
            continue
        if s[i] == "1" and i + 1 < len(s) and s[i + 1] == "0":
            v, i = 10, i + 2
        else:
            v, i = int(s[i]), i + 1
        out.append(sign * v)
        sign = 1
    return tuple(out)


def read_legacy(path):
    """Yield (row number from 1, complexity, skeleton number, disks)."""
    with open(path) as f:
        for n, row in enumerate(csv.reader(f), start=1):
            words = [row[k] for k in range(2, len(row), 3) if row[k]]
            yield n, int(row[0]), int(row[1]), [parse_word(w) for w in words]
