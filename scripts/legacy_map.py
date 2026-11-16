#!/usr/bin/env python3
"""Match every row of the 2024 data file to its name in the new tables.

Writes data/legacy_2024.csv and prints the surfaces the 2024 tables missed
(complexity <= 4, and complexity 5 without small disks).

    python3 scripts/legacy_map.py
"""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from fakesurfaces.canonical import canonical  # noqa: E402
from fakesurfaces.legacy import read_legacy  # noqa: E402
from fakesurfaces.surface import FakeSurface  # noqa: E402
from fakesurfaces.table import read_rows  # noqa: E402


def main():
    names = {}
    for n in range(1, 6):
        for r in read_rows(n):
            names[tuple(r.disks)] = r.name
    used = set()
    with open(ROOT / "data" / "legacy_2024.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["row_2024", "complexity", "skeleton", "surface"])
        for row, c, k, disks in read_legacy(ROOT / "legacy" / "fakesurfaces_2024.csv"):
            name = names[canonical(FakeSurface(disks))[1]]
            used.add(name)
            w.writerow([row, c, k, name])
    missing = [name for form, name in names.items() if name not in used
               and (not name.startswith("5.") or all(len(d) >= 3 for d in form))]
    print(f"{len(used)} rows matched; missing from 2024: {len(missing)}")
    print(" ".join(missing))


if __name__ == "__main__":
    main()
