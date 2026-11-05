#!/usr/bin/env python3
"""Shared reader for the stage-1 instrumentation dumps (stage1/dumps/<comp>/*.csv).

Dump layout: a block of `# key = value` header lines (always `# description`
and `# columns`, the comma-separated column names), then fixed-width CSV rows
written by Fortran: values are space-padded, strings may be quoted, logicals
are T/F, many files end each row with a trailing comma (one more field than
`# columns`), NaN appears as `NaN`, and floats are es26.17 (`1.23E+00`, `1.23E+001`; the old gfortran
form without the exponent letter, `1.23-315`, is also accepted).

A plain csv.DictReader does not handle this header; use this module instead:

    from dumpcsv import load, select
    meta, rows = load("stage1/dumps/clm/clm_surfalbedo.csv")
    meta["columns"]                      # list of column names
    rows[0]["albgrd_vis"]                # float
    select(rows, kind="nourban", step=2) # filter by exact column values

Values are converted to int, float, bool (T/F), or stripped str. Floats keep
full double precision, because Python's float() round-trips es26.17.

Command line, for a quick look:

    python3 stage1/analysis/dumpcsv.py <file> [-n 5] [--where kind=urban step=2] [--cols a,b,c]
"""
import re
import sys

_FORTRAN_NO_E = re.compile(r"^([+-]?\d*\.\d+)([+-]\d{3})$")


def convert(tok):
    s = tok.strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in "\"'":
        return s[1:-1].strip()
    if s in ("T", ".TRUE.", ".true."):
        return True
    if s in ("F", ".FALSE.", ".false."):
        return False
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        pass
    m = _FORTRAN_NO_E.match(s)
    if m:
        return float(f"{m.group(1)}E{m.group(2)}")
    return s


def _split(line):
    # Fortran rows never embed commas inside quoted strings, but keep quotes intact.
    out, cur, q = [], [], None
    for ch in line:
        if q:
            cur.append(ch)
            if ch == q:
                q = None
        elif ch in "\"'":
            q = ch
            cur.append(ch)
        elif ch == ",":
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur))
    return out


def load(path):
    """Return (meta, rows): meta maps header keys to strings (plus 'columns' as a
    list); rows is a list of dicts keyed by column name."""
    meta, rows, cols = {}, [], None
    with open(path) as f:
        for ln, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                key, sep, val = line[1:].partition("=")
                if sep:
                    meta[key.strip()] = val.strip()
                continue
            if cols is None:
                if "columns" not in meta:
                    raise ValueError(f"{path}: no '# columns =' header before data")
                cols = [c.strip() for c in meta["columns"].split(",")]
                meta["columns"] = cols
            toks = _split(line)
            if len(toks) == len(cols) + 1 and not toks[-1].strip():
                toks.pop()  # many dumps end each row with a trailing comma
            if len(toks) != len(cols):
                raise ValueError(f"{path}:{ln}: {len(toks)} values for {len(cols)} columns")
            rows.append({c: convert(t) for c, t in zip(cols, toks)})
    if cols is None and "columns" in meta:
        meta["columns"] = [c.strip() for c in meta["columns"].split(",")]
    return meta, rows


def select(rows, **eq):
    """Rows whose columns equal every given value (strings compared stripped)."""
    return [r for r in rows if all(r.get(k) == v for k, v in eq.items())]


def _main(argv):
    import argparse
    ap = argparse.ArgumentParser(description="Inspect a stage-1 dump CSV.")
    ap.add_argument("file")
    ap.add_argument("-n", type=int, default=5, help="rows to print (default 5)")
    ap.add_argument("--where", nargs="*", default=[], metavar="COL=VAL")
    ap.add_argument("--cols", default=None, help="comma-separated columns to print")
    a = ap.parse_args(argv)
    meta, rows = load(a.file)
    for w in a.where:
        k, _, v = w.partition("=")
        rows = select(rows, **{k: convert(v)})
    print(f"# {a.file}: {len(rows)} rows, {len(meta['columns'])} columns")
    print("# description:", meta.get("description", "")[:300])
    cols = a.cols.split(",") if a.cols else meta["columns"]
    print("# columns:", ",".join(cols))
    for r in rows[: a.n]:
        print(",".join(repr(r[c]) for c in cols))


if __name__ == "__main__":
    _main(sys.argv[1:])
