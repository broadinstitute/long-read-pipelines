#!/usr/bin/env python3
"""Numeric field-by-field comparison of two hmmIBD-style TSVs.
Strings must match exactly; numbers within rel tol 1e-4 (abs 1e-6)."""
import sys

def isnum(s):
    try:
        float(s); return True
    except ValueError:
        return False

def cmp(fa, fb, tol=1e-4, atol=1e-6):
    A = [l.rstrip("\n").split("\t") for l in open(fa)]
    B = [l.rstrip("\n").split("\t") for l in open(fb)]
    if len(A) != len(B):
        print(f"ROW COUNT DIFF {fa}:{len(A)} vs {fb}:{len(B)}"); return False
    ok = True
    diffs = 0
    maxrel = 0.0
    for i, (ra, rb) in enumerate(zip(A, B)):
        if len(ra) != len(rb):
            print(f"  line {i}: ncol {len(ra)} vs {len(rb)}"); ok = False; continue
        for j, (a, b) in enumerate(zip(ra, rb)):
            if a == b:
                continue
            if isnum(a) and isnum(b):
                fa_, fb_ = float(a), float(b)
                denom = max(abs(fa_), abs(fb_), 1.0)
                rel = abs(fa_ - fb_) / denom
                maxrel = max(maxrel, rel)
                if abs(fa_ - fb_) > atol and rel > tol:
                    if diffs < 10:
                        print(f"  line {i} col {j}: {a} vs {b} (rel {rel:.2e})")
                    diffs += 1; ok = False
            else:
                if diffs < 10:
                    print(f"  line {i} col {j}: STR {a!r} vs {b!r}")
                diffs += 1; ok = False
    print(f"{'MATCH' if ok else 'MISMATCH'}: {fa} vs {fb}  (max rel diff {maxrel:.2e}, {diffs} field diffs)")
    return ok

a = cmp(sys.argv[1], sys.argv[2])
sys.exit(0 if a else 1)
