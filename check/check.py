#!/usr/bin/env python3
# check_pdb_mapping_multi.py
# Verify a packed multi-component PDB is a concatenation of known monomer PDB atom orders.
#
# Usage:
#   python check_pdb_mapping_multi.py packed.pdb --template POLY=poly.pdb --template WAT=water.pdb
#   python check_pdb_mapping_multi.py packed.pdb --template A=a.pdb --template B=b.pdb --fields name
#   python check_pdb_mapping_multi.py packed.pdb --template ... --report-seq 50
#
# Notes:
# - Mapping assumption matches moltemplate.sh -pdb behavior: order-based coordinate overwrite.
# - This script checks if packed atoms can be parsed as repeated blocks of the provided templates.

from __future__ import annotations
import argparse
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional

PDB_REC = ("ATOM  ", "HETATM")

@dataclass
class PDBAtom:
    line_no: int
    record: str
    serial: str
    name: str
    resname: str
    chain: str
    resseq: str
    x: float
    y: float
    z: float
    element: str

def _strip(s: str) -> str:
    return s.strip()

def read_pdb_atoms(path: str) -> List[PDBAtom]:
    atoms: List[PDBAtom] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for idx, line in enumerate(f, start=1):
            if not line.startswith(PDB_REC):
                continue
            record = line[0:6]
            serial = _strip(line[6:11])
            name = _strip(line[12:16])
            resname = _strip(line[17:20])
            chain = _strip(line[21:22])
            resseq = _strip(line[22:26])
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                x = y = z = float("nan")
            element = ""
            if len(line) >= 78:
                element = _strip(line[76:78])
            if not element and name:
                element = "".join([c for c in name if c.isalpha()])[:2].upper()
            atoms.append(PDBAtom(
                line_no=idx, record=record, serial=serial, name=name,
                resname=resname, chain=chain, resseq=resseq,
                x=x, y=y, z=z, element=element
            ))
    return atoms

def make_sig(atoms: List[PDBAtom], fields: List[str], ignore_case: bool) -> List[Tuple[str, ...]]:
    out: List[Tuple[str, ...]] = []
    for a in atoms:
        parts: List[str] = []
        for f in fields:
            if f == "name":
                parts.append(a.name)
            elif f == "resname":
                parts.append(a.resname)
            elif f == "element":
                parts.append(a.element)
            elif f == "chain":
                parts.append(a.chain)
            elif f == "resseq":
                parts.append(a.resseq)
            else:
                raise ValueError(f"Unknown field: {f}")
        if ignore_case:
            parts = [p.upper() for p in parts]
        out.append(tuple(parts))
    return out

def fmt_atom(a: PDBAtom) -> str:
    return (f"line {a.line_no}: {a.record.strip()} serial={a.serial} "
            f"name={a.name} res={a.resname} chain={a.chain} resseq={a.resseq} "
            f"elem={a.element} xyz=({a.x:.3f},{a.y:.3f},{a.z:.3f})")

def first_mismatch(a: List[Tuple[str, ...]], b: List[Tuple[str, ...]]) -> Optional[int]:
    n = min(len(a), len(b))
    for i in range(n):
        if a[i] != b[i]:
            return i
    if len(a) != len(b):
        return n
    return None

def main():
    ap = argparse.ArgumentParser(description="Check packed PDB can be decomposed into repeated blocks of multiple monomer templates (order mapping for moltemplate -pdb).")
    ap.add_argument("packed_pdb", help="multi-component packed pdb (e.g., packmol output)")
    ap.add_argument("--template", action="append", required=True,
                    help="template in NAME=path.pdb form. Repeatable. Example: --template POLY=poly.pdb")
    ap.add_argument("--fields", nargs="+", default=["name", "resname", "element"],
                    help="fields used for matching: name resname element chain resseq. Default: name resname element")
    ap.add_argument("--ignore-case", action="store_true", help="case-insensitive compare")
    ap.add_argument("--report-seq", type=int, default=40, help="print first N recognized blocks in sequence (default 40)")
    ap.add_argument("--max-fail-context", type=int, default=5, help="print +/-K atoms around first failure (default 5)")
    args = ap.parse_args()

    # Load packed
    pk_atoms = read_pdb_atoms(args.packed_pdb)
    if not pk_atoms:
        raise SystemExit(f"[ERROR] No ATOM/HETATM records found in {args.packed_pdb}")
    pk_sig = make_sig(pk_atoms, args.fields, args.ignore_case)

    # Load templates
    tmpl_atoms: Dict[str, List[PDBAtom]] = {}
    tmpl_sig: Dict[str, List[Tuple[str, ...]]] = {}
    tmpl_len: Dict[str, int] = {}

    for t in args.template:
        if "=" not in t:
            raise SystemExit(f"[ERROR] Bad --template '{t}'. Use NAME=path.pdb")
        name, path = t.split("=", 1)
        name = name.strip()
        path = path.strip()
        atoms = read_pdb_atoms(path)
        if not atoms:
            raise SystemExit(f"[ERROR] No ATOM/HETATM records in template {name}: {path}")
        sig = make_sig(atoms, args.fields, args.ignore_case)
        tmpl_atoms[name] = atoms
        tmpl_sig[name] = sig
        tmpl_len[name] = len(sig)

    # Greedy parse from start
    i = 0
    seq: List[str] = []
    counts: Dict[str, int] = {k: 0 for k in tmpl_len.keys()}

    # For disambiguation: try longer templates first
    tmpl_order = sorted(tmpl_len.keys(), key=lambda k: tmpl_len[k], reverse=True)

    print("=== Multi-component PDB Mapping Check ===")
    print(f"Packed atoms: {len(pk_atoms)}  ({args.packed_pdb})")
    print(f"Compare fields: {args.fields}  ignore_case={args.ignore_case}")
    print("Templates:")
    for k in tmpl_order:
        print(f"  - {k}: {tmpl_len[k]} atoms")

    while i < len(pk_sig):
        matched = None
        for k in tmpl_order:
            L = tmpl_len[k]
            if i + L > len(pk_sig):
                continue
            if pk_sig[i:i+L] == tmpl_sig[k]:
                matched = k
                break
        if matched is None:
            print(f"\n[FAIL] Cannot match any template at packed atom index {i} (PDB line {pk_atoms[i].line_no}).")
            print("Packed atom at failure:")
            print("  ", fmt_atom(pk_atoms[i]))
            print("\nTry reasons:")
            print("  - packed.pdb order is interleaved or altered (not block-concatenated templates)")
            print("  - one template is missing / wrong (different atom order, different naming, different resname/element)")
            print("  - compare fields too strict; try --fields name (looser)")

            # Context
            K = args.max_fail_context
            lo = max(0, i-K)
            hi = min(len(pk_atoms), i+K+1)
            print(f"\nContext packed atoms [{lo}:{hi}]:")
            for j in range(lo, hi):
                mark = ">>" if j == i else "  "
                print(mark, fmt_atom(pk_atoms[j]))
            raise SystemExit(2)

        # record match
        seq.append(matched)
        counts[matched] += 1
        i += tmpl_len[matched]

    print("\n[PASS] Packed PDB can be decomposed into concatenated template blocks (order-safe for moltemplate -pdb).")
    print("Counts:")
    for k in tmpl_order:
        print(f"  {k}: {counts[k]}")

    # Print first N blocks sequence (to see ordering, grouping)
    nshow = min(args.report_seq, len(seq))
    if nshow > 0:
        print(f"\nFirst {nshow} blocks in order:")
        print("  " + " ".join(seq[:nshow]))
        if len(seq) > nshow:
            print(f"  ... (total blocks: {len(seq)})")

    # Also report whether blocks are grouped (packmol-style) or interleaved
    # Simple heuristic: count transitions
    trans = sum(1 for a, b in zip(seq, seq[1:]) if a != b)
    print(f"\nTransitions between component types: {trans} (0 means fully grouped by type).")

if __name__ == "__main__":
    main()
