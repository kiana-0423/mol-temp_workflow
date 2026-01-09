#!/usr/bin/env python3
# check_data_pdb_massmap.py
# Compare LAMMPS data atom order vs PDB atom order using element inference from Masses.
#
# Usage:
#   python check_data_pdb_massmap.py system.data packed.pdb
#   python check_data_pdb_massmap.py system.data packed.pdb --order line
#   python check_data_pdb_massmap.py system.data packed.pdb --tol 0.6 --max-report 30
#
# Notes:
# - Default assumes atom_style full-like where atom-ID and atom-type are first and third columns in Atoms section.
# - "Element" from data is inferred by matching type mass to nearest atomic weight within tolerance.

from __future__ import annotations
import argparse
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

PDB_REC = ("ATOM  ", "HETATM")

# Common atomic weights (amu). Add more if you need.
ATOMIC_WEIGHTS: Dict[str, float] = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "P": 30.974,
    "S": 32.06,
    "CL": 35.45,
    "BR": 79.904,
    "I": 126.904,

    "LI": 6.94,
    "NA": 22.990,
    "K": 39.098,
    "MG": 24.305,
    "CA": 40.078,

    "SI": 28.085,
    "AL": 26.982,
    "FE": 55.845,
    "CU": 63.546,
    "ZN": 65.38,
    "ZR": 91.224,
    "HF": 178.49,
    "Y": 88.906,
    "LA": 138.905,
    "TI": 47.867,
    "NI": 58.693,
    "CO": 58.933,
    "MN": 54.938,
    "CR": 51.996,
}

def _is_section_header(line: str) -> bool:
    # crude but works for typical data files: section headers contain letters and no leading numeric token
    s = line.strip()
    if not s:
        return False
    # common section names
    for k in ("Masses", "Atoms", "Bonds", "Angles", "Dihedrals", "Impropers",
              "Velocities", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs"):
        if s.startswith(k):
            return True
    return False

def _strip_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()

@dataclass
class DataAtom:
    atom_id: int
    atom_type: int
    line_no: int

@dataclass
class PDBAtom:
    line_no: int
    name: str
    element: str

def read_lammps_data_masses_and_atoms(path: str) -> Tuple[Dict[int, float], List[DataAtom]]:
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = list(f)

    masses: Dict[int, float] = {}
    atoms: List[DataAtom] = []

    # find "Masses" section
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith("Masses"):
            i += 1
            # skip blank/comment lines
            while i < len(lines) and not _strip_comment(lines[i]):
                i += 1
            # parse until next section header
            while i < len(lines):
                raw = lines[i]
                if _is_section_header(raw):
                    break
                s = _strip_comment(raw)
                if s:
                    parts = s.split()
                    # expecting: type mass
                    if len(parts) >= 2 and parts[0].lstrip("+-").isdigit():
                        t = int(parts[0])
                        try:
                            m = float(parts[1])
                            masses[t] = m
                        except ValueError:
                            pass
                i += 1
            continue
        i += 1

    # find "Atoms" section
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith("Atoms"):
            i += 1
            # skip blank/comment lines
            while i < len(lines) and not _strip_comment(lines[i]):
                i += 1
            # parse until next section header
            while i < len(lines):
                raw = lines[i]
                if _is_section_header(raw):
                    break
                s = _strip_comment(raw)
                if s:
                    parts = s.split()
                    # Most atom styles: first token is atom-ID (int). For full: third token is atom-type.
                    # We'll minimally require at least 3 tokens.
                    if len(parts) >= 3 and parts[0].lstrip("+-").isdigit():
                        try:
                            atom_id = int(parts[0])
                            atom_type = int(parts[2])
                            atoms.append(DataAtom(atom_id=atom_id, atom_type=atom_type, line_no=i+1))
                        except ValueError:
                            pass
                i += 1
            continue
        i += 1

    if not masses:
        raise RuntimeError("No Masses section parsed (or empty).")
    if not atoms:
        raise RuntimeError("No Atoms section parsed (or empty).")

    return masses, atoms

def infer_element_from_mass(mass: float, tol: float) -> Optional[Tuple[str, float]]:
    # return (element_symbol, delta_mass) for closest match within tol
    best_el = None
    best_delta = float("inf")
    for el, aw in ATOMIC_WEIGHTS.items():
        d = abs(mass - aw)
        if d < best_delta:
            best_delta = d
            best_el = el
    if best_el is None or best_delta > tol:
        return None
    return best_el, best_delta

def read_pdb_atoms(path: str) -> List[PDBAtom]:
    out: List[PDBAtom] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for idx, line in enumerate(f, start=1):
            if not line.startswith(PDB_REC):
                continue
            name = line[12:16].strip()
            element = ""
            if len(line) >= 78:
                element = line[76:78].strip()
            if not element:
                # infer from atom name: take leading letters, up to 2 chars
                letters = "".join([c for c in name if c.isalpha()])
                element = letters[:2].upper() if letters else ""
            element = element.upper()
            out.append(PDBAtom(line_no=idx, name=name, element=element))
    if not out:
        raise RuntimeError("No ATOM/HETATM records parsed from PDB.")
    return out

def main():
    ap = argparse.ArgumentParser(description="Check LAMMPS data vs PDB atom order by mapping atom_type masses to elements.")
    ap.add_argument("data_file", help="LAMMPS data file (contains Masses and Atoms sections)")
    ap.add_argument("pdb_file", help="PDB file (ATOM/HETATM records)")
    ap.add_argument("--order", choices=["id", "line"], default="id",
                    help="Atom ordering for data Atoms: 'id' sorts by atom-ID (default); 'line' uses file line order")
    ap.add_argument("--tol", type=float, default=0.6,
                    help="Mass-to-element matching tolerance in amu (default 0.6)")
    ap.add_argument("--max-report", type=int, default=20,
                    help="Max mismatches to print (default 20)")
    ap.add_argument("--strict", action="store_true",
                    help="If set, treat unknown element inference as mismatch too (default: unknowns are warned but not counted as mismatches)")
    args = ap.parse_args()

    masses, data_atoms = read_lammps_data_masses_and_atoms(args.data_file)
    pdb_atoms = read_pdb_atoms(args.pdb_file)

    if args.order == "id":
        data_atoms_sorted = sorted(data_atoms, key=lambda a: a.atom_id)
    else:
        data_atoms_sorted = list(data_atoms)

    n_data = len(data_atoms_sorted)
    n_pdb = len(pdb_atoms)

    print("=== data ↔ pdb order check (mass → element) ===")
    print(f"Data atoms: {n_data}  ({args.data_file})  order={args.order}")
    print(f"PDB atoms : {n_pdb}  ({args.pdb_file})")
    print(f"Mass tol  : {args.tol} amu")
    print(f"Strict    : {args.strict}")

    if n_data != n_pdb:
        print(f"[FAIL] Atom count mismatch: data={n_data} vs pdb={n_pdb}")
        return 2

    # Build type -> inferred element
    type_to_inferred: Dict[int, Tuple[Optional[str], float, float]] = {}
    # stores: (element or None, mass, delta)
    unknown_types = set()

    for t, m in masses.items():
        inf = infer_element_from_mass(m, args.tol)
        if inf is None:
            type_to_inferred[t] = (None, m, float("inf"))
            unknown_types.add(t)
        else:
            el, d = inf
            type_to_inferred[t] = (el, m, d)

    if unknown_types:
        print(f"[WARN] {len(unknown_types)} atom types could NOT be mapped to an element within tol={args.tol} amu.")
        # show a few
        show = sorted(list(unknown_types))[:10]
        for t in show:
            el, m, d = type_to_inferred[t]
            print(f"  type {t}: mass={m} -> element=UNKNOWN")
        if len(unknown_types) > 10:
            print("  ...")

    mismatches = 0
    unknown_used = 0
    reported = 0

    for i in range(n_data):
        da = data_atoms_sorted[i]
        pa = pdb_atoms[i]

        inf_el, m, d = type_to_inferred.get(da.atom_type, (None, float("nan"), float("inf")))
        if inf_el is None:
            unknown_used += 1
            if args.strict:
                mismatches += 1
                if reported < args.max_report:
                    reported += 1
                    print(f"\n[MISMATCH-UNKNOWN] index={i+1}")
                    print(f"  data: atom_id={da.atom_id} type={da.atom_type} mass={m} -> UNKNOWN  (data line {da.line_no})")
                    print(f"  pdb : element={pa.element} name={pa.name} (pdb line {pa.line_no})")
            continue

        # normalize elements like "CL"
        pdb_el = pa.element.upper()
        data_el = inf_el.upper()

        if pdb_el != data_el:
            mismatches += 1
            if reported < args.max_report:
                reported += 1
                print(f"\n[MISMATCH] index={i+1}")
                print(f"  data: atom_id={da.atom_id} type={da.atom_type} mass={m} -> {data_el} (Δ={d:.3g} amu)  (data line {da.line_no})")
                print(f"  pdb : element={pdb_el} name={pa.name} (pdb line {pa.line_no})")

    print("\n=== Summary ===")
    print(f"Unknown inferred elements used: {unknown_used} atoms (strict={args.strict})")
    print(f"Mismatches: {mismatches} / {n_data}")

    if mismatches == 0 and (not args.strict or unknown_used == 0):
        print("[PASS] data atom order is consistent with PDB order under this mass→element mapping.")
        return 0
    else:
        print("[FAIL/WARN] There are mismatches (or unknowns in strict mode).")
        print("Tips:")
        print("  - If you see many UNKNOWN types, increase --tol slightly or extend ATOMIC_WEIGHTS with your elements.")
        print("  - If 'order=id' fails but 'order=line' passes, your data Atoms lines may not be written in atom-ID order.")
        print("  - If only H/C are swapped inside polymers, element-only check may miss subtle wrong mapping; prefer type/name-level checks when possible.")
        return 2

if __name__ == "__main__":
    raise SystemExit(main())
