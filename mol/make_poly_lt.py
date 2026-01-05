import re

MOL2 = "poly_ob.mol2"
TYPES = "opls_types_by_index.txt"
OUT = "poly.lt"
OPLS_LT = "oplsaa2024.lt"   

def read_types(path):
    types = {}
    with open(path) as f:
        for line in f:
            i, t = line.split()[:2]
            m = re.search(r"(\d+)", t)
            if not m:
                raise ValueError(f"Bad type: {t}")
            types[int(i)] = int(m.group(1))
    return types

def read_bonds_from_mol2(mol2_path):
    bonds = []
    in_bond = False
    with open(mol2_path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("@<TRIPOS>BOND"):
                in_bond = True
                continue
            if s.startswith("@<TRIPOS>") and in_bond:
                break
            if in_bond and s:
                parts = s.split()
                a1 = int(parts[1]); a2 = int(parts[2])
                bonds.append((a1, a2))
    return bonds

types = read_types(TYPES)
bonds = read_bonds_from_mol2(MOL2)
n_atoms = max(types.keys())

with open(OUT, "w") as out:
    out.write(f'import "{OPLS_LT}"\n\n')
    out.write("Poly inherits OPLSAA {\n")
    out.write('  write("Data Atoms") {\n')
    out.write("    # atomID   molID   atomType   charge   x y z\n")
    for i in range(1, n_atoms + 1):
        out.write(f"    $atom:a{i:<4d}  $mol:.  @atom:{types[i]:<5d}  0.0   0 0 0\n")
    out.write("  }\n\n")
    out.write('  write("Data Bond List") {\n')
    for k, (a1, a2) in enumerate(bonds, start=1):
        out.write(f"    $bond:b{k:<4d}  $atom:a{a1}  $atom:a{a2}\n")
    out.write("  }\n")
    out.write("}\n")

print("Wrote poly.lt")

