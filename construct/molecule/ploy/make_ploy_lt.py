# save as: make_poly_lt.py
import re

MOL2 = "ploy.mol2"                  # <-- 注意你的文件名
TYPES = "opls_types_by_index.txt"
OUT  = "poly.lt"
OPLS_LT = "oplsaa2024.lt"

USE_MOL2_COORDS = False  # True=把mol2坐标写入lt；False=坐标写0 0 0，靠 moltemplate -pdb 覆盖

def read_types(path):
    types = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.strip().startswith("#"):
                continue
            i, t = line.split()[:2]
            m = re.search(r"(\d+)", t)
            if not m:
                raise ValueError(f"Bad type: {t}")
            types[int(i)] = int(m.group(1))
    return types

def read_atoms_from_mol2(mol2_path):
    """Return charges{atom_id: q} and coords{atom_id:(x,y,z)} from @<TRIPOS>ATOM"""
    charges = {}
    coords  = {}
    in_atom = False
    with open(mol2_path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("@<TRIPOS>ATOM"):
                in_atom = True
                continue
            if s.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and s:
                parts = s.split()
                # mol2 atom line typical:
                # id name x y z type subst_id subst_name charge
                aid = int(parts[0])
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                q = float(parts[-1])  # 最后一列一般是charge
                charges[aid] = q
                coords[aid] = (x, y, z)
    if not charges:
        raise RuntimeError("No atoms parsed from mol2 @<TRIPOS>ATOM section.")
    return charges, coords

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
charges, coords = read_atoms_from_mol2(MOL2)
bonds = read_bonds_from_mol2(MOL2)

# 以 mol2 的原子数为准更安全
n_atoms = max(charges.keys())

# 基本一致性检查
missing_t = [i for i in range(1, n_atoms+1) if i not in types]
missing_q = [i for i in range(1, n_atoms+1) if i not in charges]
if missing_t:
    raise RuntimeError(f"Missing types for atom IDs: {missing_t[:10]} ... (total {len(missing_t)})")
if missing_q:
    raise RuntimeError(f"Missing charges for atom IDs: {missing_q[:10]} ... (total {len(missing_q)})")

qtot = sum(charges[i] for i in range(1, n_atoms+1))
print(f"Parsed atoms: {n_atoms}, total charge = {qtot:.6f}")

with open(OUT, "w") as out:
    out.write(f'import "{OPLS_LT}"\n\n')
    out.write("Poly inherits OPLSAA {\n")
    out.write('  write("Data Atoms") {\n')
    out.write("    # atomID   molID   atomType   charge   x y z\n")
    for i in range(1, n_atoms+1):
        q = charges[i]
        if USE_MOL2_COORDS:
            x, y, z = coords[i]
            out.write(f"    $atom:a{i:<4d}  $mol:.  @atom:{types[i]:<5d}  {q: .6f}   {x: .6f} {y: .6f} {z: .6f}\n")
        else:
            out.write(f"    $atom:a{i:<4d}  $mol:.  @atom:{types[i]:<5d}  {q: .6f}   0 0 0\n")
    out.write("  }\n\n")

    out.write('  write("Data Bond List") {\n')
    for k, (a1, a2) in enumerate(bonds, start=1):
        out.write(f"    $bond:b{k:<4d}  $atom:a{a1}  $atom:a{a2}\n")
    out.write("  }\n")
    out.write("}\n")

print("Wrote", OUT)
