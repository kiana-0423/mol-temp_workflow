import parmed as pmd
from foyer import Forcefield

mol2 = "box.mol2"  
mol = pmd.load_file(mol2, structure=True)

ff = Forcefield(name="oplsaa")

typemap = ff.run_atomtyping(mol, use_residue_map=True)

with open("opls_types_by_index.txt", "w") as f:
    for atom in mol.atoms:
        t = typemap[atom.idx]["atomtype"]     
        f.write(f"{atom.idx + 1}\t{t}\n")

print("Wrote opls_types_by_index.txt")

