import parmed as pmd
from foyer import Forcefield

mol2 = "ploy.mol2"
mol = pmd.load_file(mol2, structure=True)

ff = Forcefield(name="oplsaa")
typemap = ff.run_atomtyping(mol, use_residue_map=True)

with open("opls_types_by_index.txt", "w") as f:
    for atom in mol.atoms:
        f.write(f"{atom.idx + 1}\t{typemap[atom.idx]['atomtype']}\n")
