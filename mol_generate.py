from rdkit import Chem
from rdkit.Chem import AllChem
from ase.io import read, write


mol = Chem.MolFromSmiles('O')  # dodecane
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
Chem.MolToMolFile(mol, 'water.mol')
ase_mol = read('water.mol')
ase_mol.write('water.pdb')
