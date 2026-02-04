from rdkit import Chem
from rdkit.Chem import AllChem
from ase.io import read, write


mol = Chem.MolFromSmiles('C(C(C)C1=CC=C2C(=C1)C=CC=C2C(CCCCCCCCCCCC)C)CCCCCCCCCCC')  # dodecane
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
Chem.MolToMolFile(mol, '7_c.mol')
ase_mol = read('7_c.mol')
ase_mol.write('7_c.pdb')
