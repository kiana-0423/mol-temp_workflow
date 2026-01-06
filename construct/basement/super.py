from ase.io import read, write
from ase.build import surface

atoms = read('Sio.cif')  
slab = surface(atoms, (1, 0, 0), layers=1, vacuum=0.0, periodic=True)
big = slab.repeat((41, 37, 6))
write('SiO2.cif',  big, format='cif')
write("output.pdb", big)