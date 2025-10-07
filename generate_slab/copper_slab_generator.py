from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.io import write

# Parameters
element = 'Cu'
size = (4, 4, 3)  # (x, y, z) repetitions — this gives 4x4 atoms in-plane, 3 layers
vacuum = 10.0     # 10 Å vacuum
a = 3.615         # Lattice constant for Cu in Å (can be adjusted if needed)

# Build the slab
slab = fcc111(symbol=element, size=size, a=a, vacuum=vacuum, periodic=False)

# Output info
print(f"Total atoms: {len(slab)}")  # Should be close to 48

# Save to file (e.g., for VASP, XYZ, or other codes)
write('Cu111_slab.xyz', slab)
write('Cu111_slab.vasp', slab, format='vasp', direct=True)