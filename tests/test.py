from ase import Atoms

# Example Atoms object
atoms = Atoms('H2O', positions=[(0, 0, 0), (0, 0.1, 0), (0, -0.1, 1)])
print(atoms.get_positions())

# Your permutation list
# As an example, let's say we want to swap the positions of the first and last atom.
# In this case, the permutation list would be [2, 1, 0].
permutation = [2, 1, 0]

# Reorder atoms based on your permutation
reordered_atoms = atoms[permutation]

print(reordered_atoms.get_positions())