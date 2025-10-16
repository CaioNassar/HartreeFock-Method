from atoms import Geometry, Overlap

basis_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\basis_ex.txt'
molecule_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\mol_ex.txt'


mol = Geometry(molecule_file, basis_file)


elements = list(mol.atoms.keys())
print(elements)


lap = Overlap(mol)
print(lap.dimension)

print()
