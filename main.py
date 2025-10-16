from atoms import Geometry, Overlap

basis_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\basis_ex.txt'
molecule_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\mol_ex.txt'


mol = Geometry(molecule_file, basis_file)

<<<<<<< HEAD

elements = list(mol.atoms.keys())
print(elements)


lap = Overlap(mol)
print(lap.dimension)
=======
elements = list(mol.atoms.keys())
print(elements)

print(mol.atoms['C']['basis'])


carbons = mol.atoms['C']['instance']
print(carbons)

carbons_p = mol.atoms['C']['basis']['p']
print(carbons_p)

third_carbon = mol.atoms['C']['instance'][2]['GTO'].p(1, 0, 0.3).z()
print(third_carbon)
>>>>>>> 49f6dda29c2275a187dae06ff086eb6960a26057

print()
