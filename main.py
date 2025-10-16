from atoms import Geometry

basis_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\basis_ex.txt'
molecule_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\mol_ex.txt'


mol = Geometry(molecule_file, basis_file)

elements = list(mol.atoms.keys())
print(elements)

print(mol.atoms['C']['basis'])


carbons = mol.atoms['C']['instance']
print(carbons)

carbons_p = mol.atoms['C']['basis']['p']
print(carbons_p)

third_carbon = mol.atoms['C']['instance'][2]['GTO'].p(1, 0, 0.3).z()
print(third_carbon)

print()
