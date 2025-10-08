from atoms import Geometry

basis_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\basis_ex.txt'
molecule_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\mol_ex.txt'


mol = Geometry(molecule_file, basis_file)

elements = list(mol.atom.keys())
print(elements)

hydrogen_list = mol.atom['H']
print(hydrogen_list)
print(len(hydrogen_list))


third_carbon = mol.atom['C'][2]
print(third_carbon)


C_s = third_carbon['Data'].s(-1, 0, -0.753)
print(C_s)

print()
