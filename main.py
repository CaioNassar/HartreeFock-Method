from atoms import Geometry, Overlap
import pandas as pd

basis_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\basis_ex.txt'
molecule_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\mol_ex.txt'


mol = Geometry(molecule_file, basis_file)

elements = list(mol.atoms.keys())
print(elements)
print('-'*100)

print(mol.atoms['C']['basis'])
print('-'*100)

lap = Overlap(mol)
print(lap.matrix)
print('-'*100)

elements = list(mol.atoms.keys())
print(elements)
print('-'*100)

print(mol.atoms['C']['basis'])
print('-'*100)

carbons = mol.atoms['C']['instance']
print(carbons)
print('-'*100)

carbons_p = mol.atoms['C']['basis']['p']
print(carbons_p)
print('-'*100)

third_carbon = mol.atoms['C']['instance'][2]['GTO'].p(1, 0, 0.3).z()
print(third_carbon)
print('-'*100)

lap = Overlap(mol)
df = pd.DataFrame(lap.matrix)
df.to_excel('overlap_matrix.xlsx', index=False, header=False)


print()
