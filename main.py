from atoms import Geometry, Matrix
import pandas as pd

basis_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\basis_ex.txt'
molecule_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\mol_ex.txt'

mol = Geometry(molecule_file, basis_file)

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

first_carbon = mol.atoms['C']['instance'][0]['GTO'].p(1, 0, 0.3).z()
print(first_carbon)
print('-'*100)

lap = Matrix(mol, 0)
df1 = pd.DataFrame(lap.matrix)
df1.to_excel('overlap_matrix.xlsx', index=False, header=False)

df2 = pd.DataFrame(lap.normalised_matrix)
df2.to_excel('normalised_overlap_matrix.xlsx', index=False, header=False)

kin = Matrix(mol, 1)
df3 = pd.DataFrame(kin.matrix)
df3.to_excel('kinetic_matrix.xlsx', index=False, header=False)

nuc = Matrix(mol, 2)
df4 = pd.DataFrame(nuc.matrix)
df4.to_excel('nuclear_attraction_matrix.xlsx', index=False, header=False)

rep = Matrix(mol, 3)
N = rep.matrix.shape[0]
flattened_matrix = rep.matrix.reshape(N*N, N*N)
df5 = pd.DataFrame(flattened_matrix)
df5.to_excel('repulsion_matrix_2D.xlsx', index=False, header=False)


print('end')