import json
import numpy as np
from atoms import Geometry
 
basis_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\basis_ex.txt'
xyz_path = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\xyz'

molecules = ['lif']

distances = np.linspace(0.45, 3.0, 30)
 
results = {}
for name in molecules:
    molecule_file = f'{xyz_path}\\{name}.txt'
    print(f'Calculando {name}...')
    mol = Geometry(molecule_file, basis_file)
    curve = mol.potential_energy_curve(distances)
    d_vals, e_vals = zip(*curve)
    results[name] = {
        'distances': [float(d) for d in d_vals],
        'energies': [float(e) for e in e_vals],
    }
 
output_path = 'atoms_curves.json'
with open(output_path, 'w') as f:
    json.dump(results, f, indent=2)

print(f'\nResultados salvos em {output_path}')