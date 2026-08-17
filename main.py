import matplotlib.pyplot as plt
from atoms import Geometry, Matrix, Scf
import pandas as pd
import numpy as np
from time import time
from contextlib import redirect_stdout
import io
from line_profiler import LineProfiler
import os

basis_file = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\basis_ex.txt'
xyz_path = 'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\xyz'

molecules = [arquivo.replace('.txt', '') for arquivo in os.listdir(xyz_path) if arquivo.endswith('.txt')]


for geo in molecules:
    if geo == 'br2':
        continue

    if geo == 'h2o6+':
        molecule_file = f'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\xyz\\h2o.txt'
        t_0 = time()
        mol = Geometry(molecule_file, basis_file, 6)
    elif geo == 'ne+8':
        molecule_file = f'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\xyz\\ne.txt'
        t_0 = time()
        mol = Geometry(molecule_file, basis_file, 8)
    elif geo == 'ne+6':
        molecule_file = f'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\xyz\\ne.txt'
        t_0 = time()
        mol = Geometry(molecule_file, basis_file, 6)
    elif geo == 'ne+4':
        molecule_file = f'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\xyz\\ne.txt'
        t_0 = time()
        mol = Geometry(molecule_file, basis_file, 4)
    else:
        molecule_file = f'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\xyz\\{geo}.txt'
        t_0 = time()
        mol = Geometry(molecule_file, basis_file)



    mol_scf = Scf(mol)

    # elements = list(mol.atoms.keys())
    # print(elements)
    # print('-'*100)

    # print(mol.atoms['C']['basis'])
    # print('-'*100)

    # carbons = mol.atoms['C']['instance']
    # print(carbons)
    # print('-'*100)

    # carbons_p = mol.atoms['C']['basis']['p']
    # print(carbons_p)
    # print('-'*100)

    # first_carbon = mol.atoms['C']['instance'][0]['GTO'].p(1, 0, 0.3).z()
    # print(first_carbon)
    # print('-'*100)

    # df1 = pd.DataFrame(mol_scf.F)
    # df1.to_excel('fock_matrix.xlsx', index=False, header=False)

    # df2 = pd.DataFrame(mol_scf.S.matrix)
    # df2.to_excel('normalised_overlap_matrix.xlsx', index=False, header=False)

    # df3 = pd.DataFrame(mol_scf.T.matrix)
    # df3.to_excel('kinetic_matrix.xlsx', index=False, header=False)

    # df4 = pd.DataFrame(mol_scf.V.matrix)
    # df4.to_excel('nuclear_attraction_matrix.xlsx', index=False, header=False)

    # rep = mol_scf.Rep.matrix
    # N = rep.shape[0]
    # flattened_matrix = rep.reshape(N*N, N*N)
    # df5 = pd.DataFrame(flattened_matrix)
    # df5.to_excel('repulsion_matrix_2D.xlsx', index=False, header=False)

    # df6 = pd.DataFrame(mol_scf.H_core)
    # df6.to_excel('hamiltonian_core_matrix.xlsx', index=False, header=False)

    # df7 = pd.DataFrame(mol_scf.H_huckel)
    # df7.to_excel('hamiltonian_huckel_matrix.xlsx', index=False, header=False)

    # df8 = pd.DataFrame(mol_scf.D)
    # df8.to_excel('density_guess.xlsx', index=False, header=False)

    # df9 = pd.DataFrame(mol_scf.J)
    # df9.to_excel('Coulomb_matrix.xlsx', index=False, header=False)

    # dfA = pd.DataFrame(mol_scf.K)
    # dfA.to_excel('Exchange_matrix.xlsx', index=False, header=False)

    energy = mol_scf.scf()
    t_f = time()
    print("Time in calculation: ", (t_f - t_0), "\n")

    # if mol.n_atoms != 2:
    #       continue
    # scf_dist = np.linspace(0.45, 3, 25)
    # with redirect_stdout(io.StringIO()):
    #       curve = mol.potential_energy_curve(scf_dist)

    # fig, ax = plt.subplots()
    # x_vals, y_vals = zip(*curve)
    # ax.plot(x_vals, y_vals, label=geo)

    # ax.set_xlabel('Internuclear Distance (Å)')
    # ax.set_ylabel('SCF Energy (Eh)')
    # ax.legend()

    # fig.savefig(f'C:\\Users\\Gerenciador\\Documents\\programming\\HartreeFock-Method\\scf_core_{geo}.png', dpi=200)
    # plt.close(fig)
