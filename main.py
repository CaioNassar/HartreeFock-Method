"""
main.py -- driver script: reads a molecule, runs the RHF calculation and prints
the energies and Mulliken charges.

Quick start
-----------
1. Put the geometry of the molecule in  xyz/<name>.txt  (xyz format, Angstrom).
2. Add <name> to the `molecules` list below.
3. Choose the initial density guess in INITIAL_GUESS.
4. Run:  python main.py

The commented-out blocks are optional tools (exporting matrices to
Excel, potential-energy curves). Uncomment a block together with the import
written at its top.
"""

import os
from time import time
from atoms import Geometry, Scf


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Initial density guess for the SCF:
#   'core'   -> diagonalise the core Hamiltonian H_core = T + V
#   'huckel' -> diagonalise an extended Huckel Hamiltonian (only for the
#               elements/shells that have parameters in Scf.Vsip)
#   'both'   -> run the SCF twice, once per guess, re-using the integrals, so
#               that energies and iteration counts can be compared
INITIAL_GUESS = 'huckel'

# Paths are relative to this file (main.py, basis_ex.txt and xyz/ are expected
# to live in the same folder).
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
BASIS_FILE = os.path.join(BASE_DIR, 'basis_ex.txt')
XYZ_DIR = os.path.join(BASE_DIR, 'xyz')

# Names of the files in xyz/ without ".txt"
molecules = ['ne+4']
# molecules = [f[:-4] for f in os.listdir(XYZ_DIR) if f.endswith('.txt')]
 
# Charged species reuse a neutral geometry file: name -> (xyz file, charge).
# charge > 0 removes electrons.
CHARGED = {
    'h2o6+': ('h2o', 6),
    'ne+4':  ('ne', 4),
    'ne+6':  ('ne', 6),
    'ne+8':  ('ne', 8),
}
 
guesses = ['core', 'huckel'] if INITIAL_GUESS.strip().lower() == 'both' else [INITIAL_GUESS]


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
for name in molecules:
    xyz_name, charge = CHARGED.get(name, (name, 0))
 
    t_0 = time()
    mol = Geometry(os.path.join(XYZ_DIR, f'{xyz_name}.txt'), BASIS_FILE, charge)
 
    # Integrals (S, T, V, two-electron) are computed once, even for 'both'.
    mol_scf = Scf(mol, initial_guess=guesses[0])
 
    energies = {}
    for n, guess in enumerate(guesses):
        if n > 0:
            mol_scf.set_initial_guess(guess)    # only D, J, K and F are rebuilt
        energies[guess] = mol_scf.scf()         # [one-electron, two-electron, nuclear rep., kinetic] in Eh
 
    print("Time in calculation: ", time() - t_0, "\n")
 
    # --- Optional: export matrices to Excel (pip install pandas openpyxl) ---
    # import pandas as pd
    # matrices = {'fock_matrix': mol_scf.F, 'normalised_overlap_matrix': mol_scf.S.matrix,
    #             'kinetic_matrix': mol_scf.T.matrix, 'nuclear_attraction_matrix': mol_scf.V.matrix,
    #             'hamiltonian_core_matrix': mol_scf.H_core, 'density_guess': mol_scf.D,
    #             'Coulomb_matrix': mol_scf.J, 'Exchange_matrix': mol_scf.K}
    # for fname, m in matrices.items():
    #     pd.DataFrame(m).to_excel(f'{fname}.xlsx', index=False, header=False)
    # N = mol_scf.Rep.matrix.shape[0]             # 4-index tensor -> 2-D
    # pd.DataFrame(mol_scf.Rep.matrix.reshape(N*N, N*N)).to_excel('repulsion_matrix_2D.xlsx', index=False, header=False)
 
    # --- Optional: potential-energy curve of a diatomic (pip install matplotlib) ---
    # import io, numpy as np, matplotlib.pyplot as plt
    # from contextlib import redirect_stdout
    # if mol.n_atoms == 2:
    #     with redirect_stdout(io.StringIO()):
    #         curve = mol.potential_energy_curve(np.linspace(0.45, 3, 25), initial_guess=guesses[0])
    #     fig, ax = plt.subplots()
    #     ax.plot(*zip(*curve), label=name)
    #     ax.set(xlabel='Internuclear Distance (Å)', ylabel='SCF Energy (Eh)')
    #     ax.legend()
    #     fig.savefig(os.path.join(BASE_DIR, f'scf_core_{name}.png'), dpi=200)
    #     plt.close(fig)
