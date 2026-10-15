# HartreeFock-Method

A from-scratch **restricted Hartree–Fock (RHF)** code for atoms and molecules using **contracted Cartesian-Gaussian basis functions**. It is written in Python, with the most expensive step — the two-electron repulsion integrals — implemented in C and called through `ctypes`. Results are validated against [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas).

Developed at the Institute of Physics of the University of São Paulo (IFUSP) as part of an undergraduate research project (PUB-USP).

## What it does

- Reads a geometry (`xyz` file) and a basis set (Molcas format).
- Computes the overlap (**S**), kinetic (**T**), nuclear-attraction (**V**) and two-electron repulsion $(\mu\nu|\lambda\sigma)$ integrals analytically.
- Solves the Roothaan equations $\mathbf{FC}=\mathbf{SC}\boldsymbol\varepsilon$ with a self-consistent-field (SCF) loop accelerated by Pulay's **DIIS**.
- Offers two initial guesses for the density: the **core Hamiltonian** or **extended Hückel**.
- Prints total, one-electron, two-electron, nuclear-repulsion and kinetic energies, and **Mulliken charges**.
- Scans potential-energy curves of diatomic molecules.

## Theory in short

With $\mathbf D = \mathbf C_\text{occ}\mathbf C_\text{occ}^{T}$ the density matrix of one spin:

$$
F_{\mu\nu}=H_{\mu\nu}+\sum_{\lambda\sigma}D_{\lambda\sigma}\big[2(\mu\nu|\lambda\sigma)-(\mu\sigma|\lambda\nu)\big]\qquad(\mathbf F=\mathbf H+2\mathbf J-\mathbf K)
$$

$$
E_\text{RHF}=E_\text{Nuc}+\sum_{\mu\nu}D_{\mu\nu}\big(H_{\mu\nu}+F_{\mu\nu}\big)
$$

SCF loop:

1. **Input** – geometry and basis set.
2. **Integrals** – $\mathbf S,\mathbf T,\mathbf V$ and $(\mu\nu|\lambda\sigma)$ (the latter in C).
3. **Initial guess** – density $\mathbf D$ from the core Hamiltonian or from extended Hückel.
4. **Fock matrix** – $\mathbf F=\mathbf H+2\mathbf J-\mathbf K$ built from the current $\mathbf D$.
5. **DIIS update** – extrapolate $\mathbf F$ (last 7 vectors), diagonalise, get new $\mathbf C,\mathbf D$.
6. **Convergence test** – stop when $\lVert\mathbf G\rVert_\infty<10^{-5}$, where $\mathbf G=\mathbf X^{\dagger}(\mathbf{FDS}-\mathbf{SDF})\mathbf X$ is the DIIS orbital gradient; otherwise go back to step 4.

Because the Gaussian basis is not orthonormal, the equations are solved by **canonical orthogonalisation**: $\mathbf X=\mathbf U\,\mathbf s^{-1/2}$ with $\mathbf S=\mathbf U\mathbf s\mathbf U^{T}$, discarding eigenvalues of $\mathbf S$ below $10^{-5}$ (so that $\mathbf X^{\dagger}\mathbf S\mathbf X=\mathbf 1$).

The integrals follow the Taketa–Huzinaga–O-ohata scheme as presented in Cook, *Handbook of Computational Chemistry* (1998), with the Boys function evaluated by a power series for small arguments and by the error function for large ones.

## Repository layout

| File             | Role                                                                                                                                                          |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `main.py`      | Driver script: choose the molecule(s) and the initial guess, run the SCF.                                                                                     |
| `atoms.py`     | `Geometry` (input files), `Matrix` (S, T, V and repulsion tensor), `Scf` (SCF loop, DIIS, guesses, Mulliken), plus helpers that evaluate GTOs in space. |
| `integrals.c`  | C routines for the integrals.                                                                                                  |
| `basis_ex.txt` | Basis set file (Molcas format).                                                                                                                               |
| `xyz/`         | Molecular geometries, one`<name>.txt` per molecule.                                                                                                         |

## Requirements

- Python 3 (3.12 is known to work)
- `numpy`
- a C compiler with C99 support — **gcc** or **clang** (MSVC is not supported: the C code uses variable-length arrays)

```bash
pip install numpy
```

Optional, only if you uncomment the corresponding blocks in `main.py`: `pandas` + `openpyxl` (export matrices to Excel), `matplotlib` (plot curves).

## Building the C library

```bash
# Windows (MinGW-w64)
gcc -O3 -shared -o integrals.dll integrals.c

# Linux
gcc -O3 -shared -fPIC -o integrals.so integrals.c -lm

# macOS
gcc -O3 -shared -fPIC -o integrals.dylib integrals.c -lm
```

## Input files

### Geometry (`xyz/<name>.txt`)

Standard xyz format, coordinates in **Å**:

```
2
comment line (ignored)
Li 0.0 0.0 0.0
H  0.0 0.0 1.595
```

The number of electrons is the sum of the atomic numbers minus the charge. It must be **even** (closed shell).

### Basis set (`basis_ex.txt`, Molcas format)

The parser follows these rules:

- A line containing both `/` and `inline`, whose **first token is the element symbol** (e.g. `H / inline`), starts that element's block. Blocks of elements that are not in the molecule are ignored.
- A line containing `-type` (e.g. `* S-type functions`) starts a shell. The shell label (`S`, `P`, `D`, `F`) must be the **second** token of the line.
- The next line gives `<number of primitives> <number of contracted functions>`.
- Then come one exponent per line (one line per primitive), followed by one line of contraction coefficients per primitive (one column per contracted function).
- A line containing the word `cartesian` marks the set as Cartesian (this only affects the helper that evaluates d orbitals in space, `Basis.d`).

Every element in the molecule must have a block in the file.

## Running

Edit the top of `main.py`:

```python
INITIAL_GUESS = 'core'      # 'core', 'huckel' or 'both'
molecules = ['lif']         # names of the files in xyz/ (without .txt)
```

Charged species reuse a neutral geometry file. Register them in the `CHARGED`
dictionary in `main.py` (`name -> (xyz file, charge)`):

```python
CHARGED = {'ne+4': ('ne', 4), 'h2o6+': ('h2o', 6)}
molecules = ['ne+4']
```

and run:

```bash
python main.py
```

The paths to `basis_ex.txt` and `xyz/` are relative to `main.py`.

### Choosing the initial guess

| `INITIAL_GUESS`    | Behaviour                                                                                                                                                                                                                                                                                                                                                   |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `'core'` | Density from the diagonalisation of$\mathbf H_\text{core}=\mathbf T+\mathbf V$. Works for any element.                                                                                                                                                                                                                                                    |
| `'huckel'`         | Density from an extended Hückel Hamiltonian:$H_{ii}=\text{VSIP}$, $H_{ij}=1.75\,S_{ij}(H_{ii}+H_{jj})/2$. Parameters (`Scf.Vsip`) exist only for H and He (s) and for Li, B, C, N, O, F, Ne, Br and Kr (s and p); the Kr entries are zero placeholders. With any other element or shell the code stops with a clear error — use `'core'` instead. |
| `'both'`           | Runs the SCF twice (core, then Hückel) re-using the integrals, so the two guesses can be compared.                                                                                                                                                                                                                                                         |

Both guesses normally converge to the same energy; they typically differ in the number of iterations. The summary printed at the end of each SCF states which guess was used.

### Using the code from Python

```python
from atoms import Geometry, Scf

mol = Geometry('xyz/lih.txt', 'basis_ex.txt')           # optional 3rd argument: charge
scf = Scf(mol, initial_guess='core')                    # computes the integrals once

one_e, two_e, e_nuc, kinetic = scf.scf()                # [Eh]; total = sum of the first three
print('E_total =', one_e + two_e + e_nuc)

# Re-run from another starting density without recomputing the integrals:
scf.set_initial_guess('huckel')
scf.scf()
```

`scf.scf()` returns the integer `404` if it does not converge within `max_iter` iterations (default 100). The convergence threshold on $\lVert\mathbf G\rVert_\infty$ is its first argument (`max_cond_number`, default `1e-5`).

### Potential-energy curves (diatomic molecules)

```python
import numpy as np
from atoms import Geometry

mol = Geometry('xyz/lih.txt', 'basis_ex.txt')
curve = mol.potential_energy_curve(np.linspace(0.8, 4.0, 25), initial_guess='core')   # distances in Å
# curve = [(distance, total_energy_in_Eh), ...]
```

## Output

For each molecule the SCF prints:

- **Total SCF energy** (Eh) = one-electron + two-electron + nuclear repulsion
- **One-electron energy** $=2\sum D_{\mu\nu}H^\text{core}_{\mu\nu}$ and **two-electron energy** $=\sum D_{\mu\nu}(2J-K)_{\mu\nu}$
- **Nuclear repulsion** and **kinetic energy**
- the number of SCF iterations
- **Mulliken charges**: gross population of each basis function, and total electronic and net charge of each atom

## Validation

The energies and the potential-energy curves of H₂, HF, He₂, Li₂, B₂, Br₂, C₂, N₂, O₂, F₂, LiH, LiF and CO₂ (all with 6-31G) were compared with OpenMolcas.

## Ideas for future work

- Add molecular-symmetry checks (point groups).
- Slater-type orbitals (STOs) in the integral engine.
- Porting more of the performance-critical code to C and, eventually, to GPU.

## References

1. A. Szabo and N. S. Ostlund, *Modern Quantum Chemistry*, Dover (1996).
2. D. B. Cook, *Handbook of Computational Chemistry*, Oxford University Press (1998).
3. P. Pulay, *Convergence acceleration of iterative sequences. The case of SCF iteration*, Chem. Phys. Lett. 73, 393 (1980).
4. I. Fdez. Galván *et al.*, *OpenMolcas: From Source Code to Insight*, J. Chem. Theory Comput. 15, 5925 (2019).
