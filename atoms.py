import numpy as np
from math import comb, factorial
import scipy.special as sp
from ntpath import basename, splitext
from tempfile import TemporaryDirectory
from pathlib import Path
import ctypes
import os

base_dir = os.path.dirname(os.path.abspath(__file__))
lib_path = os.path.join(base_dir, "integrals.dll")

c_lib = ctypes.CDLL(lib_path)

c_lib.calculate_B.argtypes = [
    ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,     
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, 
    ctypes.c_double, ctypes.c_double,                           
    ctypes.c_double, ctypes.c_double, ctypes.c_double,          
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
]
c_lib.calculate_B.restype = None

c_lib.calculate_A.argtypes = [
    ctypes.c_int, ctypes.c_int,                                 
    ctypes.c_double, ctypes.c_double, ctypes.c_double,          
    ctypes.c_double, ctypes.c_double,                           
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
]
c_lib.calculate_A.restype = None

c_lib.boys_recurrence.argtypes = [
    ctypes.c_int, ctypes.c_double,
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
]
c_lib.boys_recurrence.restype = None

c_lib.S_integral.argtypes = [
    ctypes.c_int, ctypes.c_int,                  
    ctypes.c_double, ctypes.c_double, ctypes.c_double 
]
c_lib.S_integral.restype = ctypes.c_double

c_lib.repulsion_integral4.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),    # lmn
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),  # A
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),  # B
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),  # C
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),  # D
    ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,                  # nums
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),  # exp1
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),  # coef1
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),  # exp2
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),  # coef2
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),  # exp3
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),  # coef3
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),  # exp4
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')   # coef4
]
c_lib.repulsion_integral4.restype = ctypes.c_double 

c_lib.normalize_N.argtypes = [
    ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int
]
c_lib.normalize_N.restype = ctypes.c_double

class Geometry:
    """
    The class Geometry stores the 3D coordinates, atomic composition and basis set of the atoms of the molecule.

    An instance of this class represents a molecule, reading its geometry.
    """
    atomic_number = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
        'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
        'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
        'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40,
        'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
        'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60,
        'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70,
        'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80,
        'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90,
        'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100,
        'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109,
        'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
        # This class attribute maps chemical symbols to their respective atomic numbers.
    }

    

    def __init__(self, path_coord, path_basis, charge=0):
        """
        Initializes the Geometry object.

        Args:
            path_coord (str): The file path to the atomic coordinates file.
            path_basis (str): The file path to the basis set.
        """
        self.path_coord = path_coord
        self.path_basis = path_basis
        self.atoms = {}  # A dictionary to cache the data read from the files.
        self.atom_list = []
        self.n_atoms = 0 # Number of atoms in the element.
        self.n_electrons = 0
        self.charge = charge
        self.cartesian = 0
        self.ang2bohr = 1.88972612590778

        # Reads the file of the atomic coordinates, sets self.atoms.
        first_file = self._read_coordinates()
        
        # Reads the file of the basis set and stores the data of the atom in the gaussians dictionary.
        if first_file:    
            second_file = self._read_basis_set()

        n_electrons = 0
        if first_file and second_file:
            # Loop through each element symbol and its atoms.
            for symbol, element_list in self.atoms.items():
                # Loop through each atom in the element.
                n_electrons += self.atomic_number[symbol]*len(element_list['instance'])
                for atom_instance in element_list['instance']:
                    # Create a Basis object for this atom and add it to the dictionary.
                    atom_instance['GTO'] = Basis(atom_instance, element_list['basis'], self.cartesian)

            self.n_electrons = n_electrons - charge

    def _read_coordinates(self):
        """
        Reads the file and stores the center coordinates of each atom.

        Return:
            bool: True if the file was read successfully, False otherwise.
        """
        
        try:
            with open(self.path_coord, 'r') as file:
                lines = list(file)
                self.n_atoms = int(lines[0])

                atomic_dict = {}
                for i in range(self.n_atoms):
                    l = lines[i+2].strip().split()
                    symbol = l[0]
                    atom_coord = {
                        'symbol': symbol,
                        'Ax': float(l[1]) * self.ang2bohr,
                        'Ay': float(l[2]) * self.ang2bohr,
                        'Az': float(l[3]) * self.ang2bohr
                    }

                    self.atom_list.append(atom_coord)
                    # Stores the data for the current atom in the dictionary.
                    if symbol not in atomic_dict:
                        atomic_dict[symbol] = {
                            'basis' : {},
                            'instance' : [atom_coord]
                        }
                    else:
                        atomic_dict[symbol]['instance'].append(atom_coord)

                ordered_atoms = sorted(atomic_dict.keys(), key=lambda item: self.atomic_number.get(item), reverse=True)
                for element in ordered_atoms:
                    self.atoms[element] = atomic_dict[element]

            return True
            
        except Exception as e:
            # Exception if happens some error.
            print(f'Something went wrong when trying to read the file: {e}')
        return False
    
    def _read_basis_set(self):
        """
        Read the file and stores info for each different element in the molecule.

        Return:
            bool: True if the file was read successfully, False otherwise.
        """

        current_element = None
        try:
            with open(self.path_basis, 'r') as file:
                lines_iterator = iter(file)
                for line in lines_iterator:
                    line = line.strip()
                    if 'cartesian' in line:
                        self.cartesian = 1
                    if '/' in line and 'inline' in line:
                        # Search the flag(/inline) to read the basis set of an atom.
                        current_element = line.split()[0]
                        # Verifies if the chemical symbol is of interest.
                        if current_element not in self.atoms:
                            current_element = None
                        continue

                    if current_element and '-type' in line.lower():
                        # Searchs the line that specifies the orbital type (e.g. 'S-type', 'p-type', etc.).
                        # For each orbital, creates an key of atoms dictionary
                        exp = []
                        coef = []

                        orbital = line.split()[1].lower().replace('-type', '')
                        
                        functions_line = next(lines_iterator).strip().split()
                        num_primitives = int(functions_line[0])
                        num_contracted = int(functions_line[1])

                        # Parse exponents and coefficients.
                        for _ in range(num_primitives):
                            exp.append(float(next(lines_iterator).strip()))

                        for _ in range(num_primitives):
                            coef.append(next(lines_iterator).strip().split())
                        
                        self.atoms[current_element]['basis'][orbital] = {
                            'number_of_primitives': num_primitives,
                            'contracted_functions': num_contracted,
                            'exponents': exp,
                            'coefficients': coef,
                        } # Dictionary that stores the parsed data.
            

            
            for element in list(self.atoms.keys()):
                if not self.atoms[element]['basis']:
                    raise ValueError(f'Configuration ERROR: The element "{element}" does NOT have the basis included in the file.')
                
            return True
    
        except Exception as e:  
            # Exception if happens some error.
            print(f'Something went wrong when trying to read the file: {e}')
            return False
            
    def potential_energy_curve(self, distances):
        if self.n_atoms != 2:
            raise ValueError('This function only works with 2 atoms')

        symbols = [atom['symbol'] for atom in self.atom_list]
        curve = []
        for distance in distances:
            z = distance/2
            xyz = f"2\n\n{symbols[0]} 0.0 0.0 {-z}\n{symbols[1]} 0.0 0.0 {z}\n"
            with TemporaryDirectory() as tmp_dir:
                tmp_path = Path(tmp_dir) / 'geometry.txt'
                tmp_path.write_text(xyz)
                mol = Geometry(str(tmp_path), self.path_basis, self.charge)
                one_elec, two_elec, e_nuc, _ = Scf(mol).scf()
            curve.append((distance, one_elec + two_elec + e_nuc))

        return curve

class Matrix:
    Quantum_numbers = {
    # l = 0
    's': [(0, 0, 0)],  # s

    # l = 1
    'p': [(1, 0, 0),  # px
          (0, 1, 0),  # py
          (0, 0, 1)], # pz

    # l = 2
    'd': [(2, 0, 0),  # dxx
          (0, 2, 0),  # dyy
          (0, 0, 2),  # dzz
          (1, 1, 0),  # dxy
          (1, 0, 1),  # dxz
          (0, 1, 1)], # dyz

    # l = 3
    'f': [(3, 0, 0),  # fxxx
          (0, 3, 0),  # fyyy
          (0, 0, 3),  # fzzz
          (2, 1, 0),  # fxxy
          (2, 0, 1),  # fxxz
          (1, 2, 0),  # fxyy
          (0, 2, 1),  # fyyz
          (1, 0, 2),  # fxzz
          (0, 1, 2),  # fyzz
          (1, 1, 1)]  # fxyz
    }

    def __init__(self, geometry, type_matrix):
        """
        Args:
        type_matrix: 0: Overlap, 1: Kinetic, 2:  Nuclear_Attraction, 3: Repulsion
        """
        self.geometry = geometry
        self.type_matrix = type_matrix

        basis_functions = []
        molecule = []

        for symbol, data in self.geometry.atoms.items():
            for center in data['instance']:
                molecule.append({
                    'symbol': symbol,
                    'center': center,
                    'basis': data['basis']
                })
                
        for atom in molecule:
            for orbital, data in atom['basis'].items():
                quant_numbers = self.Quantum_numbers.get(orbital, [])
                for k in range(data['contracted_functions']):
                    for (l, m, n) in quant_numbers:
                        function = {
                            'symbol': atom['symbol'],
                            'center': atom['center'],
                            'l': l,
                            'm': m,
                            'n': n,
                            'primitives': [],
                        }
                        for i in range(data['number_of_primitives']):
                            coef = float(data['coefficients'][i][k]) 
                            if coef != 0.0: 
                                exponent = data['exponents'][i]
                                function['primitives'].append({
                                    'exponent': exponent,
                                    'coefficient': coef * c_lib.normalize_N(exponent, l, m, n)
                                })
                        
                        basis_functions.append(function)
        self.basis_functions = basis_functions
        dimension = len(basis_functions)

        normalize = np.zeros(dimension)
        for i in range(dimension):
            normalize[i] = 1/np.sqrt(self._integral(basis_functions[i], basis_functions[i], 0))
        
        if type_matrix == 3:
            self.matrix = np.zeros((dimension, dimension, dimension, dimension))

            for i in range(dimension):
                for j in range(i, dimension):
                    for k in range(dimension):
                        for l in range(k, dimension):
                            if (k, l) >= (i, j):
                                integral  = self._integral4(basis_functions[i], basis_functions[j], basis_functions[k], basis_functions[l])

                                value = integral * normalize[i] * normalize[j] * normalize[k] * normalize[l]
                                self.matrix[i, j, k, l] = value
                                self.matrix[i, j, l, k] = value
                                self.matrix[j, i, k, l] = value
                                self.matrix[j, i, l, k] = value
                                self.matrix[k, l, i, j] = value
                                self.matrix[k, l, j, i] = value
                                self.matrix[l, k, i, j] = value
                                self.matrix[l, k, j, i] = value
                             

        else:
            self.matrix = np.zeros((dimension, dimension))

            for i in range(dimension):
                for j in range(i, dimension):
                        if type_matrix == 0 or type_matrix == 1:
                            integral = self._integral(basis_functions[i], basis_functions[j], type_matrix)

                        elif type_matrix == 2:
                            integral = 0
                            for atom in molecule:
                                integral += self._integral(basis_functions[i], basis_functions[j], type_matrix, atom)

                        value = integral * normalize[i] * normalize[j]
                        self.matrix[i, j] = value
                        self.matrix[j, i] = value
            

    def _integral4(self, eta1, eta2, eta3, eta4):
        A = np.array([eta1['center']['Ax'], eta1['center']['Ay'], eta1['center']['Az']], dtype=np.float64)
        B = np.array([eta2['center']['Ax'], eta2['center']['Ay'], eta2['center']['Az']], dtype=np.float64)
        C = np.array([eta3['center']['Ax'], eta3['center']['Ay'], eta3['center']['Az']], dtype=np.float64)
        D = np.array([eta4['center']['Ax'], eta4['center']['Ay'], eta4['center']['Az']], dtype=np.float64)

        lmn = np.array([
            eta1['l'], eta1['m'], eta1['n'],
            eta2['l'], eta2['m'], eta2['n'],
            eta3['l'], eta3['m'], eta3['n'],
            eta4['l'], eta4['m'], eta4['n']
        ], dtype=np.int32)

        e1 = np.array([p['exponent'] for p in eta1['primitives']], dtype=np.float64)
        c1 = np.array([p['coefficient'] for p in eta1['primitives']], dtype=np.float64)
        e2 = np.array([p['exponent'] for p in eta2['primitives']], dtype=np.float64)
        c2 = np.array([p['coefficient'] for p in eta2['primitives']], dtype=np.float64)
        e3 = np.array([p['exponent'] for p in eta3['primitives']], dtype=np.float64)
        c3 = np.array([p['coefficient'] for p in eta3['primitives']], dtype=np.float64)
        e4 = np.array([p['exponent'] for p in eta4['primitives']], dtype=np.float64)
        c4 = np.array([p['coefficient'] for p in eta4['primitives']], dtype=np.float64)

        value = c_lib.repulsion_integral4(
            lmn, A, B, C, D,
            len(e1), len(e2), len(e3), len(e4),
            e1, c1, e2, c2, e3, c3, e4, c4
        )
        
        return value

    def _integral(self, eta1, eta2, type_matrix = 0, nucleus = 0):
        value = 0.0
        A = eta1['center']
        B = eta2['center']

        ab2 = (A['Ax'] - B['Ax'])**2 + (A['Ay'] - B['Ay'])**2 + (A['Az'] - B['Az'])**2 
        
        for primitive1 in eta1['primitives']:
            for primitive2 in eta2['primitives']:
                C_a = primitive1['coefficient']
                C_b = primitive2['coefficient']

                integral = 0.0

                if type_matrix == 0:
                    integral += self._overlap(eta1, eta2, primitive1, primitive2, ab2)
            
                if type_matrix == 1:
                    integral += self._kinetic_energy(eta1, eta2, primitive1, primitive2, ab2)    

                if type_matrix == 2:
                    integral += self._nuclear_attraction(eta1, eta2, primitive1, primitive2, nucleus, ab2)

                value += C_a*C_b*integral

        return value
    
    def _overlap(self, eta1, eta2, p1, p2, ab2):
        A = eta1['center']
        B = eta2['center']

        alpha = p1['exponent']
        beta = p2['exponent']
        gamma = alpha + beta

        Px = (alpha*A['Ax'] + beta*B['Ax'])/gamma
        Py = (alpha*A['Ay'] + beta*B['Ay'])/gamma
        Pz = (alpha*A['Az'] + beta*B['Az'])/gamma
        
        Sx = c_lib.S_integral(eta1['l'], eta2['l'], Px - A['Ax'], Px - B['Ax'], gamma)
        Sy = c_lib.S_integral(eta1['m'], eta2['m'], Py - A['Ay'], Py - B['Ay'], gamma)
        Sz = c_lib.S_integral(eta1['n'], eta2['n'], Pz - A['Az'], Pz - B['Az'], gamma)    

        integral = np.exp(-alpha*beta*ab2/gamma)*Sx*Sy*Sz

        return integral
    
    def _kinetic_energy(self, eta1, eta2, p1, p2, ab2):
        alpha = p1['exponent']
        value = alpha*(2*(eta1['l'] + eta1['m'] + eta1['n']) + 3)*self._overlap(eta1, eta2, p1, p2, ab2)
        coord = ['l', 'm', 'n']
        
        eta1_c = eta1.copy()
        for c in coord:
            eta1_c[c] += 2
            value += -2*(alpha**2)*self._overlap(eta1_c, eta2, p1, p2, ab2)
            eta1_c[c] -= 2

        
        for c in coord:
            if eta1_c[c] >= 2:
                eta1_c[c] -= 2
                value += -0.5*((eta1_c[c] + 2)*(eta1_c[c] + 1))*self._overlap(eta1_c, eta2, p1, p2, ab2)
                eta1_c[c] += 2

        return value
        
    def _nuclear_attraction(self, eta1, eta2, p1, p2, nucleus, ab2):
        A = eta1['center']
        B = eta2['center']
        C = nucleus['center']

        alpha = p1['exponent']
        beta = p2['exponent']
        gamma = alpha + beta
        eps = 1/(4*gamma)

        P = {
            'Px' : (alpha*A['Ax'] + beta*B['Ax'])/gamma,
            'Py' : (alpha*A['Ay'] + beta*B['Ay'])/gamma,
            'Pz' : (alpha*A['Az'] + beta*B['Az'])/gamma
        }
        
        pc2 = (P['Px'] - C['Ax'])**2 + (P['Py'] - C['Ay'])**2 + (P['Pz'] - C['Az'])**2

        l1, l2 = eta1['l'], eta2['l']
        A_x = np.zeros(l1 + l2 + 1, dtype=np.float64)
        c_lib.calculate_A(l1, l2, A['Ax'], B['Ax'], C['Ax'], P['Px'], eps, A_x)

        m1, m2 = eta1['m'], eta2['m']
        A_y = np.zeros(m1 + m2 + 1, dtype=np.float64)
        c_lib.calculate_A(m1, m2, A['Ay'], B['Ay'], C['Ay'], P['Py'], eps, A_y)

        n1, n2 = eta1['n'], eta2['n']
        A_z = np.zeros(n1 + n2 + 1, dtype=np.float64)
        c_lib.calculate_A(n1, n2, A['Az'], B['Az'], C['Az'], P['Pz'], eps, A_z)

        v_max = l1 + m1 + n1 + l2 + m2 + n2
        boys = np.zeros(v_max + 1, dtype=np.float64)
        c_lib.boys_recurrence(v_max, gamma*pc2, boys)

        K = 2*np.pi*np.exp(-alpha*beta*ab2/gamma)/gamma

        value = 0.0
        for x in range(len(A_x)):
            for y in range(len(A_y)):
                for z in range(len(A_z)):
                    value += A_x[x]*A_y[y]*A_z[z]*boys[x+y+z]

        return -self.geometry.atomic_number[nucleus['symbol']]*K*value

class Basis:
    """
    The class Basis calculates Cartesians-Gaussians-Type Orbitals (GTOs).

    An instance of this class represents a single atom of the molecule.
    """

    def __init__(self, center, basis_data, cartesian_basis):
        """
        Initializes the Atom object.

        Args:
            center (dict): center coordinates of the atom.
            basis_data (dict): info about the basis set of the element.
        """

        self.Ax = center['Ax']
        self.Ay = center['Ay']
        self.Az = center['Az']
        self.basis = basis_data
        self.cartesian = cartesian_basis

    def _calculate(self, x, y, z, sel, orbital):
        """
        Calculates the value of the GTOs.

        Args:
            x, y, z (float): coordinates of evaluating point.
            sel (int): selects the contraction option.
            orbital (str): the orbital itself (s, p, d or f).
        """

        data = self.basis[orbital]
        exp = data['exponents']
        coef = []

        if 0 <= sel < data['contracted_functions']:
            # Validates the selection.      
            for c in data['coefficients']:
                coef.append(c[sel])
        else:
            raise ValueError(
                f'This selection does not exist, choose from 0 to {data['contracted_functions'] - 1}')
            # Raises if selected option is out of bounds.
        
        contracted_value = 0.0
        r2 = (self.Ax - x)**2 + (self.Ay - y)**2 + (self.Az - z)**2
        # Calculate the squared distance from the atom's center.
        
        for i in range(data['number_of_primitives']):
            if coef[i] != 0.0:
                # Skips calculation if the coefficient is zero.
                contracted_value += float(coef[i]) * \
                    np.exp(-float(exp[i]) * r2)
                # Sum of each primitive: C * e^(-a * r^2).
        return contracted_value

    def s(self, x, y, z, sel=0):
        """
        Returns the value for s-orbital.
        
        Args:
            x, y, z (float): coordinates of evaluating point.
            sel (int): selects the contraction option.
        """

        return self._calculate(x, y, z, sel, 's')

    def p(self, x, y, z, sel=0):
        """
        Returns the value for p-orbitals.
        The _calculate function returns the radial part, so it's necessary to specify the axis of projection.
        ResultP object calculates for the determined component.

        Args:
            x, y, z (float): coordinates of evaluating point.
            sel (int): selects the contractiion option.
        """
        
        radial = self._calculate(x, y, z, sel, 'p')
        return ResultP(radial, x, y, z, self.Ax, self.Ay, self.Az)
    
    def d(self, x, y, z, sel=0):
            """
            Returns the value for d-orbitals.
            ResultD object calculates for the determined component.
    
            Args:
                x, y, z (float): coordinates of evaluating point.
                sel (int): selects the contractiion option.
            """

            radial = self._calculate(x, y, z, sel, 'd')

            if self.cartesian:
                return ResultD_Cartesian(radial, x, y, z, self.Ax, self.Ay, self.Az)
            
            return ResultD_Spherical(radial, x, y, z, self.Ax, self.Ay, self.Az)

class Scf:
    eV2Eh = 1/27.211386245981
    Vsip = {
        # Valence state ionization potential, from YAeHMOP
        'H':  {'s': -13.6 * eV2Eh},
        'He': {'s': -23.4 * eV2Eh},
        'Li': {'s': -5.4 * eV2Eh, 'p': -3.5 * eV2Eh},
        'B': {'s': -15.2 * eV2Eh, 'p': -8.5 * eV2Eh},
        'C':  {'s': -21.4 * eV2Eh,  'p': -11.4 * eV2Eh},
        'N':  {'s': -26.0 * eV2Eh,  'p': -13.4 * eV2Eh},
        'O':  {'s': -32.3 * eV2Eh,  'p': -14.8 * eV2Eh},
        'F': {'s': -40.0 * eV2Eh, 'p': -18.1 * eV2Eh},
        'Ne': {'s': -43.2 * eV2Eh,  'p': -20.0 * eV2Eh},
        'Br': {'s': -22.07 * eV2Eh,  'p': -13.1 * eV2Eh},
        'Kr': {'s': 0.0,  'p': 0.0}
    }

    def __init__(self, geometry):
        self.geometry = geometry
        if self.geometry.n_electrons % 2:
            raise ValueError("An even number of electrons is necessary for the closed-shell calculation.")
        
        self.S = Matrix(self.geometry, 0)
        self.basis = self.S.basis_functions
        self.T = Matrix(self.geometry, 1)
        self.V = Matrix(self.geometry, 2)
        self.Rep = Matrix(self.geometry, 3)
        # self.H_huckel = self._extended_huckel()
        self.H_core = self.T.matrix + self.V.matrix
        self.D = self._density(self.H_core)
        self.J, self.K = self._coulomb_and_exchange(self.D)
        self.F = self.H_core + 2*self.J - self.K 

    def scf(self, max_cond_number = 1e-5, max_iter = 100):
        E_nuc = self._nuclear_energy()
        E = [0.0, 0.0]
        D = self.D

        s, U = np.linalg.eigh(self.S.matrix)
        s_crit = 1e-5
        s_mask = s > s_crit
        s = s[s_mask]
        U = U[:, s_mask]
        X = U*(1/np.sqrt(s))

        F_hist = []
        G_hist = []

        for i in range(1, max_iter):
            J, K = self._coulomb_and_exchange(D)
            F = self.H_core + 2*J - K

            # E_elec = np.sum(D*(self.H_core + F))
            One_elec = 2 * np.sum(D*(self.H_core))
            Two_elec = np.sum(D*(F - self.H_core))
            
            Fds = F @ D @ self.S.matrix
            error = Fds - Fds.T
            G = X.T @ error @ X
            G_max = np.max(np.abs(G))

            if G_max < max_cond_number:
                Kin = 2*np.sum(D*self.T.matrix)
                self.J, self.K, self.D, self.F = J, K, D, F
                if self.geometry.charge == 0:
                    print(f"----------------- For {splitext(basename(self.geometry.path_coord))[0]} molecule -----------------")
                else:
                    print(f"----------------- For {splitext(basename(self.geometry.path_coord))[0]}{self.geometry.charge}+ molecule -----------------")
                print("Total SCF energy:", One_elec.item() + Two_elec.item() + E_nuc, "\nOne-electron energy:", One_elec.item(), "\nTwo-electron energy:", Two_elec.item(), "\nNuclear repulsion energy:", E_nuc, "\nKinetic energy:", Kin, "\nConvergence after", i, "iterations")
                mulliken = self._mulliken_charges()
                return [One_elec.item(), Two_elec.item(), E_nuc, Kin]
            
            F_hist.append(F)
            G_hist.append(G)
            if len(F_hist) > 7:
                F_hist.pop(0)
                G_hist.pop(0)

            F_diis = F
            if len(F_hist) > 1:
                dim = len(F_hist)
                B = np.zeros((dim + 1, dim + 1))

                for a in range(dim):
                    for b in range(a, dim):
                        val = np.sum(G_hist[a] * G_hist[b])
                        B[a, b] = val
                        B[b, a] = val
                
                B[-1, :-1] = 1
                B[:-1, -1] = -1
                B[-1, -1] = 0
                
                rhs = np.zeros(dim + 1)
                rhs[-1] = 1
                
                try:
                    c = np.linalg.solve(B, rhs)
                    F_diis = np.zeros_like(F)
                    for idx in range(dim):
                        F_diis += c[idx] * F_hist[idx]
                except np.linalg.LinAlgError:
                    pass
            
            D = self._density(F_diis)
        return 404

    def _mulliken_charges(self):
        charge_density = np.sum(2*self.D * self.S.matrix, axis=1)
        center_j = None
        electronic_charge = 0.0
        current_symbol = ''
        print("\nMulliken charges per centre and basis function type")
        for mu, basis_function in enumerate(self.basis):
            center_i = [basis_function['center']['Ax'], basis_function['center']['Ay'], basis_function['center']['Az']]
            symbol = basis_function['symbol']

            if center_j != center_i:
                if center_j is not None:
                    charge = self.geometry.atomic_number[current_symbol] - electronic_charge
                    print(f"Total Electronic Charge = {electronic_charge:.4f} \nTotal Charge = {charge:.4f} \n")
                
                center_j = center_i
                current_symbol = symbol
                electronic_charge = 0.0
                print(current_symbol)

            q_mu = charge_density[mu]
            print(f'l={basis_function['l']} m={basis_function['m']} n={basis_function['n']}: {q_mu:.4f}')
            electronic_charge += q_mu

        if center_j is not None:
            charge = self.geometry.atomic_number[basis_function['symbol']] - electronic_charge
            print(f"Total Electronic Charge = {electronic_charge:.4f} \nTotal Charge = {charge:.4f} \n")

        return 0

    def _nuclear_energy(self):
        atoms = self.geometry.atom_list
        E_nuc = 0.0

        for i in range(self.geometry.n_atoms):
            for j in range(i + 1, self.geometry.n_atoms):
                Z_a = self.geometry.atomic_number[atoms[i]['symbol']]
                Z_b = self.geometry.atomic_number[atoms[j]['symbol']]

                dx = atoms[i]['Ax'] - atoms[j]['Ax']
                dy = atoms[i]['Ay'] - atoms[j]['Ay']
                dz = atoms[i]['Az'] - atoms[j]['Az']
                d = np.sqrt(dx**2 + dy**2 + dz**2)
                
                E_nuc += (Z_a * Z_b)/d

        return E_nuc
    
    def _extended_huckel(self):
        dimension = len(self.basis)
        K = 1.75
        H_eht = np.zeros((dimension, dimension))
       
        
        for i in range(dimension):
            symbol = self.basis[i]['symbol']
            l, m, n = self.basis[i]['l'], self.basis[i]['m'], self.basis[i]['n']
            H_eht[i,i] = self.Vsip[symbol][self._orbital_type(l,m,n)]

        for i in range(dimension):
            for j in range(i, dimension):
                H_eht[i,j] = K*self.S.matrix[i,j]*(H_eht[i,i] + H_eht[j,j])/2
                H_eht[j,i] = H_eht[i,j]

        return H_eht
    
    def _density(self, H):
        s, U = np.linalg.eigh(self.S.matrix)
        s_crit = 1e-5
        s_mask = s > s_crit
        s = s[s_mask]
        U = U[:, s_mask]
        X = U*(1/np.sqrt(s))
        H_ortho = X.T @ H @ X
        eigenvalues, C = np.linalg.eigh(H_ortho)
        Coeff = X @ C

        N_i = self.geometry.n_electrons // 2
        C_occupied = Coeff[:, :N_i]
        D = C_occupied @ C_occupied.conj().T

        return D
    
    def _coulomb_and_exchange(self, D):
        J = np.einsum('pqrs,rs->pq', self.Rep.matrix, D, optimize=True)
        K = np.einsum('psrq,rs->pq', self.Rep.matrix, D, optimize=True)
        
        return J, K
        
    def _orbital_type(self, l, m, n):
        numbers = l + m + n
        if numbers == 0:
            return 's'
        if numbers == 1:
            return 'p'
        if numbers == 2:
            return 'd'
        if numbers == 3:
            return 'f'

class ResultP:
    """
    As p-orbital's value is the product of its radial part and the coordinate
    component, the class ResultP calculates the final px, py and pz values.
    """

    def __init__(self, value, x, y, z, Ax, Ay, Az):
        """
        Initialize the ResultP object.

        Args:
            value (float): the radial part of p-orbital already calculated.
            x, y, z (float): coordinates of evaluation point.
        """
        self.value = value
        self.dx = x - Ax
        self.dy = y - Ay
        self.dz = z - Az

    def x(self):
        """
        Returns the value for px orbital
        """
        return self.value * self.dx
        
    def y(self):
        """
        Returns the value for py orbital
        """
        return self.value * self.dy
        
    def z(self):
        """
        Returns the value for pz orbital
        """
        return self.value * self.dz

class ResultD_Spherical:
    """
    Calculates the 5 Spherical d-orbitals: xx, yy, zz, xy, xz, yz.
    """

    def __init__(self, value, x, y, z, Ax, Ay, Az):
        """
        Initialize the ResultD_Spherical object.

        Args:
            value (float): the radial part of d-orbital already calculated.
            x, y, z (float): coordinates of evaluation point.
            Ax, Ay, Az (float): coordinates of the orbital center.
        """
        self.value = value
        self.dx = x - Ax
        self.dy = y - Ay
        self.dz = z - Az

        self.dx2 = self.dx ** 2
        self.dy2 = self.dy ** 2
        self.dz2 = self.dz ** 2

    def z2(self):
        angular = self.dz2 - 0.5*(self.dx2 + self.dy2)
        return self.value * angular
        
    def xz(self):
        angular = np.sqrt(3) * self.dx * self.dz
        return self.value * angular
        
    def yz(self):
        angular = np.sqrt(3) * self.dy * self.dz
        return self.value * angular
        
    def x2_y2(self):
        angular = (0.5 * np.sqrt(3)) *(self.dx2 - self.dy2)
        return self.value * angular
        
    def xy(self):
        angular = np.sqrt(3) * self.dx * self.dy
        return self.value * angular

class ResultD_Cartesian:
    """
    Calculates the 6 Cartesian d-orbitals: xx, yy, zz, xy, xz, yz.
    """
    def __init__(self, value, x, y, z, Ax, Ay, Az):
        self.value = value
        self.dx = x - Ax
        self.dy = y - Ay
        self.dz = z - Az

        self.dx2 = self.dx ** 2
        self.dy2 = self.dy ** 2
        self.dz2 = self.dz ** 2

    def xx(self):
        return self.value * self.dx2
        
    def yy(self):
        return self.value * self.dy2
        
    def zz(self):
        return self.value * self.dz2
        
    def xy(self):
        return self.value * (self.dx * self.dy)
        
    def xz(self):
        return self.value * (self.dx * self.dz)
        
    def yz(self):
        return self.value * (self.dy * self.dz)

def factorial2(n):
    value = 1
    for i in range (n, 0, -2):
        value *= i
    return value