import numpy as np
from math import  comb
from sympy import factorial, factorial2

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

    

    def __init__(self, path_coord, path_basis):
        """
        Initializes the Geometry object.

        Args:
            path_coord (str): The file path to the atomic coordinates file.
            path_basis (str): The file path to the basis set.
        """
        self.path_coord = path_coord
        self.path_basis = path_basis
        self.atoms = {}  # A dictionary to cache the data read from the files.
        self.n_atoms = 0 # Number of atoms in the element
        self.n_eletrons = 0

        # Reads the file of the atomic coordinates, sets self.atoms.
        first_file = self._read_coordinates()
        
        # Reads the file of the basis set and stores the data of the atom in the gaussians dictionary.
        if first_file:    
            second_file = self._read_basis_set()

        n_eletrons = 0
        if first_file and second_file:
            # Loop through each element symbol and its atoms.
            for symbol, element_list in self.atoms.items():
                # Loop through each atom in the element.
                n_eletrons += self.atomic_number[symbol]*len(element_list['instance'])
                for atom_instance in element_list['instance']:
                    # Create a Basis object for this atom and add it to the dictionary.
                    atom_instance['GTO'] = Basis(atom_instance, element_list['basis'])

            self.n_eletrons = n_eletrons

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
                        'Ax': float(l[1]),
                        'Ay': float(l[2]),
                        'Az': float(l[3])
                    }

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
        Docstring for __init__
        
        :param self: Description
        :param geometry: Description
        :param type_matrix: 0: Overlap, 1: Kinetic 
        """
        self.geometry = geometry
        self.type_matrix = type_matrix

        basis_functions = []
        
        molecule = []
        for symbol, data in self.geometry.atoms.items():
            for center in data['instance']:
                molecule.append({
                    'center': center,
                    'basis': data['basis']
                })
                
        
        for atom in molecule:
            
            for orbital, data in atom['basis'].items():
                quant_numbers = self.Quantum_numbers.get(orbital, [])
                
                for k in range(data['contracted_functions']):
                    
                    for (l, m, n) in quant_numbers:
                        
                        function = {
                            'center': atom['center'],
                            'l': l,
                            'm': m,
                            'n': n,
                            'primitives': [],
                        }
                        
                        for i in range(data['number_of_primitives']):
                            coef = float(data['coefficients'][i][k])
                            if coef != 0.0: 
                                function['primitives'].append({
                                    'exponent': data['exponents'][i],
                                    'coefficient': coef
                                })
                        
                        basis_functions.append(function)

        dimension = len(basis_functions)
        if type_matrix <= 2:
            self.matrix = np.zeros((dimension, dimension))
            self.normalised_matrix = np.zeros((dimension, dimension))
            
            for i in range(dimension):
                for j in range(i, dimension):
                    
                    if type_matrix <= 1:
                        value = self._integral(basis_functions[i], basis_functions[j], type_matrix)

                    if type_matrix == 0:
                        if i == j:
                            self.normalised_matrix[i, j] = 1
                        else:
                            normalised_value = self._integral(basis_functions[i], basis_functions[j], 1)
                            self.normalised_matrix[i, j] = normalised_value
                            self.normalised_matrix[j, i] = normalised_value

                    if type_matrix == 2:
                        value = 0
                        for atom in molecule:
                            value += self._integral(basis_functions[i], basis_functions[j], type_matrix, 0, atom)


                    self.matrix[i, j] = value
                    if i != j:
                        self.matrix[j, i] = value
        
        if type_matrix == 3:
            self.matrix = np.zeros((dimension, dimension, dimension, dimension))

            for i in range(dimension):
                for j in range(i, dimension):
                    for k in range(j, dimension):
                        for l in range(k, dimension):
                            value = self._integral4(basis_functions[i], basis_functions[j], basis_functions[k], basis_functions[l])

                            self.matrix[i, j, k, l] = value

    def _integral4(self, eta1, eta2, eta3, eta4):
        value = 0.0
    
        for primitive1 in eta1['primitives']:
            for primitive2 in eta2['primitives']:
                for primitive3 in eta3['primitives']:
                    for primitive4 in eta4['primitives']:
                        value += self._repulsion(eta1, eta2, eta3, eta4, primitive1, primitive2, primitive3, primitive4)
    
        return value

    def _integral(self, eta1, eta2, type_matrix = 0, normalised = 0, nucleus = 0):
        value = 0.0
    
        for primitive1 in eta1['primitives']:
            for primitive2 in eta2['primitives']:

                if type_matrix == 0:
                    value += self._overlap(eta1, eta2, primitive1, primitive2, normalised)
            
                if type_matrix == 1:
                    value += self._kinetic_energy(eta1, eta2, primitive1, primitive2)    

                if type_matrix == 2:
                    value += self._nuclear_attraction(eta1, eta2, primitive1, primitive2, nucleus)


        return value
    
    def _overlap(self, eta1, eta2, p1, p2, normalised = 0):
        value = 0.0
        A = eta1['center']
        B = eta2['center']

        ab2 = (A['Ax'] - B['Ax'])**2 + (A['Ay'] - B['Ay'])**2 + (A['Az'] - B['Az'])**2 

        C_a = p1['coefficient']
        alpha = p1['exponent']

        C_b = p2['coefficient']
        beta = p2['exponent']
        gamma = alpha + beta

        Px = (alpha*A['Ax'] + beta*B['Ax'])/gamma
        Py = (alpha*A['Ay'] + beta*B['Ay'])/gamma
        Pz = (alpha*A['Az'] + beta*B['Az'])/gamma
        
        Sx = self._S(eta1['l'], eta2['l'], A['Ax'] - Px, B['Ax'] - Px,gamma)
        Sy = self._S(eta1['m'], eta2['m'], A['Ay'] - Py, B['Ay'] - Py,gamma)
        Sz = self._S(eta1['n'], eta2['n'], A['Az'] - Pz, B['Az'] - Pz,gamma)
        
        integral = np.exp(-alpha*beta*ab2/gamma)*Sx*Sy*Sz

        if normalised == 1:
            value += C_a*C_b*self._N(alpha, eta1['l'], eta1['m'], eta1['n'])*self._N(beta, eta2['l'], eta2['m'], eta2['n'])*integral
            
        else:    
            value += C_a*C_b*integral

        return value
    
    def _kinetic_energy(self, eta1, eta2, p1, p2):
        alpha = p1['exponent']
        value = alpha*(2*(eta1['l'] + eta1['m'] + eta1['n']) + 3)*self._overlap(eta1, eta2, p1, p2)
        coord = ['l', 'm', 'n']
        
        for c in coord:
            eta1[c] += 2
            value += -2*(alpha**2)*self._overlap(eta1, eta2, p1, p2)
            eta1[c] -= 2

        
        for c in coord:
            if eta1[c] >= 2:
                eta1[c] -= 2
                value += -0.5*((eta1[c] + 2)*(eta1[c] + 1))*self._overlap(eta1, eta2, p1, p2)
                eta1[c] += 2

        return value
        
    def _nuclear_attraction(self, eta1, eta2, p1, p2, nucleus):
        A = eta1['center']
        B = eta2['center']
        C = nucleus['center']
        ab2 = (A['Ax'] - B['Ax'])**2 + (A['Ay'] - B['Ay'])**2 + (A['Az'] - B['Az'])**2

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

        A_x = self._A('l', eta1, eta2, A['Ax'], B['Ax'], C['Ax'], P['Px'], eps)
        A_y = self._A('m', eta1, eta2, A['Ay'], B['Ay'], C['Ay'], P['Py'], eps)
        A_z = self._A('n', eta1, eta2, A['Az'], B['Az'], C['Az'], P['Pz'], eps)

        v_max = eta1['l'] + eta1['m'] + eta1['n'] + eta2['l'] + eta2['m'] + eta2['n']

        boys = self._Boys_recurrence(v_max, gamma*pc2)

        K = 2*np.pi*np.exp(-alpha*beta*ab2/gamma)/gamma

        value = 0
        for x in range(len(A_x)):
            for y in range(len(A_y)):
                for z in range(len(A_z)):
                    value += A_x[x]*A_y[y]*A_z[z]*boys[x+y+z]

        return K*value
    
    def _repulsion(self, eta1, eta2, eta3, eta4, p1, p2, p3, p4):
        A = eta1['center']
        B = eta2['center']
        C = eta3['center']
        D = eta4['center']
        ab2 = (A['Ax'] - B['Ax'])**2 + (A['Ay'] - B['Ay'])**2 + (A['Az'] - B['Az'])**2
        cd2 = (C['Ax'] - D['Ax'])**2 + (C['Ay'] - D['Ay'])**2 + (C['Az'] - D['Az'])**2

        alpha1 = p1['exponent']
        beta1 = p2['exponent']
        alpha2 = p3['exponent']
        beta2 = p4['exponent']
        gamma1 = alpha1 + beta1
        gamma2 = alpha2 + beta2
        delta = 1/(4*gamma1) + 1/(4*gamma2)

        P = {
            'Px' : (alpha1*A['Ax'] + beta1*B['Ax'])/gamma1,
            'Py' : (alpha1*A['Ay'] + beta1*B['Ay'])/gamma1,
            'Pz' : (alpha1*A['Az'] + beta1*B['Az'])/gamma1
        }
        Q = {
            'Qx' : (alpha2*C['Ax'] + beta2*D['Ax'])/gamma2,
            'Qy' : (alpha2*C['Ay'] + beta2*D['Ay'])/gamma2,
            'Qz' : (alpha2*C['Az'] + beta2*D['Az'])/gamma2
        }
        p2 = (P['Px'] - Q['Qx'])**2 + (P['Py'] - Q['Qy'])**2 + (P['Pz'] - Q['Qz'])**2

        B_x = self._B('l', eta1, eta2, eta3, eta4, A['Ax'], B['Ax'], C['Ax'], D['Ax'], P['Px'], Q['Qx'], gamma1, gamma2, delta)
        B_y = self._B('m', eta1, eta2, eta3, eta4, A['Ay'], B['Ay'], C['Ay'], D['Ay'], P['Py'], Q['Qy'], gamma1, gamma2, delta)
        B_z = self._B('n', eta1, eta2, eta3, eta4, A['Az'], B['Az'], C['Az'], D['Az'], P['Pz'], Q['Qz'], gamma1, gamma2, delta)

        v_max = eta1['l'] + eta1['m'] + eta1['n'] + eta2['l'] + eta2['m'] + eta2['n'] + eta3['l'] + eta3['m'] + eta3['n'] + eta4['l'] + eta4['m'] + eta4['n']

        boys = self._Boys_recurrence(v_max, p2/(4*delta))

        omega = 2*np.pi**2*np.sqrt(np.pi/(gamma1 + gamma2))*np.exp(-(alpha1*beta1*ab2/gamma1 + alpha2*beta2*cd2/gamma2))/(gamma1*gamma2)

        value = 0
        for x in range(len(B_x)):
            for y in range(len(B_y)):
                for z in range(len(B_z)):
                    value += B_x[x]*B_y[y]*B_z[z]*boys[x+y+z]

        return omega*value

    def _B(self, lmn, eta1, eta2, eta3, eta4, A, B, C, D, P, Q, gamma1, gamma2, delta):
        
        array = np.zeros(eta1[f'{lmn}'] + eta2[f'{lmn}'] + eta3[f'{lmn}'] + eta4[f'{lmn}'] + 1)

        for index1 in range(0, eta1[f'{lmn}'] + eta2[f'{lmn}'] + 1):
            for index2 in range(0, index1//2 + 1):
                for index3 in range(0, (index1 - 2*index2)//2 + 1):
                    for index4 in range(0, eta3[f'{lmn}'] + eta4[f'{lmn}'] + 1):
                        for index5 in range(0, index4//2 + 1):
                            array[index1 + index4 - 2*(index2 + index5) - index3] += (-1)**(index4 + index3)*self._theta(index1, eta1[f'{lmn}'], eta2[f'{lmn}'], A - P, B - P, index2, gamma1)*self._theta(index4, eta3[f'{lmn}'], eta4[f'{lmn}'], C - Q, D - Q, index5, gamma2)*(2)**(2*(index2+index5-(index1 + index4)))*(delta)**(2*(index2+index5)+index3-(index1+index4))*factorial(index1 + index4 - 2*(index2 + index5))*(P  - Q)**(index1 + index4 - 2*(index3 + index2 + index5))/(factorial(index3)*factorial(index1 + index4 - 2*(index3 + index2 + index5)))

        return array

    def _theta(self, j, l, m, a, b, r, gamma):
        value = self._f_j(l, m, a, b, j)*factorial(j)*gamma**(r-j)/(factorial(r)*factorial(j - 2*r))

        return value
    
    def _A(self, lmn, eta1, eta2, A, B, C, P, eps):
        
        array = np.zeros(eta1[f'{lmn}'] + eta2[f'{lmn}'] + 1)

        for index1 in range(0, eta1[f'{lmn}'] + eta2[f'{lmn}'] + 1):
            for index2 in range(0, index1//2 + 1):
                for index3 in range(0, (index1 - 2*index2)//2 + 1):
                    array[index1 - 2*index2 - index3] += (-1)**(index1 + index3)*self._f_j(eta1[f'{lmn}'], eta2[f'{lmn}'], A - P, B - P, index1)*factorial(index1)*(C - P)**(index1 - 2*(index2 + index3))*eps**(index2 + index3)/(factorial(index2)*factorial(index3)*factorial(index1 - 2*(index2 + index3)))

        return array
    
    
    def _Boys_recurrence(self, n, x):
        exp = np.exp(-x)

        if x < n + 1/2:
            f_max = 1/(2*n+1)
        else:
            f_max = factorial2(2*n-1)*(np.sqrt(np.pi/(x**(2*n+1))))/(2**(n+1))

        array = np.zeros(n + 1)
        array[n] = f_max

        for i in range(n, 0, -1):
            array[i-1] = (2*x*array[i] + exp)/(2*i + 1)

        return array
    
    def _N(self, alpha, l, m ,n):
        N = ((4*alpha)**(l + m + n) / factorial2(2*l - 1)*factorial2(2*m - 1)*factorial2(2*n - 1))**(1/2) * (2*alpha/np.pi)**(3/4)

        return N
    
    def _S(self, l, m, Pa, Pb, g):
        value = 0.0
        
        for j in range((l + m)//2 + 1):
            value += self._f_j(l, m, Pa, Pb, 2*j) * factorial2(2*j - 1) / (2*g)**j

        value *= np.sqrt(np.pi/g)

        return value      

    def _f_j(self, l, m, a, b, j):
        """
        Calculates the coefficient of x^j in the expansion of (x + a)^l * (x + b)^m
        """
        value = 0.0

        for k in range(max(0, j - m), min(j, l) + 1):
            value += comb(l, k) * comb(m, j - k) * a**(l - k) * b**(m + k - j)

        return value

class Basis:
    """
    The class Basis calculates Cartesians-Gaussians-Type Orbitals (GTOs).

    An instance of this class represents a single atom of the molecule.
    """

    def __init__(self, center, basis_data):
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
        return ResultP(radial, x, y, z)

class SCF:
    def __init__(self, path_coord, path_basis, initial_guess=0):
        # initial_guess = 0: Core Hamiltonian
        self.geometry = Geometry(path_coord, path_basis)
        if not self.geometry.n_eletrons % 2:
            raise ValueError("An even number of electrons is necessary for the closed-shell calculation.")
        
        self.S = Matrix(self.geometry, 0)
        self.T = Matrix(self.geometry, 1)
        self.V = Matrix(self.geometry, 2)

        H_core = self.T + self.V

        closed_shell = self.geometry.n_eletrons / 2

        X = self._transform_matrix(self.S)
        H = X.T @ H_core @ X

        eig_val, C0 = np.linalg.eigh(H)
        C = X @ C0


    def _transform_matrix(self, S, s_crit=10**-5):
        s, U = np.linalg.eigh(S)

        tolerance = s > s_crit
        s_new = s[tolerance]
        U_new = U[:, tolerance]

        X = U_new @ np.diag(s_new**-0.5)

        return X


    
class ResultP:
    """
    As p-orbital's value is the product of its radial part and the coordinate
    component, the class ResultP calculates the final px, py and pz values.
    """

    def __init__(self, value, x, y, z):
        """
        Initialize the ResultP object.

        Args:
            value (float): the radial part of p-orbital already calculated.
            x, y, z (float): coordinates of evaluation point.
        """
        self.value = value
        self._x = x
        self._y = y
        self._z = z

    def x(self):
        """
        Returns the value for px orbital
        """
        return self.value * self._x
        
    def y(self):
        """
        Returns the value for py orbital
        """
        return self.value * self._y
        
    def z(self):
        """
        Returns the value for pz orbital
        """
        return self.value * self._z
