import numpy as np
from itertools import permutations, combinations, combinations_with_replacement

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

        # Reads the file of the atomic coordinates, sets self.atoms.
        first_file = self._read_coordinates()
        
        # Reads the file of the basis set and stores the data of the atom in the gaussians dictionary.
        if first_file:    
            second_file = self._read_basis_set()

        if first_file and second_file:
            # Loop through each element symbol and its atoms.
            for symbol, element_list in self.atoms.items():
                # Loop through each atom in the element.
                for atom_instance in element_list['instance']:
                    # Create a Basis object for this atom and add it to the dictionary.
                    atom_instance['GTO'] = Basis(atom_instance, element_list['basis'])
        
        

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

                for i in range(self.n_atoms):
                    l = lines[i+2].strip().split()
                    symbol = l[0]
                    atom_coord = {
                        'Ax': float(l[1]),
                        'Ay': float(l[2]),
                        'Az': float(l[3])
                    }
                    # Stores the data for the current atom in the dictionary.
                    if symbol not in self.atoms:
                        self.atoms[symbol] = {
                            'basis' : {},
                            'instance' : [atom_coord]
                        }
                    else:
                        self.atoms[symbol]['instance'].append(atom_coord)
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
                        num_functions = int(functions_line[0])
                        num_options = int(functions_line[1])

                        # Parse exponents and coefficients.
                        for _ in range(num_functions):
                            exp.append(float(next(lines_iterator).strip()))

                        for _ in range(num_functions):
                            coef.append(next(lines_iterator).strip().split())
                        
                        self.atoms[current_element]['basis'][orbital] = {
                            'number_of_functions': num_functions,
                            'options': num_options,
                            'exponents': exp,
                            'coefficients': coef,
                        } # Dictionary that stores the parsed data.
            return True
   
        except Exception as e:
            # Exception if happens some error.
            print(f'Something went wrong when trying to read the file: {e}')
        return False

class Overlap:
    def __init__(self, geometry):
        self.geometry = geometry
        self.dimension = geometry.n_atoms
        self.matrix = []
        
        elements = list(self.geometry.atoms.keys())
        elements_comb = combinations(elements, 2)
        for elements_pair in elements_comb:
            A_basis = self.geometry.atoms[elements_pair[0]]['basis']
            B_basis = self.geometry.atoms[elements_pair[1]]['basis']
            A_orbitals = list(A_basis.keys())
            B_orbitals = list(B_basis.keys())
            index_b = 0
            last_index = 0
            for i in range(len(A_orbitals)*len(B_orbitals)):
                index_a = i//len(B_orbitals)
                last_index = (i - 1)//len(B_orbitals)
                if i > 0 and index_a != last_index:
                    index_b = 0
                f_a = A_basis[A_orbitals[index_a]]
                f_b = B_basis[B_orbitals[index_b]]
                self.matrix.append(self._integral(f_a,f_b))
                index_b += 1
                

    def _integral(self, f1, f2):
        return print(f1['exponents'], f2['exponents'])

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
        if 0 <= sel < data['options']:
            # Validates the selection.      
            for c in data['coefficients']:
                coef.append(c[sel])
        else:
            raise ValueError(
                f'This selection does not exist, choose from 0 to {data['options'] - 1}')
            # Raises if selected option is out of bounds.
        
        contracted_value = 0.0
        r2 = (self.Ax - x)**2 + (self.Ay - y)**2 + (self.Az - z)**2
        # Calculate the squared distance from the atom's center.
        
        for i in range(data['number_of_functions']):
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
