import numpy as np

class Geometry:
    """
    The class Geometry stores the 3D coordinates and atomic composition of the molecule.

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
    } # This class attribute maps chemical symbols to their respective atomic numbers.

    def __init__(self, path):
        """
        Initializes the Geometry object.

        Args:
            path (str): The file path to the atomic coordinates file.
        """
        self.path = path    
        self.atom = {}
    
    def read(self):
        """
        Reads the file and stores the center coordinates of each atom,
        as well as its respective atomic number.
        """
        try:
            with open(self.path, 'r') as file:
                lines = list(file)  
                number_of_atoms = int(lines[0]) - 1
                for i in range(number_of_atoms):
                    l = lines[i+1].strip().split()
                    self.atom[i] = {
                        'atomic_number' : int(self.atomic_number[l[0]]), # Converts symbol to atomic number.
                        'Ax' : float(l[1]),
                        'Ay' : float(l[2]),
                        'Az' : float(l[3])
                    }
                    # Stores the data for the current atom in the dictionary.

            return True
        except FileNotFoundError:
            print(f"Error: The file '{self.path}' wasn't found.")
        except Exception as e:
            print(f'Something went wrong when trying to read the file: {e}')
                    
        return False

class Atom:
    """
    The class Atom calculates Cartesians-Gaussians-Type Orbitals (GTOs).

    An instance of this class represents a single atom of the molecule.
    """
    def __init__(self, path, Ax=0, Ay=0, Az=0):
        """
        Initializes the Atom object.
        
        Args:
            path (str): file path to the basis set.
            Ax, Ay, Az (float): coordinates of the atom.
            gaussians (dict): A dictionary to cache basis set data read from the file.
        """
        self.path = path
        self.Ax = Ax
        self.Ay = Ay
        self.Az = Az
        self.gaussians = {}

    def _read(self, orbital):
        """
        Read the file and stores info (mainly the exponents and coefficients) for given orbital.
        
        Args:
            orbital (str): The orbital itself (s, p, d or f).

        Returns:
            bool: True if the file was read successfully or data was already cached, False otherwise.
        """
        if orbital not in self.gaussians: 
            # Prevents redundancy, reads the file only if the orbital hasn't been already read.
            exp = [] 
            coef = []
            try:
                with open(self.path, 'r') as file:
                    lines = list(file)
                    for i, line in enumerate(lines):
                        if f'{orbital}-type' in line.lower():
                            # Search the line that specifies the orbital type (e.g. 'S-type', 'p-type', etc.).
                            num_functions = int(lines[i+1].strip().split()[0])
                            num_options = int(lines[i+1].strip().split()[1])
                            
                            # Determine line ranges and parse exponents and coefficients.
                            start_exp = i + 2
                            end_exp = start_exp + num_functions
                            for l in lines[start_exp: end_exp]:
                                exp.append(l.strip())
                            
                            start_coef = end_exp
                            end_coef = start_coef + num_functions
                            for l in lines[start_coef: end_coef]:
                                coef.append(l.strip().split())
                
                # Dictionary stores the parsed data.
                self.gaussians[orbital] = {
                    'number of functions': num_functions,
                    'options': num_options,
                    'exponents': exp,
                    'coefficients': coef
                }

                return True
            except FileNotFoundError:
                # Exception if the file doesn't exist in the path mentioned.
                print(f"Error: The file '{self.path}' wasn't found.")
            
            except Exception as e:
                # Exception if happens some other error.
                print(f'Something went wrong when trying to read the file: {e}')
                    
            return False
        else:
            return True
    
    def _calculate(self, x, y, z, sel, orbital):
        """
        After reading, calculates the value of the GTOs.
        
        Args:
            x, y, z (float): coordinates of evaluation point.
            sel (int): select the contraction option.
            orbital (str): the orbital itself (s, p, d or f).
        
        Returns:
            float: The value of the contracted orbitals.
        """
        scan = self._read(orbital)

        if not scan:
            print('Something went wrong when trying to read the file')
            return False

        data = self.gaussians[orbital]
        coef = []
        
        if 0 < sel <= data['options']:
            # Validates the selection.
            for c in data['coefficients']:
                coef.append(c[sel - 1])
        else: 
            raise ValueError(f'This selection does not exist, choose from 1 to {data['options']}')
            # Raises if selected option is out of bounds.
            
        exp = data['exponents']
        contracted_value = 0.0

        r2 = (self.Ax - x)**2 + (self.Ay - y)**2 + (self.Az - z)**2
        # Calculate the squared distance from the atom's center.
        for i in range(data['number of functions']):
            if coef[i] != 0.0:
                # Skips calculation if the coefficient is zero.
                contracted_value += float(coef[i]) * np.exp(-float(exp[i]) * r2)
                # Sum of each primitive: C * e^(-a * r^2).
        return contracted_value
    
    def s(self, x, y, z, sel=1):
        """
        Returns the value for s-orbital.
        """
        result = self._calculate(x, y, z, sel, 's')
        return result
    
    def p(self, x, y, z, sel=1):
        """
        Returns the value for p-orbitals.
        Necessary to specify the axis of projection, ResultP object calculates for the determined component.
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
