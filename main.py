from atomos import Atom
from atomos import Geometry

arquivo1 = 'C:\\Users\\caion\\OneDrive\\Documentos\\programation\\python\\ic\\atom.txt'
arquivo2 = 'C:\\Users\\caion\\OneDrive\\Documentos\\programation\\python\\ic\\mol.txt'

instancia = Atom(arquivo1)
geometry = Geometry(arquivo2)
print(instancia.s(0.5,2,1))
