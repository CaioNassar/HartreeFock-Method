from atoms import Atom
from atoms import Geometry

arquivo1 = 'C:\\Users\\caion\\OneDrive\\Documentos\\programming\\python\\ic\\atom.txt'
arquivo2 = 'C:\\Users\\caion\\OneDrive\\Documentos\\programming\\python\\ic\\mol.txt'

instancia = Atom(arquivo1, 0.5, 1, 0.3)
geometry = Geometry(arquivo2)
print(instancia.s())
