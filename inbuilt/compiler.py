import os
import sys

def compiler_core():
    path = os.getcwd()
    os.chdir(os.path.join(os.path.abspath(sys.path[0]), '%s/FORTRAN_CORE_CODE/' % path))
    os.system("./compile")
    os.chdir(path)