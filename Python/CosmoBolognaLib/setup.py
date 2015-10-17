from distutils.core import setup,Extension
from Cython.Build import cythonize
import os

name = "CosmoBolognaLib"
HOME = os.getenv("HOME")

dirLib = "../../" 
dirNumerical = HOME+"/Numerical/"
dirH = dirLib+"Headers/Lib/"
dirO = dirLib+"Headers/Objects/"
dirEH = dirLib+"Cosmology/EH/"

dirCosmology=dirLib+"Cosmology/Lib/"

FLAGS = ["-std=c++11","-fopenmp","-O3","-unroll","-ftree-loop-distribution"]
FLAGSL = ["-DSOME_DEFINE_OPT","-fopenmp","-L../../"]

include_dirs=[dirNumerical,dirLib,dirH,dirO,dirEH]

runtime_library_dirs=[dirLib]
libraries =["gomp","FUNC","STAT","COSM"]
sources = ["CosmoBolognaLib.pyx","external.cpp",
           dirCosmology+"Cosmology.cpp",
           dirCosmology+"PkXi.cpp",
	   dirCosmology+"BAO.cpp"]

l1 = '#include "Func.h"\n'
l2 = 'const string cosmobl::par::DirCosmo = "'+os.path.abspath("../../")+'/";\n'
l3 = 'const string cosmobl::par::DirLoc = "./";\n'

outf = open("external.cpp","w")
outf.write(l1+l2+l3)
outf.close()




setup(name = name,
      version = "0.0.1",
      ext_modules=cythonize(Extension(
          name,
          include_dirs=include_dirs,
          runtime_library_dirs=runtime_library_dirs,
          libraries=libraries,
          sources = sources,
          language="c++",
          extra_compile_args=FLAGS,
          extra_link_args=FLAGSL
          )))
