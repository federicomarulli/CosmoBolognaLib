#!/usr/bin/env python

from setuptools import setup, Extension
import os
import platform

from distutils.sysconfig import get_config_vars

(opt,) = get_config_vars('OPT')
os.environ['OPT'] = " ".join(flag for flag in opt.split() if flag != '-Wstrict-prototypes')

def readme():
    with open('README.rst') as f:
        return f.read()

HOME = os.getenv("HOME")

dirLib       = "CosmoBolognaLib/" 
dirH         = dirLib+"Headers/Lib/"
dirO         = dirLib+"Headers/Objects/"
dirEH        = dirLib+"External/EH/"
dir_CUBA     = "../External/Cuba-4.2/"
dir_FFTLOG   = "../External/fftlog-f90-master/"

FLAGS = ["-std=c++11", "-fopenmp", "-w"]
FLAGSL = ["-Wl,-rpath,CosmoBolognaLib/", "-LCosmoBolognaLib/", "-Wl,-rpath,../External/Cuba-4.2/", "-L../External/Cuba-4.2/"]

include_dirs = [dirLib, dirH, dirO, dirEH, dir_CUBA, dir_FFTLOG]

libraries = ["gomp", "gsl", "gslcblas", "m", "fftw3", "fftw3_omp", "cuba"]
print platform.system()

if platform.system()=='Darwin':
    os.environ["CC"] = "gcc"
    os.environ["CXX"] = "g++"
    os.environ["MPICXX"] = "mpic++"
    FLAGS=['-arch','x86_64',"-std=c++11", "-fopenmp", "-w"]
    FLAGSL=['-arch','x86_64',"-Wl,-rpath,CosmoBolognaLib/", "-LCosmoBolognaLib/", "-Wl,-rpath,../External/Cuba-4.2/", "-L../External/Cuba-4.2/"]
    libraries = ["gomp", "gsl", "gslcblas", "m", "fftw3", "cuba"]

sources = ["CBL_wrap.cxx",
           dirLib+"Func/CUBAwrapper.cpp",
           dirLib+"Func/Data1D_collection.cpp",
           dirLib+"Func/Data1D.cpp",
           dirLib+"Func/Data1D_extra.cpp",
           dirLib+"Func/Data2D.cpp",
           dirLib+"Func/Data2D_extra.cpp",
           dirLib+"Func/Data.cpp",
           dirLib+"Func/Distribution.cpp",
           dirLib+"Func/FFTlog.cpp",
           dirLib+"Func/Field3D.cpp",
           dirLib+"Func/Func.cpp",
           dirLib+"Func/FuncGrid.cpp",
           dirLib+"Func/FuncMultipoles.cpp",
           dirLib+"Func/FuncXi.cpp",
           dirLib+"Func/GSLfunction.cpp",
           dirLib+"Func/GSLwrapper.cpp",
           dirLib+"Statistics/Chain.cpp",
           dirLib+"Statistics/BaseParameter.cpp",
           dirLib+"Statistics/DerivedParameter.cpp",
           dirLib+"Statistics/Likelihood.cpp",
           dirLib+"Statistics/LikelihoodFunction.cpp",
           dirLib+"Statistics/LikelihoodParameters.cpp",
           dirLib+"Statistics/Model1D.cpp",
           dirLib+"Statistics/Model2D.cpp",
           dirLib+"Statistics/Model.cpp",
           dirLib+"Statistics/Parameter.cpp",
           dirLib+"Statistics/Sampler.cpp",
           dirLib+"External/EH/power_whu.cpp",
           dirLib+"Cosmology/Lib/BAO.cpp",
           dirLib+"Cosmology/Lib/Bias.cpp",
           dirLib+"Cosmology/Lib/Cosmology.cpp",
           dirLib+"Cosmology/Lib/DensityProfile.cpp",
           dirLib+"Cosmology/Lib/MassFunction.cpp",
           dirLib+"Cosmology/Lib/MassGrowth.cpp",
           dirLib+"Cosmology/Lib/NG.cpp",
           dirLib+"Cosmology/Lib/PkXi.cpp",
           dirLib+"Cosmology/Lib/PkXizSpace.cpp",
           dirLib+"Cosmology/Lib/RSD.cpp",
           dirLib+"Cosmology/Lib/Sigma.cpp",
           dirLib+"Cosmology/Lib/SizeFunction.cpp",
           dirLib+"Cosmology/Lib/Velocities.cpp",
           dirLib+"Cosmology/Lib/3PCF.cpp",
           dirLib+"ChainMesh/ChainMesh.cpp",
           dirLib+"Catalogue/Catalogue.cpp",
           dirLib+"Catalogue/ChainMesh_Catalogue.cpp",
           dirLib+"Catalogue/GadgetCatalogue.cpp",
           dirLib+"Catalogue/Object.cpp",
           dirLib+"Catalogue/RandomCatalogue.cpp",
           dirLib+"Catalogue/RandomCatalogueVIPERS.cpp",
           dirLib+"Catalogue/VoidCatalogue.cpp",
           dirLib+"LogNormal/LogNormal.cpp",
           dirLib+"LogNormal/LogNormalFull.cpp",
           dirLib+"Measure/TwoPointCorrelation/Pair1D.cpp",
           dirLib+"Measure/TwoPointCorrelation/Pair1D_extra.cpp",
           dirLib+"Measure/TwoPointCorrelation/Pair2D.cpp",
           dirLib+"Measure/TwoPointCorrelation/Pair2D_extra.cpp",
           dirLib+"Measure/TwoPointCorrelation/Pair.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelation1D_angular.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelation1D.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelation1D_filtered.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelation1D_monopole.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelation2D_cartesian.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelation2D.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelation2D_polar.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelation.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelation_deprojected.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelation_multipoles.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelation_projected.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelation_wedges.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelationCross1D.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelationCross1D_monopole.cpp",
           dirLib+"Measure/TwoPointCorrelation/TwoPointCorrelationCross.cpp",
           dirLib+"Measure/ThreePointCorrelation/ThreePointCorrelation_angular_connected.cpp",
           dirLib+"Measure/ThreePointCorrelation/ThreePointCorrelation_angular_reduced.cpp",
           dirLib+"Measure/ThreePointCorrelation/ThreePointCorrelation_comoving_connected.cpp",
           dirLib+"Measure/ThreePointCorrelation/ThreePointCorrelation_comoving_reduced.cpp",
           dirLib+"Measure/ThreePointCorrelation/ThreePointCorrelation.cpp",
           dirLib+"Measure/ThreePointCorrelation/Triplet.cpp",
           dirLib+"Modelling/Global/Modelling.cpp",
           dirLib+"Modelling/Cosmology/ModelFunction_Cosmology.cpp",
           dirLib+"Modelling/Cosmology/Modelling_Cosmology.cpp",
           dirLib+"Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation1D_angular.cpp",
           dirLib+"Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation2D_polar.cpp",
           dirLib+"Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation1D.cpp",
           dirLib+"Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation.cpp",
           dirLib+"Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation1D_filtered.cpp",
           dirLib+"Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation_deprojected.cpp",
           dirLib+"Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation1D_monopole.cpp",
           dirLib+"Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation_multipoles.cpp",
           dirLib+"Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation2D_cartesian.cpp",
           dirLib+"Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation_projected.cpp",
           dirLib+"Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation2D.cpp",
           dirLib+"Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation_wedges.cpp",
           dirLib+"Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation1D_angular.cpp",
           dirLib+"Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation2D_cartesian.cpp",
           dirLib+"Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation_deprojected.cpp",
           dirLib+"Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation1D.cpp",
           dirLib+"Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation2D.cpp",
           dirLib+"Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation_multipoles.cpp",
           dirLib+"Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation1D_filtered.cpp",
           dirLib+"Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation2D_polar.cpp",
           dirLib+"Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation_projected.cpp",
           dirLib+"Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation1D_monopole.cpp",
           dirLib+"Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation.cpp",
           dirLib+"Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation_wedges.cpp",
           dirLib+"Modelling/ThreePointCorrelation/ModelFunction_ThreePointCorrelation_angular_connected.cpp",
           dirLib+"Modelling/ThreePointCorrelation/ModelFunction_ThreePointCorrelation_angular_reduced.cpp",
           dirLib+"Modelling/ThreePointCorrelation/ModelFunction_ThreePointCorrelation_comoving_connected.cpp",
           dirLib+"Modelling/ThreePointCorrelation/ModelFunction_ThreePointCorrelation_comoving_reduced.cpp",
           dirLib+"Modelling/ThreePointCorrelation/ModelFunction_ThreePointCorrelation.cpp",
           dirLib+"Modelling/ThreePointCorrelation/Modelling_ThreePointCorrelation_angular_connected.cpp",
           dirLib+"Modelling/ThreePointCorrelation/Modelling_ThreePointCorrelation_angular_reduced.cpp",
           dirLib+"Modelling/ThreePointCorrelation/Modelling_ThreePointCorrelation_comoving_connected.cpp",
           dirLib+"Modelling/ThreePointCorrelation/Modelling_ThreePointCorrelation_comoving_reduced.cpp",
           dirLib+"Modelling/ThreePointCorrelation/Modelling_ThreePointCorrelation.cpp",
           dirLib+"GlobalFunc/FuncCosmology.cpp",
           dirLib+"GlobalFunc/Func.cpp",
           dirLib+"GlobalFunc/Reconstruction.cpp",
           dirLib+"GlobalFunc/SubSample.cpp",
           dirLib+"ReadParameters/ReadParameters.cpp"
       ]



CosmoBolognaLib = Extension(  "_CosmoBolognaLib",
                              language             = "c++",
                              sources              = sources,
                              include_dirs         = include_dirs,
                              libraries            = libraries,
                              extra_compile_args   = FLAGS,
                              extra_link_args      = FLAGSL )

setup(  name             = "CosmoBolognaLib",
        version          = "3.2",
        description      = "C++ libraries for cosmological calculations",
        long_description = readme(),
        author           = "Federico Marulli",
        author_email     = "federico.marulli3@unibo.it",
        url              = "http://github.com/federicomarulli/CosmoBolognaLib",
        license          = "GNU General Public License",
        ext_modules      = [CosmoBolognaLib],
        packages         = ["CosmoBolognaLib"],
        zip_safe         = False )

