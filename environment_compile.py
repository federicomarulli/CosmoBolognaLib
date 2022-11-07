import os
import sys

if len(sys.argv) > 1:
        command = str(sys.argv[1])
else:
        command = ""

if os.environ['CONDA_DEFAULT_ENV'] != "cbl":
        print("\033[1;33m ERROR: the conda environment \'cbl\' must be active!\n")
        sys.exit()

cmd = "%s "

cmd += " CC=$GCC CXX=$GXX F=$GFORTRAN CXX_OLD=$GXX "
cmd += " dir_INC_FFTW=$CONDA_BUILD_SYSROOT/include "
cmd += " dir_CCfits_root=$CONDA_PREFIX/include "
cmd += " dir_CUBA_root=$CONDA_PREFIX/include "
cmd += " dir_Eigen_root=$CONDA_PREFIX/include/eigen3 "
cmd += " dir_INC_BOOST=$CONDA_PREFIX/include "
cmd += " dir_INC_GSL=$CONDA_PREFIX/include "
cmd += " dir_LIB_FFTW=$CONDA_PREFIX/lib "
cmd += " FLAGS_LINK='-shared -Wl,-rpath,/../../lib -Wl,-rpath,./' "
cmd += "%s"

os.system(cmd%("make", command))

