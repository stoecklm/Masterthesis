#!/bin/bash

#----------------------------------------------------------------------#
export SCAFESCC_ENABLE_MPI="yes";
export SCAFESCC_ENABLE_OPENMP="yes";
export SCAFESCC_ENABLE_GMOCK="yes";
export SCAFESCC_ENABLE_DEBUG_MODE="no";
export SCAFESCC_ENABLE_CUDA="no";
export SCAFESCC_ENABLE_BOOST_MPI="yes";
export SCAFESCC_ENABLE_NETCDF="yes";
export SCAFESCC_ENABLE_ADOLC="yes";
export SCAFESCC_ENABLE_VT="no";
export SCAFESCC_ENABLE_DOXYGEN="no";
export SCAFESCC_ENABLE_GCOV="no";
#----------------------------------------------------------------------#
module purge
module load autotools
#----------------------------------------------------------------------#
module load boost/1.58.0-gnu5.1
export SCAFESCC_BOOST_ROOT="/sw/taurus/libraries/boost/1.58.0-gnu5.1";
export SCAFESCC_BOOST_LIB="$SCAFESCC_BOOST_ROOT/lib";
#----------------------------------------------------------------------#
module load bullxmpi/1.2.8.4
export SCAFESCC_MPI_ROOT="/sw/taurus/libraries/bullxmpi/1.2.8.4";
#----------------------------------------------------------------------#
module load netcdf/4.3.3.1-gcc-5.2.0
export SCAFESCC_NETCDF_ROOT="/sw/taurus/libraries/netcdf/4.3.3.1-gcc-5.2.0";
export SCAFESCC_HDF5_ROOT="/sw/taurus/libraries/hdf5/1.8.15-gcc-5.10";
#----------------------------------------------------------------------#
module load adolc/2.5.0
export SCAFESCC_ADOLC_ROOT="/sw/taurus/libraries/adolc/2.5.0";
#----------------------------------------------------------------------#
module load gmock/1.6.0
export SCAFESCC_GMOCK_ROOT="/sw/taurus/tools/gmock/1.6.0";
#----------------------------------------------------------------------#
module load doxygen/1.8.11
export SCAFESCC_DOXYGEN_ROOT="/sw/taurus/tools/doxygen/1.8.11";
#----------------------------------------------------------------------#
module load python/3.6-anaconda4.4.0
#----------------------------------------------------------------------#
module load gcc/5.3.0
export SCAFESCC_CXX='g++';
export SCAFESCC_CXXFLAGS="-DNDEBUG -Wall -Wextra -Wno-unused -O3";
#----------------------------------------------------------------------#
module list
