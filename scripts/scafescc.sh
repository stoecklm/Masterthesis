#!/bin/bash
#
# Call: ./scafescc.sh <source.cpp>
# Parameters: source.cpp: Main source file containing ScaFES related objects.
#
# Shell script for compiling and linking ScaFES applications.
# The shell script can be configured via environment variables prefixed with
# SCAFESCC. All environment variables are listed in the following
# and have to be exported.
#
# * C++11 compiler. <g++, icpc, pgc++...> default=g++
#   SCAFESCC_CXX='g++';
#
# * Compiler flags. default="-DNDEBUG -O3 -Wall"
#   SCAFESCC_CXXFLAGS="-DNDEBUG -O3 -Wall";
#
# * Enable MPI? <no/yes>. default=yes
#   SCAFESCC_ENABLE_MPI="yes";
#
# * Only needed if SCAFESCC_ENABLE_MPI="yes": default=/usr
#   Root directory of MPI installation. <path>
#   SCAFESCC_MPI_ROOT="/usr";
#
# * Root directory of Boost installation. <path> default=/usr
#   SCAFESCC_BOOST_ROOT="/usr";
#
# * Library directory of Boost installation default=/usr/lib64.
#   SCAFESCC_BOOST_LIB="/usr/lib64"; <path>
#
# * Enable Boost.MPI? <no/yes>. default=yes
#   SCAFESCC_ENABLE_BOOST_MPI="yes";
#
# * Enable OpenMP? <no/yes>. default=no
#   SCAFESCC_ENABLE_OPENMP="no";
#
# * Enable NetCDF? <no/yes>. default=no
#   SCAFESCC_ENABLE_NETCDF="no";
#
# * Only needed if SCAFESCC_ENABLE_NETCDF="yes":
#   Root directory of NetCDF installation. <path>
#   SCAFESCC_NETCDF_ROOT="/usr";
#
# * Only needed if SCAFESCC_ENABLE_NETCDF="yes":
#   Root directory of Hdf5 installation. <path> default=/usr
#   SCAFESCC_HDF5_ROOT="/usr";
#
# * Enable ADOL-C? <no/yes>. default=no
#   SCAFESCC_ENABLE_ADOLC="no";
#
# * Only needed if SCAFESCC_ENABLE_ADOLC="yes": default=/usr
#   Root directory of ADOL-C installation. <path>
#   SCAFESCC_ADOLC_ROOT="/usr";
#
# * Include directory of ScaFES? default=/usr/include
#   SCAFESCC_SCAFES_INC="/usr/include";
#
# * Library directory of ScaFES? default=/usr/lib
#   SCAFESCC_SCAFES_LIB="/usr/lib";
#
# ScaFES
# Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
# For details, see the files COPYING and LICENSE in the base directory
# of the package.

(
#------------------------------------------------------------------------------#
if [ "$SCAFESCC_CXX" == "" ] ; then
    SCAFESCC_CXX="g++";
    echo "* WARNING: Use default value SCAFESCC_CXX=$SCAFESCC_CXX."
fi

#------------------------------------------------------------------------------#
if [ "$SCAFESCC_CXXFLAGS" == "" ] ; then
    SCAFESCC_CXXFLAGS="-DNDEBUG -O3 -Wall";
    echo "* WARNING: Use default value SCAFESCC_CXXFLAGS=$SCAFESCC_CXXFLAGS."
fi

#------------------------------------------------------------------------------#
if [ "$SCAFESCC_ENABLE_MPI" == "" ] ; then
    SCAFESCC_ENABLE_MPI="yes";
    echo "* WARNING: Use default value SCAFESCC_ENABLE_MPI=$SCAFESCC_ENABLE_MPI."
fi

#------------------------------------------------------------------------------#
if [ "$SCAFESCC_ENABLE_BOOST_MPI" == "" ] ; then
    SCAFESCC_ENABLE_BOOST_MPI="yes";
    echo "* WARNING: Use default value SCAFESCC_ENABLE_BOOST_MPI=$SCAFESCC_ENABLE_BOOST_MPI."
fi

#------------------------------------------------------------------------------#
if [ "$SCAFESCC_ENABLE_ADOLC" == "" ] ; then
    SCAFESCC_ENABLE_ADOLC="no";
    echo "* WARNING: Use default value SCAFESCC_ENABLE_ADOLC=$SCAFESCC_ENABLE_ADOLC."
fi

#------------------------------------------------------------------------------#
if [ "$SCAFESCC_ENABLE_OPENMP" == "" ] ; then
    SCAFESCC_ENABLE_OPENMP="no";
    echo "* WARNING: Use default value SCAFESCC_ENABLE_OPENMP=$SCAFESCC_ENABLE_OPENMP."
fi

#------------------------------------------------------------------------------#
if [ "$SCAFESCC_ENABLE_NETCDF" == "" ] ; then
    SCAFESCC_ENABLE_NETCDF="no";
    echo "* WARNING: Use default value SCAFESCC_ENABLE_NETCDF=$SCAFESCC_ENABLE_NETCDF."
fi

#------------------------------------------------------------------------------#
if [ "$SCAFESCC_BOOST_ROOT" == "" ] ; then
    SCAFESCC_BOOST_ROOT="/usr";
    echo "* WARNING: Use default value SCAFESCC_BOOST_ROOT=$SCAFESCC_BOOST_ROOT."
fi

#------------------------------------------------------------------------------#
if [ "$SCAFESCC_BOOST_LIB" == "" ] ; then
    SCAFESCC_BOOST_LIB="/usr/lib64";
    echo "* WARNING: Use default value SCAFESCC_BOOST_LIB=$SCAFESCC_BOOST_LIB."
fi

#------------------------------------------------------------------------------#
if [ "$SCAFESCC_ADOLC_ROOT" == "" ] ; then
    SCAFESCC_ADOLC_ROOT="/usr";
    echo "* WARNING: Use default value SCAFESCC_ADOLC_ROOT=$SCAFESCC_ADOLC_ROOT."
fi

#------------------------------------------------------------------------------#
if [ "$SCAFESCC_NETCDF_ROOT" == "" ] ; then
    SCAFESCC_NETCDF_ROOT="/usr";
    echo "* WARNING: Use default value SCAFESCC_NETCDF_ROOT=$SCAFESCC_NETCDF_ROOT."
fi

#------------------------------------------------------------------------------#
if [ "$SCAFESCC_HDF5_ROOT" == "" ] ; then
    SCAFESCC_HDF5_ROOT="/usr";
    echo "* WARNING: Use default value SCAFESCC_HDF5_ROOT=$SCAFESCC_HDF5_ROOT."
fi

#------------------------------------------------------------------------------#
if [ "$SCAFESCC_SCAFES_INC" == "" ] ; then
    SCAFESCC_SCAFES_INC="/usr/include";
    echo "* WARNING: Use default value SCAFESCC_SCAFES_INC=$SCAFESCC_SCAFES_INC."
fi

#------------------------------------------------------------------------------#
if [ "$SCAFESCC_SCAFES_LIB" == "" ] ; then
    SCAFESCC_SCAFES_LIB="/usr/lib";
    echo "* WARNING: Use default value SCAFESCC_SCAFES_LIB=$SCAFESCC_SCAFES_LIB."
fi

# From here on, all variables should be defined.
set -u
set -e

#------------------------------------------------------------------------------#
declare local SCAFES_INC="$SCAFESCC_SCAFES_INC"
declare local SCAFES_LIB="$SCAFESCC_SCAFES_LIB"

declare local BOOST_INC="$SCAFESCC_BOOST_ROOT/include"
declare local BOOST_LIB="$SCAFESCC_BOOST_LIB"
declare local ADOLC_INC="$SCAFESCC_ADOLC_ROOT/include"
declare local ADOLC_LIB="$SCAFESCC_ADOLC_ROOT/lib"
declare local NETCDF_INC="$SCAFESCC_NETCDF_ROOT/include"
declare local NETCDF_LIB="$SCAFESCC_NETCDF_ROOT/lib"

#------------------------------------------------------------------------------#
declare local SCAFES_INC_DIR="-I${SCAFES_INC}"
declare local SCAFES_LDFLAGS="-L${SCAFES_LIB}"
# [KF, 06/01/2015]
# ScaFES is header-only (if NetCDF is NOT enabled).
declare local SCAFES_LIBADD=""

declare local BOOST_INC_DIR="-I${BOOST_INC}"
declare local BOOST_LDFLAGS="-L${BOOST_LIB}"
declare local BOOST_LIBADD="";
if [ ${SCAFESCC_ENABLE_BOOST} == "yes" ] ; then
    BOOST_LIBADD="${BOOST_LIBADD} -lboost_regex"
    BOOST_LIBADD="${BOOST_LIBADD} -lboost_program_options"
fi
if [ ${SCAFESCC_ENABLE_BOOST_SERIALIZATION} == "yes" ] ; then
    BOOST_LIBADD="${BOOST_LIBADD} -lboost_serialization"
fi
if [ ${SCAFESCC_ENABLE_BOOST_MPI} == "yes" ] ; then
    BOOST_LIBADD="${BOOST_LIBADD} -lboost_mpi"
fi

declare local ADOLC_INC_DIR=""
declare local ADOLC_LDFLAGS=""
declare local ADOLC_LIBADD=""
if [ ${SCAFESCC_ENABLE_ADOLC} == "yes" ] ; then
    ADOLC_INC_DIR="-I${ADOLC_INC}"
    ADOLC_LDFLAGS="-L${ADOLC_LIB}"
    ADOLC_LIBADD="-ladolc"
fi

declare local NETCDF_INC_DIR=""
declare local NETCDF_LDFLAGS=""
declare local NETCDF_LIBADD=""
if [ ${SCAFESCC_ENABLE_NETCDF} == "yes" ] ; then
   NETCDF_INC_DIR="-I${NETCDF_INC}"
   NETCDF_LDFLAGS="-L${NETCDF_LIB}"
   NETCDF_LIBADD="-lnetcdf"
fi

#------------------------------------------------------------------------------#
declare local CPPFLAGS="${SCAFES_INC_DIR}"
# [KF, 06/01/2015]
# Dirty hack at Juqueen regarding the compiler flag for the IBM compiler.
#CPPFLAGS="${CPPFLAGS} -qlanglvl=extended0x "
CPPFLAGS="${CPPFLAGS} -std=c++11"
declare local LDFLAGS="${SCAFES_LDFLAGS}"
declare local LIBADD="${SCAFES_LIBADD}"

#------------------------------------------------------------------------------#
if [ ${SCAFESCC_ENABLE_BOOST} == "yes" ] ; then
    CPPFLAGS="${CPPFLAGS} ${BOOST_INC_DIR}"
    LDFLAGS="${LDFLAGS} ${BOOST_LDFLAGS}"
    LIBADD="${LIBADD} ${BOOST_LIBADD}"
fi

#------------------------------------------------------------------------------#
if [ ${SCAFESCC_ENABLE_OPENMP} == "yes" ] ; then
    CPPFLAGS="${CPPFLAGS} -fopenmp"
    LDFLAGS="${LDFLAGS} -fopenmp"
fi

#------------------------------------------------------------------------------#
if [ ${SCAFESCC_ENABLE_ADOLC} == "yes" ] ; then
   CPPFLAGS="${CPPFLAGS} ${ADOLC_INC_DIR}"
   LDFLAGS="${LDFLAGS} ${ADOLC_LDFLAGS}"
   LIBADD="${LIBADD} ${ADOLC_LIBADD}"
fi

#------------------------------------------------------------------------------#
if [ ${SCAFESCC_ENABLE_NETCDF} == "yes" ] ; then
   CPPFLAGS="${CPPFLAGS} ${NETCDF_INC_DIR}"
   LDFLAGS="${LDFLAGS} ${NETCDF_LDFLAGS}"
   LIBADD="${LIBADD} ${NETCDF_LIBADD}"
fi

#------------------------------------------------------------------------------#
declare local MYCXX="${SCAFESCC_CXX}"
if [ ${SCAFESCC_ENABLE_MPI} == "yes" ] ; then
   MYCXX="mpic++";
fi
if [ ${SCAFESCC_ENABLE_BOOST_MPI} == "yes" ] ; then
   MYCXX="mpic++";
fi

# [KF, 06/01/2015]
# Dirty hack which solves the frontend / backend issue at Juqueen.
#MYCXX="mpig++"; # Enable this line if GCC is used.
#MYCXX="mpiclang++11"; # Enable this line if GCC is used.

#------------------------------------------------------------------------------#
execString="${CPPFLAGS} ${SCAFESCC_CXXFLAGS} ${LDFLAGS} $@ $LIBADD";
echo "+++ ${MYCXX} ${execString} +++"
${MYCXX} ${execString}
)
