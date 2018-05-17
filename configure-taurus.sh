(
set -u
set -e

################################################################################
# Must be edited by the user!
declare local nameSystem="ZIH_TAURUS_GCC_5_3_0"
declare -r local versionScafes="masterthesis"
declare local installdirScafes="${SCRATCH}/scafes-${versionScafes}-${nameSystem}"
################################################################################

#------------------------------------------------------------------------------#
case "$nameSystem" in
 ############################################################################
 ZIH_TAURUS_GCC_5_3_0 )
     #----------------------------------------------------------------------#
     module purge
     module load autotools
     #----------------------------------------------------------------------#
     export SCAFESCC_ENABLE_MPI="yes"
     export SCAFESCC_ENABLE_OPENMP="yes"
     export SCAFESCC_ENABLE_GMOCK="yes"
     export SCAFESCC_ENABLE_DEBUG_MODE="no"
     export SCAFESCC_ENABLE_CUDA="no";
     export SCAFESCC_ENABLE_BOOST_MPI="yes"
     export SCAFESCC_ENABLE_NETCDF="yes";
     export SCAFESCC_ENABLE_ADOLC="no"
     export SCAFESCC_ENABLE_VT="no";
     export SCAFESCC_ENABLE_DOXYGEN="no" # Doxygen 1.8.11 [available, 2017-06-07, MS]
     export SCAFESCC_ENABLE_GCOV="no"
     #----------------------------------------------------------------------#
     module load boost/1.58.0-gnu5.1
     #boost version 1.58.0-gnu5.1 for x86_64 architecture loaded.
     #gcc version 5.1.0 for x86_64 architecture loaded.
     #bullxmpi version 1.2.8.4 for x86_64 architecture loaded.
     #NOTE: please ignore the MXM warnings about conflicting cpu frequencies when using bullxmpi.
     #python version 2.7 for x86_64 architecture loaded.
     export SCAFESCC_BOOST_ROOT="/sw/taurus/libraries/boost/1.58.0-gnu5.1" # 1.58.0
     export SCAFESCC_BOOST_LIB="$SCAFESCC_BOOST_ROOT/lib";
     #----------------------------------------------------------------------#
     module load bullxmpi/1.2.8.4
     export SCAFESCC_MPI_ROOT="/sw/taurus/libraries/bullxmpi/1.2.8.4" # 1.2.8.4
     #----------------------------------------------------------------------#
     module load netcdf/4.3.3.1-gcc-5.2.0 # [available 28.11.2017, MS]
     # Module netcdf/4.3.3.1-gcc-5.2.0 and 2 dependencies loaded.
     # hdf5/1.8.15-gcc-5.1.0
     # gcc/5.2.0
     export SCAFESCC_NETCDF_ROOT="/sw/taurus/libraries/netcdf/4.3.3.1-gcc-5.2.0";
     export SCAFESCC_HDF5_ROOT="/sw/taurus/libraries/hdf5/1.8.15-gcc-5.10";
     #----------------------------------------------------------------------#
     module load adolc/2.5.0 # [available, 03/23/2015, KF]
     export SCAFESCC_ADOLC_ROOT="/sw/taurus/libraries/adolc/2.5.0"; # 2.5.0 With OpenMP supp.
     #----------------------------------------------------------------------#
     module load gmock/1.6.0 # [available, 14.10.2014, KF]
     export SCAFESCC_GMOCK_ROOT="/sw/taurus/tools/gmock/1.6.0"; # 1.6.0
     #----------------------------------------------------------------------#
     module load doxygen/1.8.11 # [available, 07.06.2017, MS]
     export SCAFESCC_DOXYGEN_ROOT="/sw/taurus/tools/doxygen/1.8.11"; # 1.8.11
     #----------------------------------------------------------------------#
     module load gcc/5.3.0 # [available, 2017-05-05, KF]
     export SCAFESCC_CXX='g++'
     export SCAFESCC_CXXFLAGS="-DNDEBUG -Wall -Wextra -Wno-unused -O3";
     #----------------------------------------------------------------------#
     module list
     ;;
 ############################################################################
 * )
     echo -e "\n\nError: Machine=$nameSystem not found.\n"
     ;;
esac

#------------------------------------------------------------------------------#
#wget https://tu-dresden.de/zih/forschung/ressourcen/dateien/abgeschlossene-#projekte/scafes/scafes-2.3.0.tar.gz
#tar -xzvf scafes-${versionScafes}.tar.gz
#cd scafes-${versionScafes}
declare -r local nameBuildDir="build-$nameSystem-$(date +%Y_%m_%d_%H_%M_%S)"
mkdir $nameBuildDir

./bootstrap.sh

cd $nameBuildDir

set -x

#------------------------------------------------------------------------------#
../configure \
 CXX=$SCAFESCC_CXX \
 CXXFLAGS="$SCAFESCC_CXXFLAGS" \
 --enable-dist=no \
 --with-mpi=$SCAFESCC_MPI_ROOT \
 --enable-mpi=$SCAFESCC_ENABLE_MPI \
 --enable-openmp=$SCAFESCC_ENABLE_OPENMP \
 --with-boost=$SCAFESCC_BOOST_ROOT \
 --with-boost-libdir=$SCAFESCC_BOOST_LIB \
 --with-boost-serialization \
 --enable-boost-mpi=$SCAFESCC_ENABLE_BOOST_MPI \
 --with-netcdf=$SCAFESCC_NETCDF_ROOT \
 --with-hdf5=$SCAFESCC_HDF5_ROOT \
 --enable-netcdf=$SCAFESCC_ENABLE_NETCDF \
 --enable-adolc=$SCAFESCC_ENABLE_ADOLC \
 --with-adolc=$SCAFESCC_ADOLC_ROOT \
 --enable-gmock=$SCAFESCC_ENABLE_GMOCK \
 --with-gmock=$SCAFESCC_GMOCK_ROOT \
 --enable-doc=$SCAFESCC_ENABLE_DOXYGEN \
 --prefix=$installdirScafes
set +x

#------------------------------------------------------------------------------#
make

#------------------------------------------------------------------------------#
make install

#------------------------------------------------------------------------------#
cd .. # Leave build folder.
)
