/* ScaFES
 * Copyright (c) 2011-2018, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES.hpp
 *
 *  @brief Main include file for applications using ScaFES.
 *
 * This file has to be included if a new initial boundary value problem
 * should be created.
 */

/**
 *  \mainpage ScaFES Documentation
 *  \tableofcontents
 *
 * <BLOCKQUOTE>We present ScaFES, an open-source HPC framework written in C++11 for
 * solving initial boundary value problems using explicit numerical methods
 * in time on structured grids. It is designed to be highly-scalable and
 * very user-friendly, i.e. to exploit all levels of parallelism and provide
 * easy-to-use interfaces. Besides, the numerical nomenclature is reflected
 * in a nearly one-to-one mapping.</BLOCKQUOTE>
 *
 * [Source: Martin Flehmig, Kim Feldhoff, and Ulf Markwardt:
 * "ScaFES: An Open-Source Framework for Explicit Solvers Combing
 * High-Scalability with User-Friendliness",
 * ARCS 2014 - 27th International Conference of Computing Systems,
 * Workshop-Proceedings, VDE-Verlag, 2014]
 *
 * The framework is available under an open-source license. For details,
 * please refer to the files COPYING and LICENSE in the base directory
 * of the package.
 *
 * \section example Examples
 * In order to show how a time-dependent PDE problem can be implemented
 * in ScaFES, we have chosen the following examples:
 *
 *  - \ref heatEqnFDM "Discretized heat equation problem on unit hybercube"
 *  - \ref waveEqnFDM "Discretized wave equation problem on unit hybercube (Under construction [KF], 2016-08-29)"
 *  - \ref lameNavierEqnFDM "Discretized Lame Navier problem on unit hybercube (Under construction [KF], 2016-08-29)"
 *  - \ref lameNavierEqnFDMFirstOrder "Discretized Lame Navier problem on unit hybercube rewritten as a first order system wrt. time (Under construction [KF], 2016-08-31)"
 *  - \ref EMFDTDng3D "Discretized EM (electromagnetic) problem on an ellipsoid"
 *
 * The examples can be setup and executed via bash scripts
 * located in the corresponding directories.
 * The bash scripts are based on the shell script \c scafesrun.sh
 * which can be configured via environment variables prefixed by \c SCAFESRUN.
 * First example can be executed as follows:
 *
        cd examples/HeatEqnFDM
        ./RUN.sh
 *
 * \section Implementation How to Implement a New Problem?
 * The usual way to implement a new discretised problem is as follows:
 * - Create the skeleton of the discretized problem.
 * The skeleton of the discretised problem can be setup by creating a new
 * class which is derived from the class \c ScaFES::Problem.
 * - Add all fields to the problem.
 * All fields can be added to the discrete problem via the constructor
 * of the base class \c ScaFES::Problem.
 * - Evaluate all fields at all grid nodes.
 * All known fields can be evaluated at all global inner resp. border grid nodes
 * at all time steps by implementing the methods \c evalInner resp.
 * \c evalBorder.
 * - Initialise all unknown fields at all grid nodes.
 * All unknown fields can be initialised at all global inner resp. border
 * grid nodes at the starting time by implementing the methods \c initInner
 * resp. \c initBorder.
 * -  Update all unknown fields at all grid nodes.
 * All unknown fields can be updated at all global inner resp. border
 * grid nodes at the starting time by implementing the methods \c updateInner
 * resp. \c updateBorder.
 * - Write the main program.
 *
 * \section contact Contact
 * Bugs can be reported to: scafes-commits@fusionforge.zih.tu-dresden.de
 *
 * More information about ScaFES like new releases can be found at the ScaFES web page:
 * <a href="http://www.tu-dresden.de/zih/scafes">http://www.tu-dresden.de/zih/scafes</a>
 *
 * \section funding Funding
 * ScaFES was partially supported by the Federal Ministry of Education and Research (BMBF)
 * within the project HPC-FLiS under the support code 01 IH 11 009.
 *
 ******************************************************************************/

/**
 *  \page FAQ FAQ
 *  Frequently Asked Questions
 *  \tableofcontents
 *
 * \section Required_software How can the required software by installed on current Linux distributions?
 *
 * \subsection reqsoft_opensuse Linux OpenSUSE 12.3
 *
 * Start "yast2" --> "Software Management" and install the following items:

         gcc
         gcc-c++
         make
         libboost-dev
         libboost-all-dev
         libboost-progam-options
         libboost-regex
         libboost-serialization
         libboost-iostreams
         libboost-mpi
         openmpi
         openmpi-devel
         libnetcdf-devel
         libnetcdf
         netcdf
         hdf5
         libhdf5-0
         libhdf5-hl0
         bash

 * \subsection reqsoft_ubuntu Linux Ubuntu 12.04
 * Start "Ubuntu Software Center" and install the following items:

        gcc
        gcc-c++
        make
        libboost-dev
        libboost-all-dev
        libboost-progam-options  >= 1.42.0
        libboost-regex  >= 1.42.0
        libboost-serialization >= 1.42.0
        libboost-iostreams >= 1.42.0
        libboost-mpi >= 1.42.0
        openmpi-bin
        openmpi-dev
        libnetcdf-dev
        libnetcdf
        hdf5
        libhdf5
        libhdf5-0
        libhdf5-hl0
        bash

 *
 *\remark If the packages will be installed from the Ubuntu Software Center,
 * one has to ensure that the buttom "Show technical items" located at the
 * bottom of the Ubuntu Software Center in the status bar has been activated.
 * Otherwise, not all packages which can be installed will be seen.
 *
 * \section macos_installation How to install ScaFES on macOS Sierra 10.12.5:
 * Install scafes-trunk on macOS Sierra 10.12.5 in combination with homebrew (https://brew.sh/).
 * Homebrew can be used to install missing programs, packages, libraries, ... .
 * Homebrew does use the system's default compiler, i.e. Apple LLVM version 8.1.0 (clang-802.0.42) for macOS Sierra 10.12.5.
 * Therefore ScaFES should be build with the same compiler.
 * An newer version of LLVM/clang++ compiler can be installed with homebrew if OpenMP is desired.
 * List of programs that need to be installed with homebrew:

        brew install llvm
        brew install boost@1.57 --c++11 --with-mpi --without-single
        brew install homebrew/science/adol-c
        brew install homebrew/science/hdf5
        brew install homebrew/science/netcdf
        brew install homebrew/science/parallel-netcdf
        brew install automake
        brew install doxygen
        brew install graphviz

 * this will install the following programs by dependency:

        open-mpi
        autoconf

 * call bootstrap script:

        ./bootstrap.sh

 * call configure script (should be called from a build folder with ./../configure):

        export CC='/usr/local/opt/llvm/bin/clang'
        export CFLAGS='-O3 -Wall'
        export CXX='/usr/local/opt/llvm/bin/clang++'
        export CXXFLAGS='-DNDEBUG -Wall -Wextra -Wno-unused -O3'
        export LDFLAGS='-L/usr/local/opt/llvm/lib'
        export CPPFLAGS='-I/usr/local/opt/llvm/include'
        ./configure \
        --enable-dist="no" \
        --with-mpi='/usr/local/opt/openmpi' \
        --enable-mpi="yes" \
        --enable-openmp="yes" \
        --with-boost='/usr/local/opt/boost@1.57' \
        --with-boost-libdir='/usr/local/opt/boost@1.57/lib' \
        --with-boost-serialization \
        --enable-boost-mpi="yes" \
        --with-netcdf='/usr/local/opt/parallel-netcdf' \
        --with-hdf5='/usr/local/opt/hdf5/' \
        --enable-netcdf="yes" \
        --enable-adolc="yes" \
        --with-adolc='/usr/local/opt/adol-c' \
        --enable-gmock="yes" \
        --with-gmock='/usr/local/opt/gmock' \
        --enable-doc="yes" \
        --prefix=<installdir>

 * call make and make install

        make
        make install

 * See also INSTALL for more information to build ScaFES.

 * \section Netcdf_file_correct How can be checked if the data fields were written correctly to the NetCDF data file?

         ncdump <datafile.nc>
         ncdump -h <datafile.nc>

 *
 * \section ScaFES_VT Application will end up with an run time error when using VampirTrace
 * Under OpenSuse 12.2 and VampirTrace, the following run time error occur:

         VampirTrace: FATAL: Maximum retries (10) for gethostid exceeded!

* Looking in the VampirTrace library, VampirTrace exists with a run time error
* in case of a host id equal to 0.
* For Linux OpenSUSE, the host identifier got by gethostid() is zero.
* Thus, go to the sources of VampirTrace (File vtlib/vt_pform_linux.c : 229ff)
* and replace call to gethostid() by an arbitrary integer different from zero:

        get unique numeric SMP-node identifier
        hostid_retries = 0;
        while( !vt_node_id && (hostid_retries++ < VT_MAX_GETHOSTID_RETRIES) ) {
            vt_node_id = 48; // gethostid();
        }
        if (!vt_node_id)
            vt_error_msg("Maximum retries (%i) for gethostid exceeded!",
                 VT_MAX_GETHOSTID_RETRIES);
        }

* Afterwards, compile and install VampirTrace again using make and make install.
*
*
* \section ADOL-C_and_OpenMP Compile error: Undefined reference to ADOLC_OpenMP_Handler etc.
* Case: ADOL-C and OpenMP are activated at the same time.
* Some or all of the following error messages occur during compiling:

        undefined reference to `ADOLC_parallel_doCopy'
        undefined reference to `beginParallel()'
        undefined reference to `endParallel()'

* The reason is: ADOL-C was not built with OpenMP support.
* Thus, please rebuild ADOL-C from source and enable OpenMP support.
* Add \code --enable-openmp\endcode to the configure call of ADOL-C.


\section scafes_2_1_optional_modules ScaFES 2.0.1: ptional modules like ADOL-C cannot be arbitrarily enabled / disabled.
* Some or all of the following error messages occur during compiling:

        undefined reference to `adouble::operator=(adub const&)'
        undefined reference to `adub::~adub()'
        undefined reference to `operator*(double, badouble const&)'

* scafes/ScaFES_Config.hpp was shipped with the release 2.0.1.
* \code build/scafes/ScaFES_Config.hpp\endcode will be created dynamically
* during configure.
* Header files will be included in example/HeatEqnFDM in the following order:
    \code -I$(top_builddir)/src -I$(top_srcdir)\endcode
* This is a bug which has been fixed in ScaFES version 2.1.0 by
* removing \code scafes/ScaFES_Config.hpp\endcode from the "make dist" rule.


\section scafes_2_0_1_wunused  ScaFES 2.0.1 Warnings related to unused parameters
* Some or all of the following warning messages occur during compiling:

        Warnung: unbenutzter Parameter »dfGradDep« [-Wunused-parameter]

* Names of variables have not been commented out in parameters list.
* This is harmless, especially no user action is required.
* This has been fixed in ScaFES version 2.1.0.


\section libhdf5_configure The library hdf5 cannot be found during configure run.
* Using Linux Mint, libhdf5 cannot be found during configure run
* although the correct path was given to the configure call:

        ./configure --with-hdf5-libdir=/usr/lib/x86_64-linux-gnu ...

* The package HDF5-dev has been installed via package manager,
* but the hdf5 libray is not named libhdf5.so
* Thus, set a soft link libhdf5.so to target libhdf.so.XYZ:

        cd /usr/lib/x86_64-linux-gnu/
        sudo ln -s libhdf5.so.7 libhdf5.so


\section Makefiles_not_create Makefiles in folders tests/system and tests/unit will not be created during configure call.

* ScaFES trunk has been prepared for development mode either
* - via calling the script

      ./scripts/scafespreparebuildsystem.sh no

* - or the development mode was already set in trunk but configure sets --enable-dist=yes per default.
* Thus, Makefiles will not be written (see section AC_CONFIG_FILES in configure.ac)

* Preparations of development / distribution mode must be consistent.
* Thus, either call

      BTS_ssetupsystem.sh $MACHINENAME no

* or call

      BTS_preparebuildsystem.sh no
      ./bootstrap.sh
      ./configure --enable-dist=no ...


* \section netcdf_nodatafile NetCDF is activated via "--enable-netcdf" but no NetCDF data file has been written after running the simulation.

* Write options in std::vector<ScaFES::WriteHowOften> writeToFile
* at ScaFES constructor for all data fields of the example
* are chosen as "NEVER":

          HeatEqnFDM<double, DIM> ppp(paramsCl, gg, false, nameDatafield, stencilWidth,
                                  isKnownDf, nLayers, defaultValue, writeToFile,
                                  computeError, geomparamsInit);

* Change the value of a least one data file within the parameter "writeToFile" to

        a) LIKE_GIVEN_AT_CL
        b) ALWAYS
        c) AT_START
        d) AT_END
        e) AT_START_AND_END


\section netcdf_noatalltimesteps NetCDF is activated via "--enable-netcdf", NetCDF data file has been written, but the data fields have not been written to the file at ALL time steps.
* Write options in std::vector<ScaFES::WriteHowOften> writeToFile
* at ScaFES constructor for all data fields of the example
* are not chosen as "ALWAYS" or LIKE_GIVEN_AT_CL

          HeatEqnFDM<double, DIM> ppp(paramsCl, gg, false, nameDatafield, stencilWidth,
                                  isKnownDf, nLayers, defaultValue, writeToFile,
                                  computeError, geomparamsInit);

* Change the value of data fields within the parameter "writeToFile"
* either to "ALWAYS" or to "LIKE_GIVEN_AT_CL".


\section netcdf_like_given_at NetCDF is activated via "--enable-netcdf", NetCDF data file has been written, write options in std::vector<ScaFES::WriteHowOften> writeToFile are chosen as LIKE_GIVEN_AT_CL but the data fields have not been written to the file at ALL time steps.
* Write option LIKE_GIVEN_AT_CL means that the number of snapshots must
* be given via the command line.
* Thus:
* - a) export SCAFESRUN_N_SNAPSHOTS=val and call "scafesrun.sh" afterwards
* - b) or add --nSnapshots=val as command line option to the executable
* in order to run the simulation/program.


* \section netcdf_serial NetCDF data file has been written (.nc) but VisIt shows a incorrect picture.
* Serial I/O has been used with parallel computation.
* Thus, only MPI rank 0 will write values of the data fields to the NetCDF
* data file. Hence, all other values will be equal to the default value.
* The visualization should be look similar to the one in
* ScaFES_VisIt_Serial_IO.png.
*
 * Thus,
 * - either, run simulation in serial, too
 * - or install NetCDF with parallel support: ScaFES support NetCDF versions >= 4.3 for this purpose.
 *
 *
 * \section visit_defaultvalues VisIt shows values in the range of 9.9e+36.
 * This is the default value of VisIt,
 * i.e., no value was written to the NetCDF data file.
 * You can check this with the command: ncdump OUTPUTFILE.nc
 * The output should like similar to the following line:

        , -, -, -,

 * Check other pit falls wrt. to the visualization.
 */

/**
 * \page WriteNetCDF_checklist WriteNetCDF_checklist
 * Check list for writing NetCDF data fields
 *  \tableofcontents
 *
 - Set parameter of constructor appropriate.
    * Determine which data field will be written to output file.
    * Determine at which time step the data field should be written.
    In the following example three data fields has been declared.
    The data fields F and G will never be written to output file.
    The data field U will be outputted as given on command line:

          std::vector<std::string> nameDatafield(3);
          nameDatafield[0] = "F";
          nameDatafield[1] = "G";
          nameDatafield[2] = "U";
          std::vector<ScaFES::WriteHowOften> writeToFile(3, ScaFES::WriteHowOften::LIKE_GIVEN_AT_CL);
          std::vector<bool> computeError(3);
          computeError[0] = false;
          computeError[1] = false;
          computeError[2] = true;
          HeatEqnFDM<double, DIM> ppp(paramsCl, gg, false, nameDatafield, stencilWidth,
                                  isKnownDf, nLayers, defaultValue, writeToFile,
                                  computeError, geomparamsInit);
          ppp.iterateOverTime();

    Possible values for ScaFES::WriteHowOften are:

        LIKE_GIVEN_AT_CL (configured via command line)
        NEVER
        ALWAYS
        AT_START (initialization)
        AT_END (finialization)
        AT_START_AND_END (initialization and finialization)

 - If the option <tt>LIKE_GIVEN_AT_CL</tt> was chosen, either set
     a) <tt>export SCAFESRUN_N_SNAPSHOTS=val</tt> and call <tt>scafesrun.sh</tt> afterwards
     b) or add <tt>--nSnapshots=val</tt> as command line option to the executable
    in order to run the simulation / program.
 - Case: parallel I/O?
    Is NetCDF version >= 4.3 with parallel support installed?
*/


/**
 * \page display_df_using_Visit displayDFUsingVisit
 * Displaying data fields using VisIt
 *  \tableofcontents
 *
 * Start the application (using 4 processes and without displaying the
 *  start up screen):

       visit -np 4 -o <file.nc> -nosplash

* Set up new variable "namePsi" created from variables "fnew_0_im"
  and "fnew_0_re" given in the NetCDF file.

       Controls->Expression->New:
          Name = normPsi
       Type = Scalar Mesh Variable
           Standard Editor = sqr(fnew_0_re)+sqr(fnew_0_im)
       ->Apply, ->Post

* Draw a variable:

      * Add->Pseudocolor->u_new (name of the variable).
      * Draw

* Select only a part of the whole grid:

      * Operator->Selection->Box:
           X-Minimum=119
           X-Maximum=124
           Y-Minimum=0
           Y-Maximum=389
           Z-Minimum=0
           Z-Maximum=100

* Show data in a table:

     * Add->Spreadsheet->epsNew
     * Draw
     * Check the box "Color" (in the right upper corner).
 *
 */


#ifndef SCAFES_HPP_
#define SCAFES_HPP_

#include "ScaFES_Config.hpp"

#include "ScaFES_Problem.hpp"

#endif
