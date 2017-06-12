#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
if AMHAVEVT
VTxxxCPPFLAGS = @VT_CPPFLAGS@
endif
endif

#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
if AMHAVEVT
if AMHAVEOMP
VTOPARIxxxCPPFLAGS = @VT_CPPFLAGS@ \
-vt:opari "-rcfile $(top_builddir)/src/scafes.opari.rc -table $(top_builddir)/src/scafes.opari.tab.c"
endif
endif
endif

#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
if AMHAVECUDA
CUDAxxxCOMPILER = g++
CUDAxxxCPPFLAGS = -I$(top_srcdir)/src/boost_fixed_headers
CUDAxxxLDFLAGS  = @CUDA_LDFLAGS@
# CUDAxxxCOMPILER = /sw/global/compilers/intel/2013/bin/icpc
# CUDAxxxCPPFLAGS = -I$(top_srcdir)/src/boost_fixed_headers -DJOKER
# CUDAxxxLDFLAGS  = @CUDA_LDFLAGS@  \
#                  -L/sw/global/compilers/intel/2013/composer_xe_2013/lib/intel64 \
#                  -lirc -limf
endif
endif

################################################################################
### Library.
SCAFESSRCDIR = $(top_srcdir)/src
SCAFESBUILDDIR = $(top_builddir)/src
COMMONSOURCES = ScaFES_Dummy.cpp

if !AMDISTRIBUTESOFTWARE
if AMHAVECUDA
CUDASOURCES = ScaFES_ProblemCuda.cu
endif
endif

SOURCESxxxCOMMON = $(COMMONSOURCES)
if !AMDISTRIBUTESOFTWARE
if AMHAVECUDA
SOURCESxxxCOMMON += $(CUDASOURCES)
endif
endif

#------------------------------------------------------------------------------#
SCAFESINCDIR = $(SCAFESSRCDIR)
include_HEADERS = \
   ScaFES.hpp \
   ScaFES_Ntuple.hpp \
   ScaFES_Ntuple_FreeFunctions.hpp \
   ScaFES_Complex.hpp \
   ScaFES_Complex_FreeFunctions.hpp \
   ScaFES_Parameters.hpp \
   ScaFES_Grid.hpp \
   ScaFES_GridSub.hpp \
   ScaFES_GridGlobal.hpp \
   ScaFES_Communicator.hpp \
   ScaFES_Environment.hpp \
   ScaFES_Timer.hpp \
   ScaFES_Buffer.hpp \
   ScaFES_BufferApi.hpp \
   ScaFES_BufferDummy.hpp \
   ScaFES_DataFile.hpp \
   ScaFES_DataFileApi.hpp \
   ScaFES_DataFileDummy.hpp \
   ScaFES_DataField.hpp \
   ScaFES_Problem.hpp

# Header file will be generated during run of configure.
# Header file must be installed.
# Header file must not be distributed as it contains which modules were
# enabled / disabled for the specific configure run.
# http://www.gnu.org/software/automake/manual/html_node/Headers.html
# Last accessed: 20.11.2014, [KF]
nodist_include_HEADERS=$(SCAFESBUILDDIR)/ScaFES_Config.hpp

if !AMDISTRIBUTESOFTWARE
include_HEADERS += ScaFES_ProblemCuda.hpp
endif

#------------------------------------------------------------------------------#
CPPFLAGSxxxCOMMON = \
            -I$(SCAFESSRCDIR) \
            -I$(SCAFESBUILDDIR)

#------------------------------------------------------------------------------#
LDFLAGSxxxCOMMON =

if !AMDISTRIBUTESOFTWARE
if AMHAVECUDA
LDFLAGSxxxCOMMON += $(CUDAxxxLDFLAGS)
endif
endif

#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
if AMHAVECUDA
CUDAxxxCOMPILERxxxOPTIONS= \
            $(CUDAxxxCPPFLAGS) \
            $(CFLAGS) \
            $(CPPFLAGSxxxCOMMON) \
            $(DEFAULT_INCLUDES) \
            $(INCLUDES) \
            $(AM_CPPFLAGS) \
            $(CPPFLAGS)
endif
endif

if !AMDISTRIBUTESOFTWARE
if AMHAVECUDA
.cu.o:
	$(CUDA_CC) -ccbin=$(CUDAxxxCOMPILER) \
                   -arch=sm_20 --x cu --keep --verbose \
          --compiler-options="$(CUDAxxxCOMPILERxxxOPTIONS)" -c $<
.cu.lo:
	$(top_srcdir)/cudalt.py $@ $(CUDA_CC) -ccbin=$(CUDAxxxCOMPILER) \
           -arch=sm_20 --x cu --keep --verbose \
           --compiler-options=\"$(CUDAxxxCOMPILERxxxOPTIONS)\" -c $<
endif
endif

#------------------------------------------------------------------------------#
# Libtool is used to build the library.
# Shared libraries must be installed.
lib_LTLIBRARIES       = libscafes.la
libscafes_la_SOURCES  = $(COMMONSOURCES)
if !AMDISTRIBUTESOFTWARE
libscafes_la_SOURCES  += $(CUDASOURCES)
endif
libscafes_la_CPPFLAGS = $(CPPFLAGSxxxCOMMON)
libscafes_la_LDFLAGS = -version-info @SCAFES_VERSION_LIBRARY@
libscafes_la_LDFLAGS  += $(VTxxxCPPFLAGS)
                         $(VTOPARIxxxCPPFLAGS)
libscafes_la_LIBADD   = $(LDFLAGSxxxCOMMON)


################################################################################

if !AMDISTRIBUTESOFTWARE
if AMHAVECUDA
EXTRA_DIST = boost_fixed_headers
endif
endif

#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
clean-local:
	-rm -f *.gcno *.opari.inc *.opari.rc *.pomp.cpp
endif
