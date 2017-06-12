#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
if AMHAVEVT
VT_OPARI_CPPFLAGS = @VT_CPPFLAGS@ \
-vt:opari "-rcfile $(top_builddir)/src/scafes.opari.rc -table $(top_builddir)/src/scafes.opari.tab.c"
endif
endif

#------------------------------------------------------------------------------#
SCAFES_TESTS_DIR=$(top_srcdir)/tests/unit

CPPFLAGS_COMMON= \
        -I$(top_srcdir)/src \
        -I$(top_builddir)/src \
        @NETCDF_CPPFLAGS@ \
        @BOOST_CPPFLAGS@ \
        @MPI_CPPFLAGS@ \
        @OPENMP_CXXFLAGS@ \
        @GMOCK_CPPFLAGS@ \
        @ADOLC_CPPFLAGS@
if !AMDISTRIBUTESOFTWARE
CPPFLAGS_COMMON += $(VT_OPARI_CPPFLAGS)
endif
LDADD_COMMON =
if !AMDISTRIBUTESOFTWARE
LDADD_COMMON += $(VT_OPARI_CPPFLAGS)
endif
LDADD_COMMON += \
        $(top_builddir)/src/libscafes.la \
        @BOOST_LDFLAGS@ \
        @BOOST_PROGRAM_OPTIONS_LIB@ \
        @BOOST_REGEX_LIB@ \
        @BOOST_SERIALIZATION_LIB@ \
        @MPI_LDFLAGS@ \
        @BOOST_MPI_LIB@ \
        @OPENMP_CXXFLAGS@ \
        @NETCDF_LDFLAGS@ \
        @GMOCK_LDFLAGS@ \
        @ADOLC_LDFLAGS@
DEPENDENCIES_COMMON= $(top_builddir)/src/libscafes.la

################################################################################
include_HEADERS=
bin_PROGRAMS=
EXTRA_DIST=

################################################################################
SCAFESINCDIR        = $(top_srcdir)/src
SCAFESTESTSSRCDIR   = $(top_srcdir)/tests/unit
SCAFESTESTSINCDIR   = $(SCAFESTESTSSRCDIR)

### Test programs.
check_HEADERS=
check_PROGRAMS=
TESTS=

if AMHAVEGMOCK
TESTCPPFLAGSxxxCOMMON = -I$(SCAFESINCDIR) \
                        -I$(SCAFESTESTSINCDIR) \
                         $(CPPFLAGS_COMMON) \
                         $(GMOCKxxxCPPFLAGS)
TESTLDFLAGSxxxCOMMON = $(LDADD_COMMON) \
                       $(VTOPARIxxxCPPFLAGS) \
                       $(GMOCKxxxLDFLAGS)

#------------------------------------------------------------------------------#
# This is currently deactivated as command line options (argc, argv) have to be
# passed to the constructor of the class Parameters. [KF], 2015-11-17.
# # check_PROGRAMS += ParametersTest
# # ParametersTest_SOURCES = ScaFES_ParametersTest.hpp \
# #                          ScaFES_ParametersTest.cpp \
# #                          ScaFES_mainTest.cpp
# # ParametersTest_CPPFLAGS = $(TESTCPPFLAGSxxxCOMMON)
# # ParametersTest_LDFLAGS  = $(TESTLDFLAGSxxxCOMMON)
# # TESTS += ParametersTest

#------------------------------------------------------------------------------#
check_PROGRAMS += NtupleTest
NtupleTest_SOURCES  = ScaFES_NtupleMainTest.cpp \
                      ScaFES_VectorSpaceTest.hpp
NtupleTest_CPPFLAGS = $(TESTCPPFLAGSxxxCOMMON)
NtupleTest_LDFLAGS  = $(TESTLDFLAGSxxxCOMMON)
TESTS += NtupleTest

#------------------------------------------------------------------------------#
check_PROGRAMS += GridTest
GridTest_SOURCES  = ScaFES_GridTest.cpp \
                    ScaFES_GridTest.hpp \
                    ScaFES_mainTest.cpp
GridTest_CPPFLAGS = $(TESTCPPFLAGSxxxCOMMON)
GridTest_LDFLAGS  = $(TESTLDFLAGSxxxCOMMON)
TESTS += GridTest

#------------------------------------------------------------------------------#
check_PROGRAMS += GridTest2D
GridTest2D_SOURCES  = ScaFES_GridTest2D.cpp \
                      ScaFES_GridTest2D.hpp \
                      ScaFES_mainTest.cpp
GridTest2D_CPPFLAGS = $(TESTCPPFLAGSxxxCOMMON)
GridTest2D_LDFLAGS  = $(TESTLDFLAGSxxxCOMMON)
TESTS += GridTest2D

#------------------------------------------------------------------------------#
check_PROGRAMS += GridTest1D
GridTest1D_SOURCES  = ScaFES_GridTest1D.cpp \
                      ScaFES_GridTest1D.hpp \
                    ScaFES_mainTest.cpp
GridTest1D_CPPFLAGS = $(TESTCPPFLAGSxxxCOMMON)
GridTest1D_LDFLAGS  = $(TESTLDFLAGSxxxCOMMON)
TESTS += GridTest1D

#------------------------------------------------------------------------------#
check_PROGRAMS += GridSubTest
GridSubTest_SOURCES  = ScaFES_GridSubTest.cpp \
                       ScaFES_GridSubTest.hpp \
                       ScaFES_mainTest.cpp
GridSubTest_CPPFLAGS = $(TESTCPPFLAGSxxxCOMMON)
GridSubTest_LDFLAGS  = $(TESTLDFLAGSxxxCOMMON)
TESTS += GridSubTest

#------------------------------------------------------------------------------#
check_PROGRAMS += GridGlobalTest
GridGlobalTest_SOURCES  = ScaFES_GridGlobalTest.cpp \
                       ScaFES_GridGlobalTest.hpp \
                       ScaFES_mainTest.cpp
GridGlobalTest_CPPFLAGS = $(TESTCPPFLAGSxxxCOMMON)
GridGlobalTest_LDFLAGS  = $(TESTLDFLAGSxxxCOMMON)
TESTS += GridGlobalTest

#------------------------------------------------------------------------------#
if AMHAVEBOOSTMPI
check_PROGRAMS += BufferTest
BufferTest_SOURCES  = ScaFES_BufferTest.cpp \
                       ScaFES_BufferTest.hpp \
                      ScaFES_mainTest.cpp
BufferTest_CPPFLAGS = $(TESTCPPFLAGSxxxCOMMON)
BufferTest_LDFLAGS  = $(TESTLDFLAGSxxxCOMMON)
TESTS += BufferTest
endif

endif # End of AMHAVEGMOCK.


#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
clean-local:
	-rm -f *.gcno *.opari.inc *.pomp.cpp *.core
endif

