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
SCAFES_TESTS_DIR=$(top_srcdir)/examples

CPPFLAGS_COMMON= \
        -I$(top_srcdir)/src \
        -I$(top_builddir)/src \
        @NETCDF_CPPFLAGS@ \
        @BOOST_CPPFLAGS@ \
        @MPI_CPPFLAGS@ \
        @OPENMP_CXXFLAGS@ \
        @ADOLC_CPPFLAGS@ \
        $(VTxxxCPPFLAGS) \
        $(VTOPARIxxxCPPFLAGS)

LDFLAGS_COMMON = $(VTxxxCPPFLAGS) \
                 $(VTOPARIxxxCPPFLAGS)

LDADD_COMMON = \
        $(top_builddir)/src/libscafes.la \
        @BOOST_LDFLAGS@ \
        @BOOST_PROGRAM_OPTIONS_LIB@ \
        @BOOST_REGEX_LIB@ \
        @BOOST_SERIALIZATION_LIB@ \
        @MPI_LDFLAGS@ \
        @BOOST_MPI_LIB@ \
        @OPENMP_CXXFLAGS@ \
        @NETCDF_LDFLAGS@ \
        @ADOLC_LDFLAGS@

DEPENDENCIES_COMMON= $(top_builddir)/src/libscafes.la


################################################################################
bin_PROGRAMS =
include_HEADERS =
noinst_PROGRAMS =


################################################################################
noinst_PROGRAMS += HeatEqnFDMInitFile1D
HeatEqnFDMInitFile1D_SOURCES = main1D.cpp \
                     HeatEqnFDMInitFile.hpp \
                     analyticalFunctions.hpp
HeatEqnFDMInitFile1D_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFDMInitFile1D_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFDMInitFile1D_LDADD = $(LDADD_COMMON)
HeatEqnFDMInitFile1D_DEPENDENCIES = $(DEPENDENCIES_COMMON)
EXTRA_DIST = RUN_1D.sh
################################################################################
noinst_PROGRAMS += HeatEqnFDMInitFile2D
HeatEqnFDMInitFile2D_SOURCES = main2D.cpp \
                     HeatEqnFDMInitFile.hpp \
                     analyticalFunctions.hpp
HeatEqnFDMInitFile2D_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFDMInitFile2D_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFDMInitFile2D_LDADD = $(LDADD_COMMON)
HeatEqnFDMInitFile2D_DEPENDENCIES = $(DEPENDENCIES_COMMON)
EXTRA_DIST += RUN_2D.sh
################################################################################
noinst_PROGRAMS += HeatEqnFDMInitFile3D
HeatEqnFDMInitFile3D_SOURCES = main3D.cpp \
                     HeatEqnFDMInitFile.hpp \
                     analyticalFunctions.hpp
HeatEqnFDMInitFile3D_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFDMInitFile3D_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFDMInitFile3D_LDADD = $(LDADD_COMMON)
HeatEqnFDMInitFile3D_DEPENDENCIES = $(DEPENDENCIES_COMMON)
EXTRA_DIST += RUN_3D.sh
EXTRA_DIST += RUN.sh
EXTRA_DIST += parameters.ini

#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
clean-local:
	-rm -f *.gcno *.opari.inc *.pomp.cpp *.core
endif
