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
exampledir = $(exec_prefix)/examples/MRIDataPBHEqnFDM
example_PROGRAMS =
example_SCRIPTS = RUN_HELPER.sh
example_SCRIPTS += Parameters.ini
bin_PROGRAMS =
include_HEADERS =
noinst_PROGRAMS =


################################################################################
noinst_PROGRAMS += MRIDataPBHEqnFDM3D
MRIDataPBHEqnFDM3D_SOURCES = main3D.cpp \
                             MRIDataPBHEqnFDM.hpp
MRIDataPBHEqnFDM3D_CPPFLAGS = $(CPPFLAGS_COMMON)
MRIDataPBHEqnFDM3D_LDFLAGS = $(LDFLAGS_COMMON)
MRIDataPBHEqnFDM3D_LDADD = $(LDADD_COMMON)
MRIDataPBHEqnFDM3D_DEPENDENCIES = $(DEPENDENCIES_COMMON)
#------------------------------------------------------------------------------#
EXTRA_DIST = RUN_HELPER.sh
EXTRA_DIST += Parameters.ini

#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
clean-local:
	-rm -f *.gcno *.opari.inc *.pomp.cpp *.core
endif
