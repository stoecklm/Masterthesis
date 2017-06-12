#AUTOMAKE_OPTIONS= subdir-objects

################################################################################
if !AMDISTRIBUTESOFTWARE
if AMHAVEVT
VT_OPARI_CPPFLAGS = @VT_CPPFLAGS@ \
-vt:opari "-rcfile $(top_builddir)/src/scafes.opari.rc -table $(top_builddir)/src/scafes.opari.tab.c"
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
        @OPENMP_CXXFLAGS@
if !AMDISTRIBUTESOFTWARE
CPPFLAGS_COMMON += $(VT_OPARI_CPPFLAGS)
endif

if !AMDISTRIBUTESOFTWARE
LDFLAGS_COMMON = $(VT_OPARI_CPPFLAGS)
endif
LDADD_COMMON = \
        $(top_builddir)/src/libscafes.la \
        @BOOST_LDFLAGS@ \
        @BOOST_PROGRAM_OPTIONS_LIB@ \
        @BOOST_REGEX_LIB@ \
        @BOOST_SERIALIZATION_LIB@ \
        @MPI_LDFLAGS@ \
        @BOOST_MPI_LIB@ \
        @OPENMP_CXXFLAGS@ \
        @NETCDF_LDFLAGS@
DEPENDENCIES_COMMON= $(top_builddir)/src/libscafes.la

################################################################################
include_HEADERS=
noinst_PROGRAMS=
EXTRA_DIST=

################################################################################
include_HEADERS +=
noinst_PROGRAMS += MaxwellEqnsFDTD
MaxwellEqnsFDTD_SOURCES = main.cpp
MaxwellEqnsFDTD_CPPFLAGS = $(CPPFLAGS_COMMON)
MaxwellEqnsFDTD_LDFLAGS = $(LDFLAGS_COMMON)
MaxwellEqnsFDTD_LDADD = $(LDADD_COMMON)
MaxwellEqnsFDTD_DEPENDENCIES = $(DEPENDENCIES_COMMON)

EXTRA_DIST += $(SCAFES_TESTS_DIR)/MaxwellEqnsFDTD/RUN.sh


#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
clean-local:
	-rm -f *.gcno *.opari.inc *.pomp.cpp *.core
endif
