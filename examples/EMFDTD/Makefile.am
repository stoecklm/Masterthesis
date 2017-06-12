#------------------------------------------------------------------------------#
CPPFLAGS_COMMON= \
        -I$(top_srcdir)/src \
        -I$(top_builddir)/src \
        @BOOST_CPPFLAGS@ \
        @MPI_CPPFLAGS@ \
        @OPENMP_CXXFLAGS@ \
        @HDF5_CPPFLAGS@ \
        @NETCDF_CPPFLAGS@ \
        @ADOLC_CPPFLAGS@

LDADD_COMMON = \
        -L$(top_builddir)/src -lscafes \
        @BOOST_LDFLAGS@ \
        @BOOST_PROGRAM_OPTIONS_LIB@ \
        @BOOST_REGEX_LIB@ \
        @BOOST_SERIALIZATION_LIB@ \
        @BOOST_MPI_LIB@ \
        @MPI_LDFLAGS@ \
        @OPENMP_CXXFLAGS@ \
        @HDF5_LDFLAGS@ \
        @NETCDF_LDFLAGS@ \
        @ADOLC_LDFLAGS@

DEPENDENCIES_COMMON= $(top_builddir)/src/libscafes.la

################################################################################
bin_PROGRAMS=EMFDTDsimulate
include_HEADERS=
noinst_PROGRAMS=
EXTRA_DIST=

include_HEADERS += \
                   EMFDTD.hpp \
                   physical_constants.hpp \
                   CPMLParams.hpp \
                   conversion.hpp \
                   sources/ModulatedGaussianSource.hpp \
                   geometry/IndicatorEllipsoid.hpp \
                   stencils/CompleteCPMLCurlStencil.hpp

################################################################################
noinst_PROGRAMS += EMFDTDsimulate
EMFDTDsimulate_SOURCES = simulate.cpp
EMFDTDsimulate_CPPFLAGS = -DEMFDTD_WRITE_DATA $(CPPFLAGS_COMMON)
EMFDTDsimulate_LDFLAGS = $(LDFLAGS_COMMON)
EMFDTDsimulate_LDADD = $(LDADD_COMMON)
EMFDTDsimulate_DEPENDENCIES = $(DEPENDENCIES_COMMON)
