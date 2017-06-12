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
include_HEADERS=
bin_PROGRAMS=
noinst_HEADERS=
check_HEADERS=
check_PROGRAMS=
EXTRA_DIST=

################################################################################
noinst_HEADERS +=
check_PROGRAMS += HeatEqnFdmStencilWidthOne1D
HeatEqnFdmStencilWidthOne1D_SOURCES = \
    HeatEqnFdmStencilWidthOne.hpp \
    main1D.cpp
HeatEqnFdmStencilWidthOne1D_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFdmStencilWidthOne1D_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFdmStencilWidthOne1D_LDADD = $(LDADD_COMMON)
HeatEqnFdmStencilWidthOne1D_DEPENDENCIES = $(DEPENDENCIES_COMMON)
#------------------------------------------------------------------------------#
check_PROGRAMS += HeatEqnFdmStencilWidthOne2D
HeatEqnFdmStencilWidthOne2D_SOURCES = \
    HeatEqnFdmStencilWidthOne.hpp \
    main2D.cpp
HeatEqnFdmStencilWidthOne2D_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFdmStencilWidthOne2D_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFdmStencilWidthOne2D_LDADD = $(LDADD_COMMON)
HeatEqnFdmStencilWidthOne2D_DEPENDENCIES = $(DEPENDENCIES_COMMON)
#------------------------------------------------------------------------------#
check_PROGRAMS += HeatEqnFdmStencilWidthOne3D
HeatEqnFdmStencilWidthOne3D_SOURCES = \
    HeatEqnFdmStencilWidthOne.hpp \
    main3D.cpp
HeatEqnFdmStencilWidthOne3D_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFdmStencilWidthOne3D_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFdmStencilWidthOne3D_LDADD = $(LDADD_COMMON)
HeatEqnFdmStencilWidthOne3D_DEPENDENCIES = $(DEPENDENCIES_COMMON)

#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
clean-local:
	-rm -f *.gcno *.opari.inc *.pomp.cpp *.core
endif

