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
check_PROGRAMS=
EXTRA_DIST=

################################################################################
noinst_HEADERS +=
check_PROGRAMS += HeatEqnFdmStencilWidthVar1D
HeatEqnFdmStencilWidthVar1D_SOURCES = \
    HeatEqnFdmStencilWidthVar.hpp \
    main1D.cpp
HeatEqnFdmStencilWidthVar1D_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFdmStencilWidthVar1D_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFdmStencilWidthVar1D_LDADD = $(LDADD_COMMON)
HeatEqnFdmStencilWidthVar1D_DEPENDENCIES = $(DEPENDENCIES_COMMON)
#------------------------------------------------------------------------------#
noinst_HEADERS +=
check_PROGRAMS += HeatEqnFdmStencilWidthVar2D
HeatEqnFdmStencilWidthVar2D_SOURCES = \
    HeatEqnFdmStencilWidthVar.hpp \
    main2D.cpp
HeatEqnFdmStencilWidthVar2D_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFdmStencilWidthVar2D_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFdmStencilWidthVar2D_LDADD = $(LDADD_COMMON)
HeatEqnFdmStencilWidthVar2D_DEPENDENCIES = $(DEPENDENCIES_COMMON)
#------------------------------------------------------------------------------#
noinst_HEADERS +=
check_PROGRAMS += HeatEqnFdmStencilWidthVar3D
HeatEqnFdmStencilWidthVar3D_SOURCES = \
    HeatEqnFdmStencilWidthVar.hpp \
    main3D.cpp
HeatEqnFdmStencilWidthVar3D_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFdmStencilWidthVar3D_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFdmStencilWidthVar3D_LDADD = $(LDADD_COMMON)
HeatEqnFdmStencilWidthVar3D_DEPENDENCIES = $(DEPENDENCIES_COMMON)
#------------------------------------------------------------------------------#
noinst_HEADERS +=
check_PROGRAMS += HeatEqnFdmStencilWidthVar4D
HeatEqnFdmStencilWidthVar4D_SOURCES = \
    HeatEqnFdmStencilWidthVar.hpp \
    main4D.cpp
HeatEqnFdmStencilWidthVar4D_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFdmStencilWidthVar4D_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFdmStencilWidthVar4D_LDADD = $(LDADD_COMMON)
HeatEqnFdmStencilWidthVar4D_DEPENDENCIES = $(DEPENDENCIES_COMMON)

#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
clean-local:
	-rm -f *.gcno *.opari.inc *.pomp.cpp *.core
endif

