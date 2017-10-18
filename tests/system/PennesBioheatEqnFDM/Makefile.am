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
noinst_PROGRAMS=
check_HEADERS=
check_PROGRAMS=
EXTRA_DIST=

################################################################################
noinst_PROGRAMS += PennesBioheatEqnFDM
PennesBioheatEqnFDM_SOURCES = main.cpp \
                              PennesBioheatEqnFDM.hpp \
                              analyticalSolutions.hpp \
                              boundaryConditions.hpp
PennesBioheatEqnFDM_CPPFLAGS = $(CPPFLAGS_COMMON)
PennesBioheatEqnFDM_LDFLAGS = $(LDFLAGS_COMMON)
PennesBioheatEqnFDM_LDADD = $(LDADD_COMMON)
PennesBioheatEqnFDM_DEPENDENCIES = $(DEPENDENCIES_COMMON)
################################################################################
noinst_PROGRAMS += PennesBioheatEqnFDM3DConstantDirichlet
PennesBioheatEqnFDM3DConstantDirichlet_SOURCES = main3DConstantDirichlet.cpp \
                                                 PennesBioheatEqnFDM.hpp \
                                                 analyticalSolutions.hpp \
                                                 boundaryConditions.hpp
PennesBioheatEqnFDM3DConstantDirichlet_CPPFLAGS = $(CPPFLAGS_COMMON)
PennesBioheatEqnFDM3DConstantDirichlet_LDFLAGS = $(LDFLAGS_COMMON)
PennesBioheatEqnFDM3DConstantDirichlet_LDADD = $(LDADD_COMMON)
PennesBioheatEqnFDM3DConstantDirichlet_DEPENDENCIES = $(DEPENDENCIES_COMMON)
################################################################################
noinst_PROGRAMS += PennesBioheatEqnFDM3DLinearDirichlet
PennesBioheatEqnFDM3DLinearDirichlet_SOURCES = main3DLinearDirichlet.cpp \
                                               PennesBioheatEqnFDM.hpp \
                                               analyticalSolutions.hpp \
                                               boundaryConditions.hpp
PennesBioheatEqnFDM3DLinearDirichlet_CPPFLAGS = $(CPPFLAGS_COMMON)
PennesBioheatEqnFDM3DLinearDirichlet_LDFLAGS = $(LDFLAGS_COMMON)
PennesBioheatEqnFDM3DLinearDirichlet_LDADD = $(LDADD_COMMON)
PennesBioheatEqnFDM3DLinearDirichlet_DEPENDENCIES = $(DEPENDENCIES_COMMON)

#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
clean-local:
	-rm -f *.gcno *.opari.inc *.pomp.cpp *.core
endif
