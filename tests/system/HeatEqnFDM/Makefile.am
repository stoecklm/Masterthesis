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
noinst_PROGRAMS += HeatEqnFDM
HeatEqnFDM_SOURCES = main.cpp \
                     HeatEqnFDM.hpp \
                     analyticalSolutions.hpp \
                     boundaryConditions.hpp
HeatEqnFDM_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFDM_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFDM_LDADD = $(LDADD_COMMON)
HeatEqnFDM_DEPENDENCIES = $(DEPENDENCIES_COMMON)
################################################################################
noinst_PROGRAMS += HeatEqnFDM3DConstantDirichlet
HeatEqnFDM3DConstantDirichlet_SOURCES = main3DConstantDirichlet.cpp \
                                        HeatEqnFDM.hpp \
                                        analyticalSolutions.hpp \
                                        boundaryConditions.hpp
HeatEqnFDM3DConstantDirichlet_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFDM3DConstantDirichlet_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFDM3DConstantDirichlet_LDADD = $(LDADD_COMMON)
HeatEqnFDM3DConstantDirichlet_DEPENDENCIES = $(DEPENDENCIES_COMMON)
################################################################################
noinst_PROGRAMS += HeatEqnFDM3DLinearDirichlet
HeatEqnFDM3DLinearDirichlet_SOURCES = main3DLinearDirichlet.cpp \
                                      HeatEqnFDM.hpp \
                                      analyticalSolutions.hpp \
                                      boundaryConditions.hpp
HeatEqnFDM3DLinearDirichlet_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFDM3DLinearDirichlet_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFDM3DLinearDirichlet_LDADD = $(LDADD_COMMON)
HeatEqnFDM3DLinearDirichlet_DEPENDENCIES = $(DEPENDENCIES_COMMON)
################################################################################
noinst_PROGRAMS += HeatEqnFDM3DConstantNeumann
HeatEqnFDM3DConstantNeumann_SOURCES = main3DConstantNeumann.cpp \
                                      HeatEqnFDM.hpp \
                                      analyticalSolutions.hpp \
                                      boundaryConditions.hpp
HeatEqnFDM3DConstantNeumann_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFDM3DConstantNeumann_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFDM3DConstantNeumann_LDADD = $(LDADD_COMMON)
HeatEqnFDM3DConstantNeumann_DEPENDENCIES = $(DEPENDENCIES_COMMON)
################################################################################
noinst_PROGRAMS += HeatEqnFDM3DLinearNeumann
HeatEqnFDM3DLinearNeumann_SOURCES = main3DLinearNeumann.cpp \
                                    HeatEqnFDM.hpp \
                                    analyticalSolutions.hpp \
                                    boundaryConditions.hpp
HeatEqnFDM3DLinearNeumann_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFDM3DLinearNeumann_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFDM3DLinearNeumann_LDADD = $(LDADD_COMMON)
HeatEqnFDM3DLinearNeumann_DEPENDENCIES = $(DEPENDENCIES_COMMON)
################################################################################
noinst_PROGRAMS += HeatEqnFDM3DConstantCauchy
HeatEqnFDM3DConstantCauchy_SOURCES = main3DConstantCauchy.cpp \
                                     HeatEqnFDM.hpp \
                                     analyticalSolutions.hpp \
                                     boundaryConditions.hpp
HeatEqnFDM3DConstantCauchy_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFDM3DConstantCauchy_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFDM3DConstantCauchy_LDADD = $(LDADD_COMMON)
HeatEqnFDM3DConstantCauchy_DEPENDENCIES = $(DEPENDENCIES_COMMON)
################################################################################
noinst_PROGRAMS += HeatEqnFDM3DLinearCauchy
HeatEqnFDM3DLinearCauchy_SOURCES = main3DLinearCauchy.cpp \
                                   HeatEqnFDM.hpp \
                                   analyticalSolutions.hpp \
                                   boundaryConditions.hpp
HeatEqnFDM3DLinearCauchy_CPPFLAGS = $(CPPFLAGS_COMMON)
HeatEqnFDM3DLinearCauchy_LDFLAGS = $(LDFLAGS_COMMON)
HeatEqnFDM3DLinearCauchy_LDADD = $(LDADD_COMMON)
HeatEqnFDM3DLinearCauchy_DEPENDENCIES = $(DEPENDENCIES_COMMON)
################################################################################
noinst_PROGRAMS += testAnalyticalSolutions
testAnalyticalSolutions_SOURCES = testAnalyticalSolutions.cpp \
                                  analyticalSolutions.hpp
testAnalyticalSolutions_CPPFLAGS = $(CPPFLAGS_COMMON)
testAnalyticalSolutions_LDFLAGS = $(LDFLAGS_COMMON)
testAnalyticalSolutions_LDADD = $(LDADD_COMMON)
testAnalyticalSolutions_DEPENDENCIES = $(DEPENDENCIES_COMMON)

#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
clean-local:
	-rm -f *.gcno *.opari.inc *.pomp.cpp *.core
endif
