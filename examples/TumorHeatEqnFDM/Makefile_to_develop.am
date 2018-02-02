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
exampledir = $(exec_prefix)/examples/TumorHeatEqnFDM
example_PROGRAMS =
example_SCRIPTS = RUN_HELPER.sh
example_SCRIPTS += RUN_SCALINGTEST_HELPER.sh
example_SCRIPTS += RUN_SCALINGTEST_HYBRID.sh
example_SCRIPTS += RUN_SCALINGTEST_MPI.sh
example_SCRIPTS += RUN_SCALINGTEST_OPENMP.sh
example_SCRIPTS += plotSurface.plt.in
example_SCRIPTS += plotSurface.sh.in
example_SCRIPTS += Parameters.ini
example_SCRIPTS += cases/Bousselham/Figure-5A.ini
example_SCRIPTS += cases/Bousselham/Figure-5B.ini
example_SCRIPTS += cases/Bousselham/Figure-5C.ini
example_SCRIPTS += cases/Bousselham/Figure-7A.ini
example_SCRIPTS += cases/Bousselham/Figure-7B.ini
example_SCRIPTS += cases/Bousselham/Figure-7C.ini
example_SCRIPTS += cases/ScalingTest/5A_120x120x10.ini
example_SCRIPTS += cases/ScalingTest/5A_120x120x50.ini
example_SCRIPTS += cases/ScalingTest/5A_120x120x100.ini
example_SCRIPTS += cases/ScalingTest/5A_240x240x10.ini
example_SCRIPTS += cases/ScalingTest/5A_240x240x50.ini
example_SCRIPTS += cases/ScalingTest/5A_240x240x100.ini
example_SCRIPTS += cases/ScalingTest/5A_240x240x192.ini
example_SCRIPTS += cases/ScalingTest/5A_480x480x192.ini
example_SCRIPTS += cases/ScalingTest/7C_120x120x10.ini
example_SCRIPTS += cases/ScalingTest/7C_120x120x50.ini
example_SCRIPTS += cases/ScalingTest/7C_120x120x100.ini
example_SCRIPTS += cases/ScalingTest/7C_240x240x10.ini
example_SCRIPTS += cases/ScalingTest/7C_240x240x50.ini
example_SCRIPTS += cases/ScalingTest/7C_240x240x100.ini
example_SCRIPTS += cases/ScalingTest/7C_240x240x192.ini
example_SCRIPTS += cases/ScalingTest/7C_480x480x192.ini
bin_PROGRAMS =
include_HEADERS =
noinst_PROGRAMS =


################################################################################
noinst_PROGRAMS += TumorHeatEqnFDM1D
TumorHeatEqnFDM1D_SOURCES = main1D.cpp \
                            TumorHeatEqnFDM.hpp
TumorHeatEqnFDM1D_CPPFLAGS = $(CPPFLAGS_COMMON)
TumorHeatEqnFDM1D_LDFLAGS = $(LDFLAGS_COMMON)
TumorHeatEqnFDM1D_LDADD = $(LDADD_COMMON)
TumorHeatEqnFDM1D_DEPENDENCIES = $(DEPENDENCIES_COMMON)
#------------------------------------------------------------------------------#
noinst_PROGRAMS += TumorHeatEqnFDM2D
TumorHeatEqnFDM2D_SOURCES = main2D.cpp \
                            TumorHeatEqnFDM.hpp
TumorHeatEqnFDM2D_CPPFLAGS = $(CPPFLAGS_COMMON)
TumorHeatEqnFDM2D_LDFLAGS = $(LDFLAGS_COMMON)
TumorHeatEqnFDM2D_LDADD = $(LDADD_COMMON)
TumorHeatEqnFDM2D_DEPENDENCIES = $(DEPENDENCIES_COMMON)
#------------------------------------------------------------------------------#
noinst_PROGRAMS += TumorHeatEqnFDM3D
TumorHeatEqnFDM3D_SOURCES = main3D.cpp \
                            TumorHeatEqnFDM.hpp
TumorHeatEqnFDM3D_CPPFLAGS = $(CPPFLAGS_COMMON)
TumorHeatEqnFDM3D_LDFLAGS = $(LDFLAGS_COMMON)
TumorHeatEqnFDM3D_LDADD = $(LDADD_COMMON)
TumorHeatEqnFDM3D_DEPENDENCIES = $(DEPENDENCIES_COMMON)
#------------------------------------------------------------------------------#
EXTRA_DIST = RUN_HELPER.sh
EXTRA_DIST += RUN_SCALINGTEST_HELPER.sh
EXTRA_DIST += RUN_SCALINGTEST_HYBRID.sh
EXTRA_DIST += RUN_SCALINGTEST_MPI.sh
EXTRA_DIST += RUN_SCALINGTEST_OPENMP.sh
EXTRA_DIST += results/scafes-scalingtests.plt \
        results/ZIH_TAURUS-scaling_strong.plt \
        results/ZIH_TAURUS-scaling_weak.dat \
        results/ZIH_TAURUS-scaling_strong.dat \
        results/ZIH_TAURUS-scaling_weak.plt \
        results/README
EXTRA_DIST += plotSurface.plt.in
EXTRA_DIST += plotSurface.sh.in
EXTRA_DIST += Parameters.ini
EXTRA_DIST += cases/Bousselham/Figure-5A.ini
EXTRA_DIST += cases/Bousselham/Figure-5B.ini
EXTRA_DIST += cases/Bousselham/Figure-5C.ini
EXTRA_DIST += cases/Bousselham/Figure-7A.ini
EXTRA_DIST += cases/Bousselham/Figure-7B.ini
EXTRA_DIST += cases/Bousselham/Figure-7C.ini
EXTRA_DIST += cases/ScalingTest/5A_120x120x10.ini
EXTRA_DIST += cases/ScalingTest/5A_120x120x50.ini
EXTRA_DIST += cases/ScalingTest/5A_120x120x100.ini
EXTRA_DIST += cases/ScalingTest/5A_240x240x10.ini
EXTRA_DIST += cases/ScalingTest/5A_240x240x50.ini
EXTRA_DIST += cases/ScalingTest/5A_240x240x100.ini
EXTRA_DIST += cases/ScalingTest/5A_240x240x192.ini
EXTRA_DIST += cases/ScalingTest/5A_480x480x192.ini
EXTRA_DIST += cases/ScalingTest/7C_120x120x10.ini
EXTRA_DIST += cases/ScalingTest/7C_120x120x50.ini
EXTRA_DIST += cases/ScalingTest/7C_120x120x100.ini
EXTRA_DIST += cases/ScalingTest/7C_240x240x10.ini
EXTRA_DIST += cases/ScalingTest/7C_240x240x50.ini
EXTRA_DIST += cases/ScalingTest/7C_240x240x100.ini
EXTRA_DIST += cases/ScalingTest/7C_240x240x192.ini
EXTRA_DIST += cases/ScalingTest/7C_480x480x192.ini


#------------------------------------------------------------------------------#
if !AMDISTRIBUTESOFTWARE
clean-local:
	-rm -f *.gcno *.opari.inc *.pomp.cpp *.core
endif
