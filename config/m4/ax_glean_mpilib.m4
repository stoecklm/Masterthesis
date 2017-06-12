# ===========================================================================
#          http://www.gnu.org/software/autoconf-archive/ax_mpi.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_GLEAN_MPILIB([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile programs that use MPI
#   (Message Passing Interface), a standard API for parallel process
#   communication (see http://www-unix.mcs.anl.gov/mpi/)
#
#   Dependent on the current language, it will unwind the the MPI frontend
#   compiler's --show-me output, and feed this into MPI_LDFLAGS and
#   MPI_CPPFLAGS, and subsequently call AC_SUBST on both.
#   Going this way allows to pick out all MPI compiler flags (since it
#   is just another library after all), without needing to override the
#   compiler: just use CC instead of MPICC etc.
#
# LICENSE
#
#   Copyright (c) 2012 Sebastian Hegler <sebastian.hegler@tu-dresden.de>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_GLEAN_MPILIB], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

AC_LANG_CASE(
[C], [
	AC_REQUIRE([AC_PROG_CC])
	AC_ARG_VAR(MPICC,[MPI C compiler command])
	AC_CHECK_PROGS(MPICC, mpigcc mpicc hcc mpxlc_r mpxlc mpcc cmpicc, $CC)
	ax_mpi_save_compiler="$MPICC"
],
[C++], [
	AC_REQUIRE([AC_PROG_CXX])
	AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
        dnl echo "BEFORE MPICXX=$MPICXX"
        dnl echo "BEFORE CXX=$CXX"
        dnl Order is important for Juqueen. [KF, 04/22/2015]
        dnl If gcc should be used for cross-compiling on Juqueen then mpig++ must be 
        dnl set to first position.
        dnl If xlc should be used then mpixlcxx should be used.
        AC_CHECK_PROGS(MPICXX, mpig++ mpic++ mpicxx mpiCC hcp mpxlC_r mpxlC mpCC cmpic++, $CXX)
	ax_mpi_save_compiler="$MPICXX"
        dnl echo "AFTER MPICXX=$MPICXX"
],
[Fortran 77], [
	AC_REQUIRE([AC_PROG_F77])
	AC_ARG_VAR(MPIF77,[MPI Fortran 77 compiler command])
	AC_CHECK_PROGS(MPIF77, mpif77 hf77 mpxlf_r mpxlf mpf77 cmpifc, $F77)
	ax_mpi_save_compiler="$MPIF77"
],
[Fortran], [
	AC_REQUIRE([AC_PROG_FC])
	AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
	AC_CHECK_PROGS(MPIFC, mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c, $FC)
	ax_mpi_save_compiler="$MPIFC"
])


AC_MSG_CHECKING([whether MPI wrapper compiler supports --show-me])
dnl call --show-me
ax_mpi_showme=`$ax_mpi_save_compiler -show 2>&1`

MPI_LDFLAGS=''
MPI_CPPFLAGS=''


for ax_mpi_show in $ax_mpi_showme ; do
	case ${ax_mpi_show:0:3} in
		-Wl)
 			MPI_LDFLAGS=$MPI_LDFLAGS' '$ax_mpi_show
			;;
		*)
			;;
	esac
	case ${ax_mpi_show:0:2} in
		-L)
 			MPI_LDFLAGS=$MPI_LDFLAGS' '$ax_mpi_show
			;;
		-l)
			MPI_LDFLAGS=$MPI_LDFLAGS' '$ax_mpi_show
			;;
		-I)
			MPI_CPPFLAGS=$MPI_CPPFLAGS' '$ax_mpi_show
			;;
		*)
			;;
	esac
done

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPI_LDFLAGS"; then
        $2
	AC_MSG_RESULT([no])
else
    ifelse([$1], , :, [$1])
    AC_DEFINE(HAVE_MPI,[1],[Indicates presence of MPI headers and lib])

    AC_MSG_RESULT([yes])

    AC_SUBST(MPI_LDFLAGS)
    AC_SUBST(MPI_CPPFLAGS)
    echo "MPI_LDLFLAGS=$MPI_LDFLAGS"
    echo "MPI_CPPFLAGS=$MPI_CPPFLAGS"
fi
])dnl AX_GLEAN_MPILIB

