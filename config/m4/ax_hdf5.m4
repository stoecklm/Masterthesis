dnl ######################################################################
dnl
dnl Purpose: Determine the locations of C++ Interface for
dnl          the hdf5 includes and libraries
dnl
dnl ######################################################################

dnl ######################################################################
dnl
dnl Allow the user to specify an overall hdf5 directory.  If specified,
dnl we look for include and lib under this.
dnl
dnl ######################################################################
dnl SYNOPSIS
dnl
dnl AX_HDF5([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])

AC_DEFUN([AX_HDF5],[

dnl guess from env, or use given value
AC_ARG_WITH([hdf5],
	AS_HELP_STRING([--with-hdf5@<:@=DIR@:>@],[location of hdf5 installation, default $H5_HOME]),
	[HDF5_DIR="$withval"],[HDF5_DIR="$H5_HOME"])

dnl override include dir if explicitly given
AC_ARG_WITH([hdf5-incdir],
	AS_HELP_STRING([--with-hdf5-incdir@<:@=DIR@:>@],[location of hdf5 includes, overrides hdf5-dir]),
	[HDF5_INCDIR="$withval"],[HDF5_INCDIR="$HDF5_DIR"/include])
dnl override lib dir if explicitly given
AC_ARG_WITH([hdf5-libdir],
	AS_HELP_STRING([--with-hdf5-libdir@<:@=DIR@:>@],[location of hdf5 library, overrides hdf5-dir]),
	[HDF5_LIBDIR="$withval"],[HDF5_LIBDIR="$HDF5_DIR"/lib])


CPPFLAGS_SAVED="$CPPFLAGS"
HDF5_CPPFLAGS="-I$HDF5_INCDIR"
CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
export CPPFLAGS

dnl ################ C interface
#AC_CHECK_LIB([hdf5],[H5open],
#        [AC_CHECK_HEADER(
#		[hdf5.h],
#		[ifelse([$1], , :, [$1])],
#		[ifelse([$2], , :, [$2])])
#	],
#       [ifelse([$2], , :, [$2])],
#        [-L$HDF5_LIBDIR]
#)
#AC_DEFINE(HAVE_HDF5,[1],[Indicates presence of basic HDF5 library])





AC_CHECK_LIB([hdf5],[H5open],
        [AC_CHECK_HEADER([hdf5.h],[ac_hdf5_available=yes],[ac_hdf5_available=no],
                [#include <hdf5.h>])],
        [ac_hdf5_available=no],
        [-L$HDF5_LIBDIR]
)
if test $ac_hdf5_available = yes ; then

        ifelse([$1], , , [$1])
        AC_DEFINE(HAVE_HDF5,[1],[Indicates presence of basic HDF5 library])

        export HDF5_CPPFLAGS
        export HDF5_LDFLAGS

        AC_SUBST(HDF5_CPPFLAGS)
        AC_SUBST(HDF5_LDFLAGS)
else
        ifelse([$2], , , [$2])
fi





dnl ################ C LT interface
HDF5_LDFLAGS="-L$HDF5_LIBDIR -lhdf5_hl -lhdf5"
LDFLAGS="$LDFLAGS_SAVED $HDF5_LDFLAGS"
export LDFLAGS

AC_CHECK_LIB([hdf5_hl],[H5LTread_dataset],
        [AC_CHECK_HEADER([hdf5_hl.h],[h5l_ok=yes],[h5l_ok=no],
			[#include <hdf5.h>])],
        [h5l_ok=no],
        [$HDF5_LDFLAGS]
)

if test $h5l_ok = yes ; then
	AC_DEFINE(HAVE_HDF5_LT,[1],[Indicates presence of HDF5 LT library])
	HDF5_LT="-lhdf5_hl"
else
	AC_MSG_WARN(*** HDF5 Lite library not found ***)
	HDF5_LT=""
fi

LDFLAGS="$LDFLAGS_SAVED"
HDF5_LDFLAGS="-L$HDF5_LIBDIR $HDF5_LT -lhdf5"

export HDF5_LDFLAGS
export HDF5_CPPFLAGS

AC_SUBST(HDF5_CPPFLAGS)
AC_SUBST(HDF5_LDFLAGS)

LDFLAGS="$LDFLAGS_SAVED"
export LDFLAGS
CPPFLAGS="$CPPFLAGS_SAVED"
export CPPFLAGS

])

