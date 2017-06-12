# SYNOPSIS
#
#   AX_NETCDF([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Test for NETCDF.
#
#   This macro calls:
#
#  AC_SUBST(NETCDF_CPPFLAGS)
#  AC_SUBST(NETCDF_LDFLAGS)
#
#   And sets:
#
#  HAVE_NETCDF
#
# Requires HDF5 as a prerequisite.


AC_DEFUN([AX_NETCDF],[

AC_ARG_WITH(
    [netcdf],
    AS_HELP_STRING(
        [--with-netcdf@<:@=DIR@:>@],
        [home dir of netcdf library, default '$NETCDF_HOME']),
    [ac_netcdf_dir="$withval"],
    [ac_netcdf_dir="$NETCDF_HOME"]
)

AC_LANG_PUSH([C])

CPPFLAGS_SAVED="$CPPFLAGS"
NETCDF_CPPFLAGS="-I$ac_netcdf_dir/include"
CPPFLAGS="$CPPFLAGS $NETCDF_CPPFLAGS $HDF5_CPPFLAGS"
export CPPFLAGS

NETCDF_LDFLAGS="-L$ac_netcdf_dir/lib -lnetcdf"
LDFLAGS_SAVED="$LDFLAGS"
LDFLAGS="$LDFLAGS $NETCDF_LDFLAGS $HDF5_LDFLAGS"
export LDFLAGS

AC_CHECK_LIB([netcdf],[nc_def_var],
    [AC_CHECK_HEADER([netcdf.h],[ac_netcdf_available=yes],
                     [ac_netcdf_available=no],
        [#include <netcdf.h>])],
    [ac_netcdf_available=no],
    [$NETCDF_LDFLAGS]
)

AC_CHECK_LIB([netcdf],[nc_def_var],
    [AC_CHECK_HEADER([netcdf_par.h],[ac_netcdfpar_available=yes],
                     [ac_netcdfpar_available=no], [#include <netcdf_par.h>])],
    [ac_netcdfpar_available=no],
    [$NETCDF_LDFLAGS]
)


LDFLAGS="$LDFLAGS_SAVED"
export LDFLAGS
CPPFLAGS="$CPPFLAGS_SAVED"
export CPPFLAGS

AC_LANG_POP([C])

if test $ac_netcdf_available = yes ; then

    ifelse([$1], , , [$1])
    AC_DEFINE(HAVE_NETCDF,[1],[Indicates presence of NETCDF headers and lib])
    if test $ac_netcdfpar_available = yes ; then
        AC_DEFINE(HAVE_NETCDF_PAR,[1],
                 [Indicates presence of parallel NETCDF headers and lib])
    fi

    export NETCDF_CPPFLAGS
    export NETCDF_LDFLAGS

    AC_SUBST(NETCDF_CPPFLAGS)
    AC_SUBST(NETCDF_LDFLAGS)
else
    ifelse([$2], , , [$2])
fi

])

