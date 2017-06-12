# SYNOPSIS
#
#   AX_MPI([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Test for MPI.
#
#   This macro calls:
#
#  AC_SUBST(MPI_CPPFLAGS)
#  AC_SUBST(MPI_LDFLAGS)
#
#   and sets the preprocessor definition using AC_DEFINE
#   if MPI is available:
#
#  HAVE_MPI
#
# Additionally, it defines a variable named
#
#  ac_mpi_available
#
# in order to indicate the presence of MPI.


AC_DEFUN([AX_MPI],[

AC_ARG_WITH(
    [mpi],
    AS_HELP_STRING(
        [--with-mpi@<:@=DIR@:>@],
        [home dir of mpi library, default '$MPI_HOME']),
    [ac_mpi_dirs="$withval"],
    [ac_mpi_dirs="$MPI_HOME"' /usr /usr/local /opt/local ']
)

ac_mpi_available=no

AC_LANG_PUSH([C++])

for ac_mpi_iterate in $ac_mpi_dirs ; do

        CPPFLAGS_SAVED="$CPPFLAGS"
        MPI_CPPFLAGS="-I$ac_mpi_iterate/include"
        CPPFLAGS="$CPPFLAGS $MPI_CPPFLAGS"
        export CPPFLAGS

        MPI_LDFLAGS="-L$ac_mpi_iterate/lib -lmpi++ -lmpi"
        LDFLAGS_SAVED="$LDFLAGS"
        LDFLAGS="$LDFLAGS $MPI_LDFLAGS"
        export LDFLAGS

        AC_MSG_CHECKING([whether MPI is available in $ac_mpi_iterate])

AC_CHECK_LIB([mpi],[MPI_Init],
    [AC_CHECK_HEADER([mpi.h],[ac_mpi_available=yes],[ac_mpi_available=no],
        [#include <mpi.h>])],
    [ac_mpi_available=no],
    [$MPI_LDFLAGS]
)


        LDFLAGS="$LDFLAGS_SAVED"
        export LDFLAGS
        CPPFLAGS="$CPPFLAGS_SAVED"
        export CPPFLAGS

        if test $ac_mpi_available = yes ; then
                AC_MSG_RESULT([yes])
                break
        else
                AC_MSG_RESULT([no])
        fi
done

AC_LANG_POP([C++])

if test $ac_mpi_available = yes ; then
    ifelse([$1], , :, [$1])
    AC_DEFINE(HAVE_MPI,[1],[Indicates presence of MPI headers and lib])

    export MPI_CPPFLAGS
    export MPI_LDFLAGS

    AC_SUBST(MPI_CPPFLAGS)
    AC_SUBST(MPI_LDFLAGS)
else
    ifelse([$2], , :, [$2])
fi

])

