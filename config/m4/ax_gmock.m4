# SYNOPSIS
#
#   AX_GMOCK([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Test for Gmock.
#
#   This macro calls:
#
#  AC_SUBST(GMOCK_CPPFLAGS)
#  AC_SUBST(GMOCK_LDFLAGS)
#
#   and sets the preprocessor definition using AC_DEFINE
#   if VampirTrace is available:
#
#  HAVE_GMOCK
#
# Additionally, it defines a variable named
#
#  ac_gmock_available
# 
# in order to indicate the presence of GMOCK.


AC_DEFUN([AX_GMOCK],[

AC_ARG_WITH(
	[gmock],
	AS_HELP_STRING(
		[--with-gmock@<:@=DIR@:>@],
		[home dir of Gmock library, default '$GMOCK_HOME']),
	[ac_gmock_dirs="$withval"],
	[ac_gmock_dirs="$GMOCK_HOME"' /usr /usr/local /opt/local ']
)

ac_gmock_available=no

AC_LANG_PUSH([C++])

for ac_gmock_iterate in $ac_gmock_dirs ; do

        CPPFLAGS_SAVED="$CPPFLAGS"
        GMOCK_CPPFLAGS="-I$ac_gmock_iterate/include"
        CPPFLAGS="$CPPFLAGS $GMOCK_CPPFLAGS"
        export CPPFLAGS

        GMOCK_LDFLAGS="-L$ac_gmock_iterate/lib -pthread -lgtest -lgmock"
        LDFLAGS_SAVED="$LDFLAGS"
        LDFLAGS="$LDFLAGS $GMOCK_LDFLAGS"
        export LDFLAGS

        AC_MSG_CHECKING([whether GMOCK is available in $ac_gmock_iterate])

        LD_LIBRARY_PATH_SAVED=$LD_LIBRARY_PATH
        LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ac_gmock_iterate/lib
        export LD_LIBRARY_PATH

	dnl Test program for gmock: include headers.
        AC_RUN_IFELSE(
		[AC_LANG_PROGRAM(
			[[@%:@include <gmock/gmock.h>] ],
			[ ])
		],
		[ac_gmock_available=yes],
		[ac_gmock_available=no]
	)

        LDFLAGS="$LDFLAGS_SAVED"
        export LDFLAGS
        CPPFLAGS="$CPPFLAGS_SAVED"
        export CPPFLAGS
        LD_LIBRARY_PATH="$LD_LIBRARY_PATH_SAVED"
        export LD_LIBRARY_PATH

        if test $ac_gmock_available = yes ; then
                AC_MSG_RESULT([yes])
                break
        else
                AC_MSG_RESULT([no])
        fi
done

AC_LANG_POP([C++])

if test $ac_gmock_available = yes ; then
	ifelse([$1], , :, [$1])
	AC_DEFINE(HAVE_GMOCK,[1],[Indicates presence of GMOCK headers and lib])

	export GMOCK_CPPFLAGS
	export GMOCK_LDFLAGS

	AC_SUBST(GMOCK_CPPFLAGS)
	AC_SUBST(GMOCK_LDFLAGS)
else
	ifelse([$2], , :, [$2])
fi

])

