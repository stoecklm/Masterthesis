# SYNOPSIS
#
#   AX_ADOLC([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Test for ADOLC.
#
#   This macro calls:
#
#  AC_SUBST(ADOLC_CPPFLAGS)
#  AC_SUBST(ADOLC_LDFLAGS)
#
#   and sets the preprocessor definition using AC_DEFINE
#   if ADOL-C is available:
#
#  HAVE_ADOLC
#
# Additionally, it defines a variable named
#
#  ac_adolc_available
#
# in order to indicate the presence of ADOL-C.


AC_DEFUN([AX_ADOLC],[

AC_ARG_WITH(
	[adolc],
	AS_HELP_STRING(
		[--with-adolc@<:@=DIR@:>@],
		[home dir of adolc library, default '$ADOLC_HOME']),
	[ac_adolc_dirs="$withval"],
	[ac_adolc_dirs="$ADOLC_HOME"' /usr /usr/local /opt/local ']
)

ac_adolc_available=no

AC_LANG_PUSH([C++])

for ac_adolc_iterate in $ac_adolc_dirs ; do

        CPPFLAGS_SAVED="$CPPFLAGS"
        ADOLC_CPPFLAGS="-I$ac_adolc_iterate/include"
        CPPFLAGS="$CPPFLAGS $ADOLC_CPPFLAGS"
        export CPPFLAGS

        ADOLC_LDFLAGS="-L$ac_adolc_iterate/lib -ladolc"
        LDFLAGS_SAVED="$LDFLAGS"
        LDFLAGS="$LDFLAGS $ADOLC_LDFLAGS"
        export LDFLAGS

        AC_MSG_CHECKING([whether ADOLC is available in $ac_adolc_iterate])

	dnl Test program for adolc: include headers and declare an adouble.
        AC_RUN_IFELSE(
		[AC_LANG_PROGRAM(
			[[@%:@include <adolc/adolc.h>] ],
			[[double x,y; adouble dydx;] ])
		],
		[ac_adolc_available=yes],
		[ac_adolc_available=no]
	)

        LDFLAGS="$LDFLAGS_SAVED"
        export LDFLAGS
        CPPFLAGS="$CPPFLAGS_SAVED"
        export CPPFLAGS

        if test $ac_adolc_available = yes ; then
                AC_MSG_RESULT([yes])
                break
        else
                AC_MSG_RESULT([no])
        fi
done

AC_LANG_POP([C++])

if test $ac_adolc_available = yes ; then
	ifelse([$1], , :, [$1])
	AC_DEFINE(HAVE_ADOLC,[1],[Indicates presence of ADOLC headers and lib])

	export ADOLC_CPPFLAGS
	export ADOLC_LDFLAGS

	AC_SUBST(ADOLC_CPPFLAGS)
	AC_SUBST(ADOLC_LDFLAGS)
else
	ifelse([$2], , :, [$2])
fi

])

