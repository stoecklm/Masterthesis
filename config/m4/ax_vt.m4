# SYNOPSIS
#
#   AX_VT([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Test for VampirTrace.
#
#   This macro calls:
#
#  AC_SUBST(VT_CPPFLAGS)
#  AC_SUBST(VT_LDFLAGS)
#
#   and sets the preprocessor definition using AC_DEFINE
#   if VampirTrace is available:
#
#  HAVE_VT
#
# Additionally, it defines a variable named
#
#  ac_vt_available
# 
# in order to indicate the presence of VT.


AC_DEFUN([AX_VT],[

AC_ARG_WITH(
	[vt],
	AS_HELP_STRING(
		[--with-vt@<:@=DIR@:>@],
		[home dir of VampirTrace library, default '$VT_HOME']),
	[ac_vt_dirs="$withval"],
	[ac_vt_dirs="$VT_HOME"' /usr /usr/local /opt/local ']
)

ac_vt_available=no

AC_LANG_PUSH([C++])

for ac_vt_iterate in $ac_vt_dirs ; do

        CPPFLAGS_SAVED="$CPPFLAGS"
        VT_CPPFLAGS="-I$ac_vt_iterate/include"
        CPPFLAGS="$CPPFLAGS $VT_CPPFLAGS"
        export CPPFLAGS

        VT_LDFLAGS="-L$ac_vt_iterate/lib -lvt"
        LDFLAGS_SAVED="$LDFLAGS"
        LDFLAGS="$LDFLAGS $VT_LDFLAGS"
        export LDFLAGS

        AC_MSG_CHECKING([whether VT is available in $ac_vt_iterate])

	dnl Test program for vt: include headers.
        AC_RUN_IFELSE(
		[AC_LANG_PROGRAM(
			[[@%:@include <vampirtrace/vt_user.h>] ],
			[ ])
		],
		[ac_vt_available=yes],
		[ac_vt_available=no]
	)

        LDFLAGS="$LDFLAGS_SAVED"
        export LDFLAGS
        CPPFLAGS="$CPPFLAGS_SAVED"
        export CPPFLAGS

        if test $ac_vt_available = yes ; then
                AC_MSG_RESULT([yes])
                break
        else
                AC_MSG_RESULT([no])
        fi
done

AC_LANG_POP([C++])

if test $ac_vt_available = yes ; then
	ifelse([$1], , :, [$1])
	AC_DEFINE(HAVE_VT,[1],[Indicates presence of VT headers and lib])

	export VT_CPPFLAGS
	export VT_LDFLAGS

	AC_SUBST(VT_CPPFLAGS)
	AC_SUBST(VT_LDFLAGS)
else
	ifelse([$2], , :, [$2])
fi

])

