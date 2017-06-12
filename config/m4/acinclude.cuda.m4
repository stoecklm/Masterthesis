# SYNOPSIS
#
#   AX_CUDA([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Test for CUDA.
#
#   This macro calls:
#
#  AC_SUBST(CUDA_CFLAGS)
#  AC_SUBST(CUDA_LDFLAGS)
#  AC_SUBST(CUDA_CC)
#
#   And sets:
#
#  HAVE_CUDA
#

AC_DEFUN([AX_CUDA],[

AC_ARG_WITH(
        [cuda],
        AS_HELP_STRING(
                [--with-cuda@<:@=DIR@:>@],
                [home dir of cuda library, default '$CUDA_HOME']),
        [ac_cuda_dir="$withval"],
        [ac_cuda_dir="$CUDA_HOME"]
)

ac_cuda_available=no

AC_LANG_PUSH([C])

CFLAGS_SAVED="$CFLAGS"
CUDA_CFLAGS="-I$ac_cuda_dir/include"
CFLAGS="$CUDA_CFLAGS"
export CFLAGS

CUDA_LDFLAGS="-L$ac_cuda_dir/lib64 -L$ac_cuda_dir/lib -lcudart"
LDFLAGS_SAVED="$LDFLAGS"
LDFLAGS="$CUDA_LDFLAGS"
export LDFLAGS

CC_SAVED="$CC"
CUDA_CC="$ac_cuda_dir/bin/nvcc"
CC="$CUDA_CC"
export CC

AC_CHECK_LIB([cudart],[cudaDriverGetVersion],
        [AC_CHECK_HEADER([cuda.h],[ac_cuda_available=yes],[ac_cuda_available=no],
                [#include <cuda.h>])],
        [ac_cuda_available=no],
        [$CUDA_LDFLAGS]
)

LDFLAGS="$LDFLAGS_SAVED"
export LDFLAGS
CFLAGS="$CFLAGS_SAVED"
export CFLAGS
CC="$CC_SAVED"
export CC

AC_LANG_POP([C])

if test $ac_cuda_available = yes ; then

        ifelse([$1], , , [$1])
        AC_DEFINE(HAVE_CUDA,[1],[Indicates presence of CUDA headers and lib])

        export CUDA_CFLAGS
        export CUDA_LDFLAGS
        export CUDA_CC

        AC_SUBST(CUDA_CFLAGS)
        AC_SUBST(CUDA_LDFLAGS)
        AC_SUBST(CUDA_CC)
else
        ifelse([$2], , , [$2])
fi

])

