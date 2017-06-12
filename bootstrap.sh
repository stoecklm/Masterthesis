#! /bin/sh

(

#-----------------------------------------------------------------------------#
### Get full paths of autotools.
ACLOCAL=`which aclocal`
AUTOCONF=`which autoconf`
AUTOHEADER=`which autoheader`
AUTOMAKE=`which automake`
#-----------------------------------------------------------------------------#
### MacOSX (and other strange) platform fixes.
if [ ! -z `which libtoolize` ]; then
   LIBTOOLIZE=`which libtoolize`
else
   if [ ! -z `which glibtoolize` ]; then
      LIBTOOLIZE=`which glibtoolize`
   fi
fi

#-----------------------------------------------------------------------------#
### Run bootstrap process tightly, so as to catch errors.
$ACLOCAL -I config/m4/acarchive -I config/m4  && \
$LIBTOOLIZE --force && \
$AUTOHEADER -Wall --force && \
$AUTOMAKE --add-missing --gnu -Wall && \
$AUTOCONF

)

