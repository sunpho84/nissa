# Check for NTL
# Bradford Hovinen, 2001-06-13
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_NTL ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for Victor Shoup's NTL (Number Theory Library) and define
dnl NTL_CFLAGS and NTL_LIBS

AC_DEFUN([LB_CHECK_NTL],
[

AC_ARG_WITH(ntl-prefix,[  --with-ntl-prefix=PFX   Prefix where NTL is installed (optional)],
[ntl_prefix="$withval"],[ntl_prefix=""])

min_ntl_version=ifelse([$1], ,4.0,$1)
AC_MSG_CHECKING(for NTL >= $min_ntl_version)

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LDFLAGS=${LDFLAGS}
BACKUP_LIBS=${LIBS}

if test "$ntl_prefix" != ""; then
   CXXFLAGS+="-I${ntl_prefix}/include"
   LDFLAGS+="-L${ntl_prefix}"
fi
LIBS+=" -lntl"

dnl Check for existence

AC_TRY_LINK(
[#include <NTL/ZZ.h>
using namespace NTL;],
[ZZ a;],
[
AC_TRY_RUN(
[#include <NTL/version.h>
#include <iostream>
int main () { if (NTL_MAJOR_VERSION < 4) return -1; else return 0; }
],[
AC_MSG_RESULT(found)
AC_DEFINE([HAVE_NTL], [], [NTL found])

ifelse([$2], , :, [$2])
],[
AC_MSG_RESULT(not found)
echo "Sorry, your NTL version is too old. Disabling."

CXXFLAGS=${BACKUP_CXXFLAGS}
LDFLAGS=${BACKUP_LDFLAGS}
LIBS=${BACKUP_LIBS}

ifelse([$3], , :, [$3])
])
],
[
AC_MSG_RESULT(not found)
if test "$ntl_prefix" != ""; then
	AC_MSG_WARN(NTL >= 4.0 was not found. Please double-check the directory you gave.)
fi

unset NTL_CXXFLAGS
unset NTL_LDFLAGS
unset NTL_LIBS

CXXFLAGS=${BACKUP_CXXFLAGS}
LDFLAGS=${BACKUP_LDFLAGS}
LIBS=${BACKUP_LIBS}

ifelse([$3], , :, [$3])
])

])
