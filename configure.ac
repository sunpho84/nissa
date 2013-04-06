AC_INIT([nissa], [1.0], [francesco.sanfilippo@th.u-psud.fr])

#include all config flags
config_flags="$*"
AC_DEFINE_UNQUOTED([CONFIG_FLAGS],["$config_flags"],[Flags passed to configure])

AM_INIT_AUTOMAKE([1.0 no-define foreign])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([src/config.h:src/config.h.in])

# Get SVN Revision - idea from xfce-power-manager (C) 2007 by B.Tarricone
dnl get svn revision
AC_MSG_CHECKING([code revision])
SVN=`which svn`
SVN_REVISION=0
if test -n "$SVN" -a -x "$SVN"; then
    SVN_REVISION=`LC_ALL=C svn info -R "$srcdir" | awk 'BEGIN { x=0 } /^Revision: / { if($2>x) x=$2; } END { printf "%d\n",x }'`
    if test "$SVN_REVISION" -ne "0" ; then
        AC_MSG_RESULT([from svn info: $SVN_REVISION])
    fi
fi
if test "$SVN_REVISION" -eq 0 -a -e "$srcdir/ChangeLog"; then
    SVN_REVISION=`awk 'BEGIN { x=0 } /Revision: / { if($2>x) x=$2; } END { printf "%d\n",x }' "$srcdir/ChangeLog"`
    if test "$SVN_REVISION" -ne "0" ; then
        AC_MSG_RESULT([from ChangeLog: $SVN_REVISION])
    fi
fi
if test "$SVN_REVISION" -eq "0" ; then
    AC_MSG_RESULT([unknown, using 0])
fi
AC_DEFINE_UNQUOTED([SVN_REVISION], "$SVN_REVISION", [subversion revision number])
AC_SUBST([SVN_REVISION])

#bgq
AC_ARG_ENABLE(bgq,
	AS_HELP_STRING([--enable-bgq],[Enable bgq stuff (ony for BG/Q)]),
	enable_bgq="${enableval}",
	enable_bgq="no")
if test "$enable_bgq" == "yes";then AC_DEFINE([BGQ],1,[Enable bgq]) fi
AM_CONDITIONAL([BGQ],[test "$enable_bgq" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling bgq... $enable_bgq])

#repropducible run
AC_ARG_ENABLE(reproducible-run,
	AS_HELP_STRING([--enable-reproducible-run],[Enable reproducibilty of runs (reductions, etc)]),
	enable_reproducible_run="${enableval}",
	enable_reproducible_run="yes")
if test "$enable_reproducible_run" == "yes";then AC_DEFINE([REPRODUCIBLE_RUN],1,[Enable reproducible run]) fi
AM_CONDITIONAL([REPRODUCIBLE_RUN],[test "$enable_reproducible_run" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling reproducibile run... $enable_reproducible_run])

#overwrite default flags
if test x${ac_cv_env_CXXFLAGS_set} != xset; then
  CXXFLAGS=''
fi 

#check gcc
AC_PROG_CC
AC_PROG_CXX

#check library creating package
AC_PROG_RANLIB

#check mpi.h
#AC_CHECK_HEADERS([mpi.h], , [AC_MSG_ERROR([MPI not found])] )

#create makefile
AC_CONFIG_FILES([Makefile])

#create config.h
AC_OUTPUT
