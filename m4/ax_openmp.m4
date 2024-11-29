# usage: AX_OPENMP

# overload AC_OPENMP adding a further check
# set OPENMP_CFLAGS and all derivatives, and have_openmp to yes or no

#edit from autoconf header

AC_DEFUN([AX_OPENMP], [
		      AC_REQUIRE([_AC_OPENMP_SAFE_WD])
		AC_ARG_ENABLE([openmp],
		   AS_HELP_STRING([--disable-openmp], [do not use OpenMP]))dnl

		  OPENMP_[]_AC_LANG_PREFIX[]FLAGS=
		  if test "$enable_openmp" != no; then
		    AC_CACHE_CHECK([for $[]_AC_CC[] option to support OpenMP],
		      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp],
		      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='not found'
		      for ac_option in '' -fopenmp '-Xcompiler -fopenmp' -xopenmp -openmp -mp -omp -qsmp=omp -homp \
		                       -Popenmp --openmp; do
		
		        ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
		        _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $ac_option"
		        AC_COMPILE_IFELSE([_AC_LANG_OPENMP],
		          [AC_LINK_IFELSE([_AC_LANG_OPENMP],
		            [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp=$ac_option],
		            [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='unsupported'])])
		        _AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS
		
		        if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != 'not found'; then
		          break
		        fi
		      done
		      if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" = 'not found'; then
		        ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='unsupported'
		      elif test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" = ''; then
		        ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='none needed'
		      fi
		      dnl _AC_OPENMP_SAFE_WD checked that these files did not exist before we
		      dnl started probing for OpenMP support, so if they exist now, they were
		      dnl created by the probe loop and it's safe to delete them.
		      rm -f penmp mp])
		    if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != 'unsupported' && \
		       test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != 'none needed'; then
		      OPENMP_[]_AC_LANG_PREFIX[]FLAGS="$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp"
		    fi
		  fi
		  AC_SUBST([OPENMP_]_AC_LANG_PREFIX[FLAGS])
		
		      
		      SAVE_CFLAGS="$CFLAGS"
		      SAVE_CPPFLAGS="$CPPFLAGS"
		      SAVE_CXXFLAGS="$CXXFLAGS"
		      
		      CFLAGS="$CFLAGS $OPENMP_CFLAGS"
		      CPPFLAGS="$CPPFLAGS $OPENMP_CPPFLAGS"
		      CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
		      
		      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <omp.h>]], [[ return omp_get_num_threads (); ]])],[have_openmp="yes"],[have_openmp="no"])

		      CFLAGS="$SAVE_CFLAGS"
		      CPPFLAGS="$SAVE_CPPFLAGS"
		      CXXFLAGS="$SAVE_CXXFLAGS"
		      
		      if test "$have_openmp" == yes
		      then
		      	 AC_DEFINE([HAVE_OPENMP],1,"Define to 1 if you have OpenMP support")
		      fi		      
		      
		      AC_MSG_RESULT([checking for OpenMP... ${have_openmp}])])
