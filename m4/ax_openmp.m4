# usage: AX_OPENMP

# overload AC_OPENMP adding a further check
# set OPENMP_CFLAGS and all derivatives, and have_openmp to yes or no


AC_DEFUN([AX_OPENMP], [
OPENMP_[]_AC_LANG_PREFIX[]FLAGS=
  AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--disable-openmp], [do not use OpenMP])])
  if test "$enable_openmp" != no; then
    AC_CACHE_CHECK([for $CC option to support OpenMP],
      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp],
      [AC_LINK_IFELSE([_AC_LANG_OPENMP],
         [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='none needed'],
         [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='unsupported'
          dnl Try these flags:
          dnl   GCC >= 4.2           -fopenmp
          dnl   SunPRO C             -xopenmp
          dnl   Intel C              -openmp
          dnl   SGI C, PGI C         -mp
          dnl   Tru64 Compaq C       -omp
          dnl   IBM C (AIX, Linux)   -qsmp=omp
          dnl If in this loop a compiler is passed an option that it doesn't
          dnl understand or that it misinterprets, the AC_LINK_IFELSE test
          dnl will fail (since we know that it failed without the option),
          dnl therefore the loop will continue searching for an option, and
          dnl no output file called 'penmp' or 'mp' is created.
          for ac_option in -fopenmp -xopenmp -openmp -mp -omp -qsmp=omp "-Xcompiler -fopenmp"; do
            ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
            _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $ac_option"
            AC_LINK_IFELSE([_AC_LANG_OPENMP],
              [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp=$ac_option])
            _AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS
            if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != unsupported; then
              break
            fi
          done])])
    case $ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp in #(
      "none needed" | unsupported)
        ;; #(
      *)
        OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp ;;
    esac
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
		      
		      if test "$enable_openmp" != no
		      then
			if test "$have_openmp" != yes
			then
				AC_MSG_ERROR(["Asked to enable OpenMP when the compiler is not supporting it"])
			fi
		      	 AC_DEFINE([USE_OPENMP],1,"Define to 1 if enabling OpenMP support")
			 CFLAGS="$CFLAGS $OPENMP_CFLAGS"
			 CPPFLAGS="$CPPFLAGS $OPENMP_CPPFLAGS"
			 CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
		      fi		      
		      
		      AC_MSG_RESULT([checking for OpenMP... ${have_openmp}])])
