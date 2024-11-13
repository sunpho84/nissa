# usage: AX_OPENMP

# overload AC_OPENMP adding a further check
# set OPENMP_CFLAGS and all derivatives, and have_openmp to yes or no


AC_DEFUN([AX_OPENMP], [
		      AC_OPENMP
		      
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
