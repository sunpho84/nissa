AC_INIT([nissa], [1.0], [francesco.sanfilippo@th.u-psud.fr])
AM_INIT_AUTOMAKE([1.0 no-define foreign])
AC_CONFIG_SRCDIR([./src])
AC_CONFIG_HEADERS([./src/config.h])

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
