AC_INIT([nissa],[1.0],[francesco.sanfilippo@th.u-psud.fr])

#overwrite default flags
if test x${ac_cv_env_CXXFLAGS_set} != xset; then
  CXXFLAGS=''
fi 

#include all config flags
config_flags="$*"
AC_DEFINE_UNQUOTED([CONFIG_FLAGS],["$config_flags"],[Flags passed to configure])

#get configuration time
AC_MSG_CHECKING([configuration date and time])
CONFIG_DATE=`date`
AC_DEFINE_UNQUOTED([CONFIG_TIME], "$CONFIG_TIME", [time when configured])
AC_MSG_RESULT([$CONFIG_DATE])

#parameters definition
AM_PROG_CC_C_O
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign])
AC_CONFIG_SRCDIR([src],[src/IO])
AC_CONFIG_HEADERS([src/config.h:src/config.h.in])

#check gcc
AC_PROG_CXX

#check library creating package
AC_PROG_RANLIB

#check yacc and lex
AC_PROG_YACC
AM_PROG_LEX

#checking for functions, headers and other stuff
#AC_CHECK_FUNCS([memset])
#AC_CHECK_FUNCS([mkdir])
#AC_CHECK_FUNCS([pow])
#AC_CHECK_FUNCS([sqrt])
#AC_CHECK_FUNCS([strcasecmp])
#AC_CHECK_FUNCS([strdup])
#AC_CHECK_FUNCS([strncasecmp])
#AC_CHECK_FUNCS([strstr])
#AC_CHECK_HEADERS([inttypes.h])
#AC_CHECK_HEADERS([libintl.h])
#AC_CHECK_HEADERS([malloc.h])
#AC_CHECK_HEADERS([stddef.h])
#AC_CHECK_HEADERS([stdint.h])
#AC_CHECK_HEADERS([stdlib.h])
#AC_CHECK_HEADERS([string.h])
#AC_CHECK_HEADERS([unistd.h])
#AC_C_INLINE
#AC_HEADER_STDBOOL
#AC_TYPE_SIZE_T
#AC_FUNC_ALLOCA
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
#AC_TYPE_INT16_T
#AC_TYPE_INT32_T
#AC_TYPE_INT8_T
#AC_TYPE_UINT16_T
#AC_TYPE_UINT32_T
#AC_TYPE_UINT64_T
#AC_TYPE_UINT8_T

# Get SVN Revision - idea from xfce-power-manager (C) 2007 by B.Tarricone
dnl get svn revision
AC_MSG_CHECKING([code revision])
SVN=`which svn`
SVN_REVISION=0
if test -n "$SVN" -a -x "$SVN"; then
    SVN_REVISION=`LC_ALL=C svn info -R "$srcdir" 2> /dev/null| awk 'BEGIN { x=0 } /^Revision: / { if($2>x) x=$2; } END { printf "%d\n",x }'`
    if test "$SVN_REVISION" -ne "0" ; then
        AC_MSG_RESULT([from svn info: $SVN_REVISION])
    fi
fi
if test "$SVN_REVISION" -eq "0" ; then
   SVN_REVISION=`if test -f svn_revision;then cat svn_revision;else echo 0;fi`
   AC_MSG_RESULT([reading from svn_revision: $SVN_REVISION])
fi
if test "$SVN_REVISION" -eq "0" ; then
   AC_MSG_RESULT([unknown, using 0])
fi
if test "$SVN_REVISION" != "0" ; then echo "$SVN_REVISION" > svn_revision;fi
AC_DEFINE_UNQUOTED([SVN_REVISION], "$SVN_REVISION", [subversion revision number])
AC_SUBST([SVN_REVISION])

#spi
AC_ARG_ENABLE(spi,
	AS_HELP_STRING([--enable-spi],[Enable spi stuff (ony for BG/Q)]),
	enable_spi="${enableval}",
	enable_spi="no")
if test "$enable_spi" == "yes";then AC_DEFINE([SPI],1,[Enable spi]) fi
AM_CONDITIONAL([SPI],[test "$enable_spi" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling spi... $enable_spi])

#MPI I/O
AC_ARG_ENABLE(mpi-io,
	AS_HELP_STRING([--disable-mpi-io],[Enable MPI I/O]),
	enable_mpi_io="${enableval}",
	enable_mpi_io="yes")
if test "$enable_mpi_io" == "yes";then AC_DEFINE([USE_MPI_IO],1,[Enable MPI I/O]) fi
AM_CONDITIONAL([USE_MPI_IO],[test "$enable_mpi_io" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling MPI I/O... $enable_mpi_io])

#threads
AC_ARG_ENABLE(threads,
	AS_HELP_STRING([--disable-threads],[Enable threads]),
	enable_threads="${enableval}",
	enable_threads="yes")
if test "$enable_threads" == "yes";then AC_DEFINE([USE_THREADS],1,[Enable threads]) fi
AM_CONDITIONAL([USE_THREADS],[test "$enable_threads" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling threads... $enable_threads])

#bgq
AC_ARG_ENABLE(bgq,
	AS_HELP_STRING([--enable-bgq],[Enable bgq stuff (ony for BG/Q)]),
	enable_bgq="${enableval}",
	enable_bgq="no")
if test "$enable_bgq" == "yes";then AC_DEFINE([BGQ],1,[Enable bgq]) AC_DEFINE([USE_VNODES],1,[Enable virtual node parallelization]) fi
AM_CONDITIONAL([BGQ],[test "$enable_bgq" == "yes" ],[true],[false])
AM_CONDITIONAL([USE_VNODES],[test "$enable_bgq" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling bgq... $enable_bgq])

#bgq emulation mode
AC_ARG_ENABLE(bgq-emu,
	AS_HELP_STRING([--enable-bgq_emu],[Enable bgq_emu stuff (ony for BG/Q)]),
	enable_bgq_emu="${enableval}",
	enable_bgq_emu="no")
if test "$enable_bgq_emu" == "yes";then AC_DEFINE([BGQ_EMU],1,[Enable bgq emulation]) fi
AM_CONDITIONAL([BGQ_EMU],[test "$enable_bgq_emu" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling bgq emulation... $enable_bgq_emu])

#repropducible run
AC_ARG_ENABLE(reproducible-run,
	AS_HELP_STRING([--enable-reproducible-run],[Enable reproducibilty of runs (reductions, etc)]),
	enable_reproducible_run="${enableval}",
	enable_reproducible_run="yes")
if test "$enable_reproducible_run" == "yes";then AC_DEFINE([REPRODUCIBLE_RUN],1,[Enable reproducible run]) fi
AM_CONDITIONAL([REPRODUCIBLE_RUN],[test "$enable_reproducible_run" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling reproducibile run... $enable_reproducible_run])

#thread debug
AC_ARG_ENABLE(thread-debug,
	AS_HELP_STRING([--enable-thread-debug],[Enable thread debugging]),
	enable_thread_debug="${enableval}",
	enable_thread_debug="no")
if test "$enable_thread_debug" == "yes";then AC_DEFINE([THREAD_DEBUG],1,[Enable thread_debug]) fi
AM_CONDITIONAL([THREAD_DEBUG],[test "$enable_thread_debug" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling thread debug... $enable_thread_debug])

#benchmark communications
AC_ARG_ENABLE(comm-bench,
	AS_HELP_STRING([--enable-comm-bench],[Enable benchmark of communications]),
	enable_comm_bench="${enableval}",
	enable_comm_bench="no")
if test "$enable_comm_bench" == "yes";then AC_DEFINE([COMM_BENCH],1,[Enable comm_bench]) fi
AM_CONDITIONAL([COMM_BENCH],[test "$enable_comm_bench" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling communication benchmark... $enable_comm_bench])

#max verbosity lv
AC_ARG_ENABLE(max-verbosity-lv,
        AS_HELP_STRING([--enable-max-verbosity-lv=lv],[Set max verbosity lv]),
        max_verbosity_lv="${enableval}",
        max_verbosity_lv=2)
if ! printf "%f" "$max_verbosity_lv" >/dev/null 2>&1; then AC_MSG_ERROR([String $max_verbosity_lv is not a number]); fi
AC_DEFINE_UNQUOTED([MAX_VERBOSITY_LV],$max_verbosity_lv,[Max_verbosity_lv])
AC_MSG_RESULT([enabling max verbosity lv... $max_verbosity_lv])

#check mpi.h
#AC_CHECK_HEADERS([mpi.h], , [AC_MSG_ERROR([MPI not found])] )

#create makefile
#AC_CONFIG_FILES([Makefile])
#AC_CONFIG_FILES([work_in_prog/bgq/Makefile])
#AC_CONFIG_FILES([work_in_prog/opt_2pts/parser/Makefile])

#create config.h
AC_OUTPUT
AC_OUTPUT(Makefile src/Makefile src/IO/Makefile)
