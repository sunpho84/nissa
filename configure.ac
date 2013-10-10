AC_INIT([nissa],[1.0],[francesco.sanfilippo@th.u-psud.fr])

#include all config flags
config_flags="$*"
AC_DEFINE_UNQUOTED([CONFIG_FLAGS],["$config_flags"],[Flags passed to configure])

#get configuration time
AC_MSG_CHECKING([configuration date and time])
CONFIG_DATE=`date`
AC_DEFINE_UNQUOTED([CONFIG_TIME], "$CONFIG_TIME", [time when configured])
AC_MSG_RESULT([$CONFIG_DATE])

#parameters definition
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_DISABLE_SHARED
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_MKDIR_P
AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.
AC_FUNC_ALLOCA
AC_CHECK_HEADERS([float.h inttypes.h libintl.h malloc.h stddef.h stdint.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT8_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset mkdir pow sqrt strcasecmp strdup strncasecmp strndup strstr])

#cuda bins
ARCH=`uname -m`
if [[ $ARCH == "x86_64" ]];
then
  CUDA_SUFFIX="64"
else
  CUDA_SUFFIX=""
fi

#cuda
AC_ARG_ENABLE(cuda,
	AS_HELP_STRING([--enable-cuda],[Enable cuda]),
	enable_cuda="${enableval}",
	enable_cuda="no")
if test "$enable_cuda" == "yes";then 
AC_DEFINE([CUDA],1,[Enable cuda]) 
AX_CHECK_CUDA
fi
AM_CONDITIONAL([CUDA],[test "$enable_cuda" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling cuda... $enable_cuda])

#spi
AC_ARG_ENABLE(spi,
	AS_HELP_STRING([--enable-spi],[Enable spi stuff (ony for BG/Q)]),
	enable_spi="${enableval}",
	enable_spi="no")
if test "$enable_spi" == "yes";then AC_DEFINE([SPI],1,[Enable spi]) fi
AM_CONDITIONAL([SPI],[test "$enable_spi" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling spi... $enable_spi])

#parser generation
AC_ARG_ENABLE(parser-generation,
	AS_HELP_STRING([--enable-parser-generation],[Enable generation of parsers]),
	enable_parser_generation="${enableval}",
	enable_parser_generation="yes")
	if test "$enable_parser_generation" == "yes"
	then
		#check yacc and lex
		AC_PROG_YACC
		AM_PROG_LEX
	fi
AM_CONDITIONAL([PARSER_GENERATION],[test "$enable_parser_generation" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling parser generation... $enable_parser_generation])

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

#MPI
AC_ARG_ENABLE(MPI,
	AS_HELP_STRING([--disable-MPI],[Enable MPI]),
	enable_MPI="${enableval}",
	enable_MPI="yes")
if test "$enable_MPI" == "yes";then AC_DEFINE([USE_MPI],1,[Enable MPI]) fi
AM_CONDITIONAL([USE_MPI],[test "$enable_MPI" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling MPI... $enable_MPI])

#bgq
AC_ARG_ENABLE(bgq,
	AS_HELP_STRING([--enable-bgq],[Enable bgq stuff (ony for BG/Q)]),
	enable_bgq="${enableval}",
	enable_bgq="no")
if test "$enable_bgq" == "yes"
then
	 AC_DEFINE([BGQ],1,[Enable bgq])
	 AC_DEFINE([USE_VNODES],1,[Enable virtual node parallelization])
fi
AM_CONDITIONAL([BGQ],[test "$enable_bgq" == "yes" ],[true],[false])
AM_CONDITIONAL([USE_VNODES],[test "$enable_bgq" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling bgq... $enable_bgq])

#bgq emulation mode
AC_ARG_ENABLE(bgq-emu,
	AS_HELP_STRING([--enable-bgq_emu],[Enable bgq_emu stuff (ony for BG/Q)]),
	enable_bgq_emu="${enableval}",
	enable_bgq_emu="no")
if test "$enable_bgq_emu" == "yes"
then
	AC_DEFINE([BGQ_EMU],1,[Enable bgq emulation])
fi
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

#benchmarks
AC_ARG_ENABLE(bench,
	AS_HELP_STRING([--enable-bench],[Enable internal benchmarks]),
	enable_bench="${enableval}",
	enable_bench="yes")
if test "$enable_bench" == "yes";then AC_DEFINE([BENCH],1,[Enable bench]) fi
AM_CONDITIONAL([BENCH],[test "$enable_bench" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling internal benchmarks... $enable_bench])

#max verbosity lv
AC_ARG_ENABLE(max-verbosity-lv,
        AS_HELP_STRING([--enable-max-verbosity-lv=lv],[Set max verbosity lv]),
        max_verbosity_lv="${enableval}",
        max_verbosity_lv=2)
if ! printf "%f" "$max_verbosity_lv" >/dev/null 2>&1; then AC_MSG_ERROR([String $max_verbosity_lv is not a number]); fi
AC_DEFINE_UNQUOTED([MAX_VERBOSITY_LV],$max_verbosity_lv,[Max_verbosity_lv])
AC_MSG_RESULT([enabling max verbosity lv... $max_verbosity_lv])

#create config.hpp
AC_CONFIG_FILES(
	Makefile
	src/Makefile
	src/base/Makefile
	src/bgq/Makefile
	src/cuda/Makefile
	src/communicate/Makefile
	src/dirac_operators/Makefile
	src/geometry/Makefile
	src/hmc/Makefile
	src/inverters/Makefile
	src/io/Makefile
	src/linalgs/Makefile
	src/new_types/Makefile
	src/operations/Makefile
	src/routines/Makefile
	projects/Makefile
	projects/generate_confs/Makefile
	projects/semileptonic/Makefile)
AC_OUTPUT
