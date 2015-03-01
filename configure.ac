AC_INIT([nissa],[1.0],[fr.sanfilippo@gmail.com])

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
AC_CANONICAL_HOST

# Checks for programs
AC_PROG_CXX
AC_PROG_RANLIB

# Checks for libraries
AC_SEARCH_LIBS(__gmpz_init,gmp)

# Check ntl library
LB_CHECK_NTL

#eigen
AC_CHECK_HEADER([eigen3/Eigen/Dense],[AC_DEFINE([HAVE_EIGEN],[1],[Defined if Eigen found])])

#spi
AC_ARG_ENABLE(spi,
	AS_HELP_STRING([--enable-spi],[Enable spi stuff (ony for BG/Q)]),
	enable_spi="${enableval}",
	enable_spi="no")
if test "$enable_spi" == "yes";then AC_DEFINE([SPI],1,[Enable spi]) fi
AM_CONDITIONAL([SPI],[test "$enable_spi" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling spi... $enable_spi])

#high-prec arithmetic
AC_DEFINE([NATIVE_HIGH_PREC],[0],["Native implementation"])
AC_DEFINE([GMP_HIGH_PREC],[1],["GMP library"])
AC_ARG_WITH(high-prec,
	AS_HELP_STRING([--with-high-prec],[Select the high precision type (native=256 bit float, gmp [default if found])]),
	with_high_prec="${withval}",
	if test "$ac_cv_search___gmpz_init" != "no"; then with_high_prec="gmp";else with_high_prec="native";fi)
	case "$with_high_prec" in
     	     native) AC_DEFINE_UNQUOTED([HIGH_PREC],[NATIVE_HIGH_PREC],[Enable native high-precision]) ;;
     	     gmp) if test "$ac_cv_search___gmpz_init" != ""; then AC_DEFINE_UNQUOTED([HIGH_PREC],[GMP_HIGH_PREC],[Enable gmp high-precision]) AM_LDFLAGS="-lgmp -lgmpxx";else AC_MSG_ERROR(["Unable to find gmp lib"]);fi ;;
	     *) AC_MSG_ERROR(["Unkwnown high precision type ${withval}"])
	esac
AC_MSG_RESULT([with high precision type... $with_high_prec])

#cgm-debug
AC_ARG_ENABLE(cgm-debug,
	AS_HELP_STRING([--enable-cgm-debug],[Enable extremely verbose cgm inverter]),
	enable_cgm_debug="${enableval}",
	enable_cgm_debug="no")
if test "$enable_cgm_debug" == "yes";then AC_DEFINE([CGM_DEBUG],1,[Enable debugging cgm inverter]) fi
AC_MSG_RESULT([enabling cgm debug... $enable_cgm_debug])

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
if test "$enable_thread_debug" = "yes";then AC_DEFINE([THREAD_DEBUG],1,[Enable thread_debug]) fi
if test "$enable_thread_debug" = "2";then AC_DEFINE([THREAD_DEBUG],2,[Enable thread_debug]) fi
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
AC_MSG_RESULT([enabling max verbosity lv... $max_verbosity_lv])
AC_DEFINE_UNQUOTED([MAX_VERBOSITY_LV],$max_verbosity_lv,[Max_verbosity_lv])

#ndim
AC_ARG_ENABLE(ndim,
        AS_HELP_STRING([--enable-ndir=n],[Set number of dimensions]),
        ndim="${enableval}",
        ndim=4)
if ! printf "%f" "$ndim" >/dev/null 2>&1; then AC_MSG_ERROR([String $ndim is not a number]); fi
AC_MSG_RESULT([enabling number of dimension... $ndim])
AC_DEFINE_UNQUOTED([NDIM],$ndim,[Ndim])

#create config.hpp
AC_CONFIG_FILES(
	Makefile
	src/Makefile
	projects/Makefile
	projects/eight_BK/Makefile
	projects/g/Makefile
	projects/generate_confs/Makefile
	projects/ib/Makefile
	projects/semileptonic/Makefile
	projects/static_potential/Makefile
	projects/x/Makefile
	tools/Makefile)

AC_OUTPUT
