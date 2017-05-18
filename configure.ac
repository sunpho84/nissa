AC_INIT([nissa],[1.0],[fr.sanfilippo@gmail.com])

#include all config flags
config_flags="$*"
AC_DEFINE_UNQUOTED([CONFIG_FLAGS],["$config_flags"],[Flags passed to configure])

#silent automake
AM_SILENT_RULES([yes])

#get configuration time
AC_MSG_CHECKING([configuration date and time])
CONFIG_TIME=`date`
AC_DEFINE_UNQUOTED([CONFIG_TIME], "$CONFIG_TIME", [time when configured])
AC_MSG_RESULT([$CONFIG_TIME])

#parameters definition
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

AC_LANG([C++])

#checks for programs
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

AX_CXXFLAGS_WARN_ALL

#set mpi
LX_FIND_MPI
CPPFLAGS="$AM_CPPFLAGS $MPI_CXXFLAGS $CPPFLAGS"
CXXFLAGS="$AM_CXXFLAGS $MPI_CXXFLAGS $CXXFLAGS"
CFLAGS="$AM_CFLAGS $MPI_CFLAGS $CFLAGS"
LDFLAGS="$AM_LDFLAGS $MPI_CXXLDFLAGS $LDFLAGS"
LIBS="$MPI_CXXLDFLAGS $LIBS"

#checks for headers
AC_CHECK_HEADERS([fftw3.h],[fftw_found_headers=yes])
AC_CHECK_HEADERS([gmpxx.h],[gmp_found_headers=yes])

#checks for libraries
AC_SEARCH_LIBS(__gmpz_init,gmp)
AC_SEARCH_LIBS(fftw_plan_dft_r2c_2d,fftw3)

#check ntl library
#LB_CHECK_NTL(,,AC_MSG_ERROR(["Unable to find NTL library"]))

#spi
AC_ARG_ENABLE(spi,
	AS_HELP_STRING([--enable-spi],[Enable spi stuff (ony for BG/Q)]),
	enable_spi="${enableval}",
	enable_spi="no")
if test "$enable_spi" == "yes"
then
	AC_DEFINE([SPI],1,[Enable spi])
fi
AM_CONDITIONAL([SPI],[test "$enable_spi" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling spi... $enable_spi])

#fast fourier
AC_DEFINE([NATIVE_FFT],[0],["Native implementation"])
AC_DEFINE([FFTW_FFT],[1],["FFTW library"])
AC_ARG_WITH(fft,
	AS_HELP_STRING([--with-fft],[Select the Fast Fourier transform implementation (native, fftw [NOT default if found])]),
	with_fft="${withval}",
	if test "$ac_cv_search_fftw_plan_dft_r2c_2d" != "no" -a x"$fftw_found_headers" == "xyes"
	then
		with_fft="fftw"
	else
		with_fft="native"
	fi)
case "$with_fft" in
     native) AC_DEFINE_UNQUOTED([FFT_TYPE],[NATIVE_FFT],[Enable native Fast Fourier transform]) ;;
     fftw) if test "$ac_cv_search_fftw_plan_dft_r2c_2d" != "" -a x"$fftw_found_headers" == "xyes"
    	   then
		AC_DEFINE_UNQUOTED([FFT_TYPE],[FFTW_FFT],[Enable FFTW])
		LDFLAGS="-lfftw3 -lfftw3_threads $LDFLAGS"
	   else
		AC_MSG_ERROR(["Unable to find fftw3 lib"])
     	   fi ;;
     *) AC_MSG_ERROR(["Unkwnown Fast Fourier type ${withval}"])
esac
AC_MSG_RESULT([with Fast Fourier type... ${with_fft}])

#add Eigen path, search it and enable
AC_ARG_WITH(eigen,
	AS_HELP_STRING([--with-eigen],[Specify where to find Eigen/Dense]),
	with_eigen=${withval}
	CPPFLAGS="$CPPFLAGS -I$with_eigen",
	with_eigen=no)
AC_MSG_RESULT([with Eigen... ${with_eigen}])
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([Eigen/Dense])
AC_LANG_POP([C++])
AC_ARG_ENABLE(eigen,
	AS_HELP_STRING([--enable-eigen],[Enable Eigen]),
	enable_eigen="${enableval}",
	if test "$ac_cv_header_Eigen_Dense" == "no"
	then
		enable_eigen="no"
	else
		enable_eigen="yes"
	fi)
if test "$enable_eigen" == "yes"
then
	if test "$ac_cv_header_Eigen_Dense" == "no"
	then
		AC_MSG_ERROR(["Unable to find Eigen"])
	else
		AC_DEFINE_UNQUOTED([USE_EIGEN],[1],[Enable Eigen])
	fi
fi
AC_MSG_RESULT(enabling Eigen ... ${enable_eigen})

#lime
AC_ARG_WITH(lime,
	AS_HELP_STRING([--with-lime[=dir]], [Specify where to find lime]),
	with_lime="${withval}"
	LDFLAGS="-L${with_lime}/lib/ $LDFLAGS",
	with_lime=no)
AC_SEARCH_LIBS(limeReaderNextRecord,lime)
AC_MSG_RESULT(with lime ... ${with_lime})

#threads
AC_ARG_ENABLE(threads,
	AS_HELP_STRING([--disable-threads],[Enable threads]),
	enable_threads="${enableval}",
	enable_threads="yes")
if test "$enable_threads" == "yes"
then
	#check openmp
	AC_OPENMP
	CFLAGS="$CXXFLAGS $OPENMP_CFLAGS"
	CPPFLAGS="$CPPFLAGS $OPENMP_CPPFLAGS"
	CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
	AC_DEFINE([USE_THREADS],1,[Enable threads])
fi
AM_CONDITIONAL([USE_THREADS],[test "$enable_threads" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling threads... $enable_threads])

#MPI
AC_ARG_ENABLE(MPI,
	AS_HELP_STRING([--disable-MPI],[Enable MPI]),
	enable_MPI="${enableval}",
	enable_MPI="yes")
if test "$enable_MPI" == "yes"
then
	AC_DEFINE([USE_MPI],1,[Enable MPI])
fi
AM_CONDITIONAL([USE_MPI],[test "$enable_MPI" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling MPI... $enable_MPI])

#tmLQCD
AC_ARG_WITH(tmLQCD,
	AS_HELP_STRING([--with-tmLQCD[=dir]],[Specify where to find tmLQCD]),
	with_tmlqcd="${withval}",
	with_tmlqcd="no")
AM_CONDITIONAL([USE_TMLQCD],[test "x${with_tmlqcd}" != "xno" ],[true],[false])
if test "x${with_tmlqcd}" != "xno"
then
	if test "$ac_cv_search_limeReaderNextRecord" == "no"
	then
		AC_MSG_ERROR(["Need lime to compile"])
	fi
	AC_DEFINE_UNQUOTED([USE_TMLQCD],[1],[Enable tmLQCD usage])
	LDFLAGS="-L${with_tmlqcd}/lib $LDFLAGS"
	LIBS="-lwrapper -lhmc -lmeas -lmonomial -loperator -lsolver -linit -llinalg -lhmc -lxchange -lrational -lio -lblas -llapack -llime $LIBS"
	CPPFLAGS="-I${with_tmlqcd}/include $CPPFLAGS"
fi
AC_MSG_RESULT([with tmLQCD ... ${with_tmlqcd}])

#DDalphaAMG
AC_ARG_WITH(DDalphaAMG,
	AS_HELP_STRING([--with-DDalphaAMG[=dir]],[Specify where to find DDalphaAMG]),
	with_ddalphaamg="${withval}",
	with_ddalphaamg="no")
AM_CONDITIONAL([USE_DDALPHAAMG],[test "x${with_ddalphaamg}" != "xno" ],[true],[false])
if test "x${with_ddalphaamg}" != "xno"
then
	if test "$ac_cv_search_limeReaderNextRecord" == "no"
	then
		AC_MSG_ERROR(["Need lime to compile"])
	fi
	LDFLAGS="-L${with_ddalphaamg}/lib $LDFLAGS"
	CPPFLAGS="-I${with_ddalphaamg}/include $CPPFLAGS"
	AC_CHECK_LIB(DDalphaAMG,DDalphaAMG_get_communicator,,AC_MSG_ERROR(["Unable to find libDDalphaAMG.a"]))
	AC_CHECK_HEADERS([DDalphaAMG.h],,AC_MSG_ERROR(["Unable to find DDalphaAMG.h"]))
	AC_DEFINE_UNQUOTED([USE_DDALPHAAMG],[1],[Enable DDalphaAMG usage])
fi
AC_MSG_RESULT([with DDalphaAMG ... ${with_ddalphaamg}])

#high-prec arithmetic
AC_DEFINE([NATIVE_HIGH_PREC],[0],["Native implementation"])
AC_DEFINE([GMP_HIGH_PREC],[1],["GMP library"])
AC_ARG_WITH(high-prec,
	AS_HELP_STRING([--with-high-prec],[Select the high precision type (native=256 bit float, gmp [default if found])]),
	with_high_prec="${withval}",
	if test "$ac_cv_search___gmpz_init" != "no"; then with_high_prec="gmp";else with_high_prec="native";fi)
	case "$with_high_prec" in
     	     native) AC_DEFINE_UNQUOTED([HIGH_PREC],[NATIVE_HIGH_PREC],[Enable native high-precision]) ;;
     	     gmp) if test "$ac_cv_search___gmpz_init" != "" -a  x"$gmp_found_headers" == "xyes"; then AC_DEFINE_UNQUOTED([HIGH_PREC],[GMP_HIGH_PREC],[Enable gmp high-precision]) AM_LDFLAGS="-lgmp -lgmpxx";else AC_MSG_ERROR(["Unable to find gmp lib"]);fi ;;
	     *) AC_MSG_ERROR(["Unkwnown high precision type ${withval}"])
	esac
AC_MSG_RESULT([with high precision type... $with_high_prec])

#cgm-debug
AC_ARG_ENABLE(cgm-debug,
	AS_HELP_STRING([--enable-cgm-debug],[Enable extremely verbose cgm inverter]),
	enable_cgm_debug="${enableval}",
	enable_cgm_debug="no")
if test "$enable_cgm_debug" == "yes"
then
	AC_DEFINE([CGM_DEBUG],1,[Enable debugging cgm inverter])
fi
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
if test "$enable_mpi_io" == "yes"
then
	AC_DEFINE([USE_MPI_IO],1,[Enable MPI I/O])
fi
AM_CONDITIONAL([USE_MPI_IO],[test "$enable_mpi_io" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling MPI I/O... $enable_mpi_io])

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
	if test "$enable_bgq" == "yes"
	then
		enable_bgq_emu="no"
	else
	enable_bgq_emu="yes"
	fi)
if test "$enable_bgq_emu" == "yes"
then
	AC_DEFINE([BGQ_EMU],1,[Enable bgq emulation])
fi
AM_CONDITIONAL([BGQ_EMU],[test "$enable_bgq_emu" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling bgq emulation... $enable_bgq_emu])

#repropducible run
AC_ARG_ENABLE(reproducible-run,
	AS_HELP_STRING([--enable-reproducible-run],[Enable reproducibilty of runs (default: no)]),
	enable_reproducible_run="${enableval}",
	enable_reproducible_run="no")
if test "$enable_reproducible_run" == "yes"
then
   	AC_DEFINE([REPRODUCIBLE_RUN],1,[Enable reproducible run])
fi
AM_CONDITIONAL([REPRODUCIBLE_RUN],[test "$enable_reproducible_run" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling reproducibile run... $enable_reproducible_run])

#thread debug
AC_ARG_ENABLE(thread-debug,
	AS_HELP_STRING([--enable-thread-debug],[Enable thread debugging]),
	enable_thread_debug="${enableval}",
	enable_thread_debug="no")
if test "$enable_thread_debug" = "yes"
then
	AC_DEFINE([THREAD_DEBUG],1,[Enable thread_debug])
	fi
if test "$enable_thread_debug" = "2"
then
	AC_DEFINE([THREAD_DEBUG],2,[Enable thread_debug])
fi
AM_CONDITIONAL([THREAD_DEBUG],[test "$enable_thread_debug" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling thread debug... $enable_thread_debug])

#max verbosity lv
AC_ARG_ENABLE(max-verbosity-lv,
        AS_HELP_STRING([--enable-max-verbosity-lv=lv],[Set max verbosity lv]),
        max_verbosity_lv="${enableval}",
        max_verbosity_lv=2)
if ! printf "%f" "$max_verbosity_lv" >/dev/null 2>&1; then AC_MSG_ERROR([String $max_verbosity_lv is not a number]); fi
AC_DEFINE_UNQUOTED([MAX_VERBOSITY_LV],$max_verbosity_lv,[Max_verbosity_lv])
AC_MSG_RESULT([enabling max verbosity lv... $max_verbosity_lv])

#number of dimension
AC_ARG_ENABLE(ndim,
        AS_HELP_STRING([--enable-ndim=n],[Set number of dimensions [Default: 4]]),
        ndim="${enableval}",
        ndim=4)
if ! printf "%f" "$ndim" >/dev/null 2>&1; then AC_MSG_ERROR([String $ndim is not a number]); fi
AM_CONDITIONAL([NDIM2], [test x$ndim = x2])
AC_DEFINE_UNQUOTED([NDIM],$ndim,[Number of dimensions])
AC_MSG_RESULT([enabling number of dimension... $ndim])

#number of colors
AC_ARG_ENABLE(ncol,
        AS_HELP_STRING([--enable-ncol=n],[Set number of colors [Default: 3]]),
        ncol="${enableval}",
        ncol=3)
if ! printf "%f" "$col" >/dev/null 2>&1; then AC_MSG_ERROR([String $ncol is not a number]); fi
AM_CONDITIONAL([NCOL3], [test x$ncol = x3])
AC_DEFINE_UNQUOTED([NCOL],$ncol,[Number of colors])
AC_MSG_RESULT([enabling number of colors... $ncol])

#macros to concatenate
AC_DEFINE([_CONCAT(X,Y)],[X##Y],[Internally concatenation])
AC_DEFINE([CONCAT(X,Y)],[_CONCAT(X,Y)],[Concatenate two pieces to produce a new token])
AC_DEFINE([CONCAT2(s1,s2)],[CONCAT(s1,s2)],[Wrapper to beat CPP])
AC_DEFINE([CONCAT3(s1,s2,s3)],[CONCAT(CONCAT2(s1,s2),s3)],[Concatenate three])
AC_DEFINE([NAME2(s1,s2)],[CONCAT3(s1,_,s2)],[Link with a _])
AC_DEFINE([NAME3(s1,s2,s3)],[NAME2(CONCAT3(s1,_,s2),s3)],[Name with two _])
AC_DEFINE([NAME4(s1,s2,s3,s4)],[NAME3(CONCAT3(s1,_,s2),s3,s4)],[Name with four _])
AC_DEFINE([REORDER_BARRIER()],[__asm volatile ("")],[Avoid spilling])
AC_DEFINE([MACRO_GUARD(...)],[do{__VA_ARGS__}while(0)],[Guard semicolon])

#create config.hpp
AC_CONFIG_FILES(
	Makefile
	src/Makefile
	projects/Makefile
	projects/cpn/Makefile
	projects/generate_confs/Makefile
	projects/ib/Makefile
	projects/nucleons/Makefile
	projects/reno_const/Makefile
	projects/semileptonic/Makefile
	tools/Makefile)

AC_OUTPUT
