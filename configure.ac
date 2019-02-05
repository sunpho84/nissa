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
AM_INIT_AUTOMAKE([1.14 no-define foreign subdir-objects])
AM_CONFIG_HEADER(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

#checks for programs
AC_PROG_CC
AC_PROG_FC
AC_PROG_CXX
AC_PROG_RANLIB
AC_LANG([C++])

LIBRARY_RESULT="---"
SUMMARY_RESULT="---"

#check if we enable fortran
AC_ARG_ENABLE(fortran,
	AS_HELP_STRING([--enable-fortran],[Enable fortran support]),
	enable_fortran="${enableval}",
	enable_fortran="no")
AC_MSG_RESULT([enabling fortran... $enable_fortran])

#check if we enable splitting MPI autosearch
AC_ARG_ENABLE(MPI-autosearch,
	AS_HELP_STRING([--enable-MPI-autosearch],[Enable MPI autosearch]),
	enable_MPI_autosearch="${enableval}",
	enable_MPI_autosearch="yes")
AC_MSG_RESULT([enabling MPI autosearch... $enable_MPI_autosearch])

#split fortran flags into LDFLAGS and LIBS
if test "$enable_fortran" == "yes"
then
	AC_F77_LIBRARY_LDFLAGS
	fortran_ldflags=` echo "$FLIBS" | grep -o -- '\(^\| \)-L\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)' | tr '\n' ' ' | sed 's/^[[ \t]]*//;s/[[ \t]]*$//' | sed 's/  +/ /g'`
	fortran_libs=`       echo "$FLIBS" | grep -o -- '\(^\| \)-l\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)' | tr '\n' ' ' | sed 's/^[[ \t]]*//;s/[[ \t]]*$//' | sed 's/  +/ /g'`
	LDFLAGS="$fortran_ldflags $LDFLAGS"
	LIBS="$fortran_libs $LIBS"
	FLIBS=""
fi

#AX_CXXFLAGS_WARN_ALL
AC_MSG_CHECKING(whether compiler understands -Wall)
OLD_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="$CXXFLAGS -Wall"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[])],
  AC_MSG_RESULT(yes),
  AC_MSG_RESULT(no)
  CXXFLAGS="$OLD_CXXFLAGS")

#blas and lapack
if test "$enable_fortran" == "yes"
then
	AX_BLAS([blas_found=yes],[blas_found=no])
	AX_LAPACK([lapack_found=yes],[lapack_found=no])
	LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS"
	LIBRARY_RESULT="$LIBRARY_RESULT
	blas                : $BLAS_LIBS
	lapack              : $LAPACK_LIBS"
fi

#set MPI for fortran
if test "$enable_fortran" == "yes" -a "$enable_MPI_autosearch" == "yes"
then
	AC_LANG_PUSH([Fortran])
	LX_FIND_MPI
	AC_LANG_POP([Fortran])
	AC_MSG_RESULT([with Fortran MPI link options... $MPI_FLDFLAGS])
	LIBS="$MPI_FLIBS $LIBS"
	LDFLAGS="$MPI_FLDFLAGS $LDFLAGS"
fi

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
SUMMARY_RESULT="$SUMMARY_RESULT
Threads enabled     : $enable_threads"

#set MPI for c++
if test "$enable_MPI_autosearch" == "yes"
then
	AC_LANG_PUSH([C++])
	LX_FIND_MPI
	AC_LANG_POP([C++])
	AC_MSG_RESULT([with C++ MPI link options... $MPI_CXXLDFLAGS])
fi

#MPI
AC_ARG_ENABLE(MPI,
	AS_HELP_STRING([--disable-MPI],[Enable MPI]),
	enable_MPI="${enableval}",
	enable_MPI="yes")
AM_CONDITIONAL([USE_MPI],[test "$enable_MPI" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling MPI... $enable_MPI])
SUMMARY_RESULT="$SUMMARY_RESULT
MPI enabled         : $enable_MPI"

if test "$enable_MPI" == "yes"
then
	# Check if we have MPI
	if test "$have_CXX_mpi" == "yes"
	then
		CPPFLAGS="$MPI_CXXFLAGS $CPPFLAGS"
		LIBS="$MPI_CXXLIBS $LIBS"
		LDFLAGS="$MPI_CXXLDFLAGS $LDFLAGS"
		AC_DEFINE([USE_MPI],1,[Enable MPI])
	else
		AC_MSG_ERROR(["Unable to find MPI"])
	fi
fi

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

##############################################

AX_SUBPACKAGE(ssl,openssl/ssl.h,crypto ssl,MD5_Final,SSL)

#################### fftw ####################

AX_SUBPACKAGE(fftw,fftw3.h,fftw3,fftw_plan_dft_r2c_2d,FFTW)

#fast fourier type
AC_DEFINE([NATIVE_FFT],[0],["Native implementation"])
AC_DEFINE([FFTW_FFT],[1],["FFTW library"])
AC_ARG_WITH(fft-type,
	AS_HELP_STRING([--with-fft-type],[Select the Fast Fourier transform implementation (native, fftw [NOT default if found])]),
	with_fft_type="${withval}",
	if test "$enable_fftw" != "no"
	then
		with_fft_type="fftw"
	else
		with_fft_type="native"
	fi)
case "$with_fft_type" in
     native) AC_DEFINE_UNQUOTED([FFT_TYPE],[NATIVE_FFT],[Enable native Fast Fourier transform]) ;;
     fftw) if test "$fftw_found" != "no"
    	   then
		AC_DEFINE_UNQUOTED([FFT_TYPE],[FFTW_FFT],[Enable FFTW])
	   else
		AC_MSG_ERROR(["Unable to find fftw3 lib"])
     	   fi ;;
     *) AC_MSG_ERROR(["Unkwnown Fast Fourier type ${withval}"])
esac
AC_MSG_RESULT([with Fast Fourier type... ${with_fft_type}])
SUMMARY_RESULT="$SUMMARY_RESULT
FFT type            : $with_fft_type"

#################### eigen ####################

AX_SUBPACKAGE(eigen,eigen3/Eigen/Dense,,,EIGEN)

#enable eigen for linear algebra
AC_ARG_ENABLE(eigen_everywhere,
	AS_HELP_STRING([--enable-eigen-everywhere],[Enable Eigen Everywhere (default: no)]),
	enable_eigen_everywhere="${enableval}",
	enable_eigen_everywhere="no")
AC_MSG_RESULT(enabling Eigen everywhere ... ${enable_eigen_everywhere})

#check activability
if test "$enable_eigen_everywhere" == "yes"
then
	if test "$eigen_found" == "no"
	then
		AC_MSG_ERROR(["Unable to find Eigen"])
	else
		AC_DEFINE_UNQUOTED([USE_EIGEN_EVERYWHERE],[1],[Enable Eigen everywhere])
	fi
fi

#################### parpack ####################

AX_SUBPACKAGE(parpack,arpack/arpack.h arpack/parpack.h,arpack parpack,cnaupd_c pcnaupd_c,PARPACK)

#parpack-debug
AC_ARG_ENABLE(parpack-debug,
	AS_HELP_STRING([--enable-parpack-debug],[Enable parpack debug]),
	enable_parpack_debug="${enableval}",
	enable_parpack_debug="no")
if test "$enable_parpack_debug" == "yes"
then
	AC_DEFINE([ENABLE_PARPACK_DEBUG],1,[Enable debugging parpack])
fi
AC_MSG_RESULT([enabling parpack debug... $enable_parpack_debug])

#################### lime ####################

AX_SUBPACKAGE(lime,lime.h,lime,limeReaderNextRecord,LIME)

#################### tmLQCD ####################

#search for the library
TMLQCD_LIBS="-lwrapper -lhmc -lmeas -lmonomial -loperator -lsolver -linit -llinalg -lhmc -lxchange -lrational -lio"
OLD_LIBS=$LIBS
LIBS="$TMLQCD_LIBS $LIBS"
AX_SUBPACKAGE(tmLQCD,,tmLQCD,tmLQCD_invert_init,TMLQCD)
LIBS=$OLD_LIBS

if test "$enable_tmLQCD" == "yes"
then
	LIBS="$LIBS $TMLQCD_LIBS"
fi

#################### DDalphaAMG ####################

AX_SUBPACKAGE(DDalphaAMG,DDalphaAMG.h,DDalphaAMG,DDalphaAMG_get_communicator,DDALPHAAMG)

#################### gmp ####################

AX_SUBPACKAGE(gmp,gmp.h,gmp,__gmpz_init,GMP)

#high-prec arithmetic
AC_DEFINE([NATIVE_HIGH_PREC],[0],["Native implementation"])
AC_DEFINE([GMP_HIGH_PREC],[1],["GMP library"])
AC_ARG_WITH(high-prec,
	AS_HELP_STRING([--with-high-prec],[Select the high precision type (native=256 bit float, gmp [default if found])]),
	with_high_prec_type="${withval}",
	if test "$gmp_found" != "no"
	then
		with_high_prec_type="gmp"
	else
		with_high_prec_type="native"
	fi)
case "$with_high_prec_type" in
     native) AC_DEFINE_UNQUOTED([HIGH_PREC_TYPE],[NATIVE_HIGH_PREC],[Enable native high-precision]) ;;
     gmp) if test "$gmp_found" != "no"
	  then
		AC_DEFINE_UNQUOTED([HIGH_PREC_TYPE],[GMP_HIGH_PREC],[Enable gmp high-precision])
	  else
		AC_MSG_ERROR(["Unable to find gmp lib"]);
	  fi ;;
     *) AC_MSG_ERROR(["Unkwnown high precision type ${withval}"])
esac
AC_MSG_RESULT([with high precision type... $with_high_prec_type])
SUMMARY_RESULT="$SUMMARY_RESULT
High prec type      : $with_high_prec_type"

#################### ####################

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

#huepages
AC_ARG_ENABLE(hugepages,
	AS_HELP_STRING([--enable-hugepages],[Enable usage of hugepages]),
	enable_hugepages="${enableval}",
	enable_hugepages="no")
if test "$enable_hugepages" == "yes"
then
	AC_DEFINE([USE_HUGEPAGES],1,[Enable hugepages])
fi
AM_CONDITIONAL([USE_HUGEPAGES],[test "$enable_hugepags" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling hugepags... $enable_hugepages])
SUMMARY_RESULT="$SUMMARY_RESULT
Enable hugepages    : $enable_hugepages"

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
SUMMARY_RESULT="$SUMMARY_RESULT
Number of colors    : $ncol"

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

#result of conf to screen
AC_MSG_RESULT([
--------------------------------------------------
Configuration summary for $PACKAGE_NAME ($VERSION)
--------------------------------------------------
CC                  : $CC
CPPFLAGS            : $CPPFLAGS
CXXFLAGS            : $CXXFLAGS
CXX                 : $CXX
---
LDFLAGS             : $LDFLAGS
LIBS                : $LIBS
---
Installation prefix : $prefix
$LIBRARY_RESULT
$SUMMARY_RESULT
--------------------------------------------------
Configuration OK
--------------------------------------------------
])

#create config.hpp
AC_CONFIG_FILES(
	Makefile
	src/Makefile
	projects/Makefile
	tools/Makefile)

AC_OUTPUT
