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

#checks for programs
AC_PROG_CC
AC_PROG_FC
AC_PROG_CXX
AC_PROG_RANLIB
AC_LANG([C++])

LIBRARY_RESULT="---"
SUMMARY_RESULT="---"

#split fortran flags into LDFLAGS and LIBS
AC_F77_LIBRARY_LDFLAGS
fortran_ldflags=` echo "$FLIBS" | grep -o -- '\(^\| \)-L\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)' | tr '\n' ' ' | sed 's/^[[ \t]]*//;s/[[ \t]]*$//' | sed 's/  +/ /g'`
fortran_libs=`       echo "$FLIBS" | grep -o -- '\(^\| \)-l\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)' | tr '\n' ' ' | sed 's/^[[ \t]]*//;s/[[ \t]]*$//' | sed 's/  +/ /g'`
LDFLAGS="$fortran_ldflags $LDFLAGS"
LIBS="$fortran_libs $LIBS"
FLIBS=""

#AX_CXXFLAGS_WARN_ALL
AC_MSG_CHECKING(whether compiler understands -Wall)
OLD_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="$CXXFLAGS -Wall"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[])],
  AC_MSG_RESULT(yes),
  AC_MSG_RESULT(no)
  CXXFLAGS="$OLD_CXXFLAGS")

#blas and lapack
AX_BLAS([blas_found=yes],[blas_found=no])
AX_LAPACK([lapack_found=yes],[lapack_found=no])
LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS"
LIBRARY_RESULT="$LIBRARY_RESULT
BLAS                : $BLAS_LIBS
LAPACK              : $LAPACK_LIBS"

#set mpi for fortran
AC_LANG_PUSH([Fortran])
LX_FIND_MPI
AC_LANG_POP([Fortran])
AC_MSG_RESULT([with Fortran MPI link options... $MPI_FLDFLAGS])
LIBS="$MPI_FLIBS $LIBS"
LDFLAGS="$MPI_FLDFLAGS $LDFLAGS"

#set mpi for c++
AC_LANG_PUSH([C++])
LX_FIND_MPI
AC_LANG_POP([C++])
AC_MSG_RESULT([with C++ MPI link options... $MPI_CXXLDFLAGS])
CPPFLAGS="$MPI_CXXFLAGS $CPPFLAGS"
LIBS="$MPI_CXXLIBS $LIBS"
LDFLAGS="$MPI_CXXLDFLAGS $LDFLAGS"

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

#################### fftw ####################

#fftw parameters
AC_ARG_WITH(fftw,
	AS_HELP_STRING([--with-fftw],[Specify where to find installation of fftw]),
	with_fftw=${withval}
	CPPFLAGS="$CPPFLAGS -I$with_fftw/include"
	LFFLAGS="$LDFLAGS -I$with_fftw/lib",
	with_fftw=no)
AC_CHECK_HEADERS([fftw3.h],[fftw_found_header=yes],[fftw_found_header=no])
AC_SEARCH_LIBS([fftw_plan_dft_r2c_2d],[fftw3],[fftw_found_library=yes],[fftw_found_library=no])
if test "$fftw_found_header" != "no"  -a "$fftw_found_library" != "no"
then
	fftw_found=yes
else
	fftw_found=no
fi
LIBRARY_RESULT="$LIBRARY_RESULT
FFTW                : $ac_cv_search_fftw_plan_dft_r2c_2d"

#fast fourier type
AC_DEFINE([NATIVE_FFT],[0],["Native implementation"])
AC_DEFINE([FFTW_FFT],[1],["FFTW library"])
AC_ARG_WITH(fft-type,
	AS_HELP_STRING([--with-fft-type],[Select the Fast Fourier transform implementation (native, fftw [NOT default if found])]),
	with_fft_type="${withval}",
	if test "$fftw_found" != "no"
	then
		with_fft_type="fftw"
	else
		with_fft_type="native"
	fi)
case "$with_fft_type" in
     native) AC_DEFINE_UNQUOTED([WITH_FFT_TYPE],[NATIVE_FFT],[Enable native Fast Fourier transform]) ;;
     fftw) if test "$fftw_found" != "no"
    	   then
		AC_DEFINE_UNQUOTED([WITH_FFT_TYPE],[FFTW_FFT],[Enable FFTW])
	   else
		AC_MSG_ERROR(["Unable to find fftw3 lib"])
     	   fi ;;
     *) AC_MSG_ERROR(["Unkwnown Fast Fourier type ${withval}"])
esac
AC_MSG_RESULT([with Fast Fourier type... ${with_fft_type}])
SUMMARY_RESULT="$SUMMARY_RESULT
FFT type            : $with_fft_type"

#################### eigen ####################

#introduce eigen flags
AC_ARG_WITH(eigen,
	AS_HELP_STRING([--with-eigen],[Specify where to find Eigen/Dense]),
	with_eigen=${withval}
	CPPFLAGS="$CPPFLAGS -I$with_eigen",
	with_eigen=no)
AC_MSG_RESULT([with Eigen... ${with_eigen}])

#search eigen header
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([Eigen/Dense],[eigen_header_found=yes],[eigen_header_found=no])
AC_LANG_POP([C++])
eigen_found=$eigen_header_found

#introdue enable eigen
AC_ARG_ENABLE(eigen,
	AS_HELP_STRING([--enable-eigen],[Enable Eigen (automatically enabled if found)]),
	enable_eigen="${enableval}",
	enable_eigen="${eigen_found}")
AC_MSG_RESULT([enabling eigen ... ${enable_eigen}])

#check activability for eigen
if test "$enable_eigen" == "yes"
then
	if test "$eigen_found" == "no"
	then
		AC_MSG_ERROR(["Cannot enable eigen, library not found, please provide/improve hint using --with-eigen flag"])
	fi
	AC_DEFINE([USE_EIGEN],1,[Enable eigen])
fi
AM_CONDITIONAL([USE_EIGEN],[test "$enable_eigen" == "yes" ],[true],[false])
SUMMARY_RESULT="$SUMMARY_RESULT
Eigen enabled       : $enable_eigen"

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

#introduce parpack flags
AC_ARG_WITH(parpack,
	AS_HELP_STRING([--with-parpack[=dir]], [Specify where to find parpack]),
	with_parpack="${withval}"
	CPPFLAGS="-I${with_parpack}/include/ $CPPFLAGS"
	LDFLAGS="-L${with_parpack}/lib/ $LDFLAGS",
	with_parpack=no)
AC_MSG_RESULT(with parpack ... ${with_parpack})

#search for arpack and parpack library
AC_CHECK_HEADERS([arpack/parpack.h],[parpack_found_headers=yes],[parpack_found_headers=no])
AC_SEARCH_LIBS([cnaupd_c],[arpack],[arpack_found_library=yes],[arpack_found_library=no])
AC_SEARCH_LIBS([pcnaupd_c],[parpack],[parpack_found_library=yes],[parpack_found_library=no])
if test "$parpack_found_headers" != "no" -a "$arpack_found_library" != "no" -a "$parpack_found_library" != "no"
then
	parpack_found=yes
else
	parpack_found=no
fi

#result of the search
LIBRARY_RESULT="$LIBRARY_RESULT
ARPACK              : $ac_cv_search_cnaupd_c
PARPACK             : $ac_cv_search_pcnaupd_c"

#introduce enable parpack
AC_ARG_ENABLE(parpack,
	AS_HELP_STRING([--enable-parpack],[Enable parpack (automatically enabled if found)]),
	enable_parpack="${enableval}",
	enable_parpack="${parpack_found}")
AC_MSG_RESULT([enabling parpack ... ${enable_parpack}])

#check activability for parpack
if test "$enable_parpack" == "yes"
then
	if test "$parpack_found" == "no"
	then
		AC_MSG_ERROR(["Cannot enable parpack, library not found, please provide/improve hint using --with-parpack flag"])
	fi
	AC_DEFINE([USE_PARPACK],1,[Enable parpack])
fi
AM_CONDITIONAL([USE_PARPACK],[test "$enable_parpack" == "yes" ],[true],[false])
SUMMARY_RESULT="$SUMMARY_RESULT
Parpack enabled     : $enable_parpack"

#################### lime ####################

#introduce lime flags
AC_ARG_WITH(lime,
	AS_HELP_STRING([--with-lime[=dir]], [Specify where to find lime]),
	with_lime="${withval}"
	CPPFLAGS="-I${with_lime}/include/ $CPPFLAGS"
	LDFLAGS="-L${with_lime}/lib/ $LDFLAGS",
	with_lime=no)
AC_MSG_RESULT(with lime ... ${with_lime})
	
#search for lime header and library
AC_CHECK_HEADERS([lime.h],[lime_found_headers=yes],[lime_found_header=no])
AC_SEARCH_LIBS([limeReaderNextRecord],[lime],[lime_found_library=yes],[lime_found_library=no])
if test "$lime_found_header" != "no"  -a "$lime_found_library" != "no"
then
	lime_found=yes
else
	lime_found=no
fi
LIBRARY_RESULT="$LIBRARY_RESULT
LIME                : $ac_cv_search_limeReaderNextRecord"

#introduce enable lime
AC_ARG_ENABLE(lime,
	AS_HELP_STRING([--enable-lime],[Enable lime (automatically enabled if found)]),
	enable_lime="${enableval}",
	enable_lime="${lime_found}")
AC_MSG_RESULT([enabling lime ... ${enable_lime}])

#check activability for lime
if test "$enable_lime" == "yes"
then
	if test "$lime_found" == "no"
	then
		AC_MSG_ERROR(["Cannot enable lime, library not found, please provide/improve hint using --with-lime flag"])
	fi
	AC_DEFINE([USE_LIME],1,[Enable lime])
fi
AM_CONDITIONAL([USE_LIME],[test "$enable_lime" == "yes" ],[true],[false])
SUMMARY_RESULT="$SUMMARY_RESULT
Lime enabled        : $enable_lime"

#################### threads ####################

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
SUMMARY_RESULT="$SUMMARY_RESULT
MPI enabled         : $enable_threads"

#################### tmLQCD ####################

#introduce tmLQCD flags
AC_ARG_WITH(tmLQCD,
	AS_HELP_STRING([--with-tmLQCD[=dir]],[Specify where to find tmLQCD]),
	with_tmLQCD="${withval}"
	CPPFLAGS="-I${with_tmLQCD}/include $CPPFLAGS"
	LDFLAGS="-L${with_tmLQCD}/lib $LDFLAGS",
	with_tmLQCD="no")

#search for the library
TMLQCD_LIBS="-lwrapper -lhmc -lmeas -lmonomial -loperator -lsolver -linit -llinalg -lhmc -lxchange -lrational -lio"
OLD_LIBS=$LIBS
LIBS="$TMLQCD_LIBS $LIBS"
AC_SEARCH_LIBS(tmLQCD_invert_init)
LIBS=$OLD_LIBS
if test "$ac_cv_search_tmLQCD_invert_init" != "no"
then
	tmLQCD_found=yes
else
	tmLQCD_found=no
fi

#introduce enable tmLQCD
AC_ARG_ENABLE(tmLQCD,
	AS_HELP_STRING([--enable-tmLQCD],[Enable tmLQCD (automatically enabled if found)]),
	enable_tmLQCD="${enableval}",
	enable_tmLQCD="${tmLQCD_found}")
AC_MSG_RESULT([enabling tmLQCD ... ${enable_tmLQCD}])

#check activability for tmLQCD
if test enable_tmLQCD == "yes"
then
	if test "$tmLQCD_found" == "no"
	then
		AC_MSG_ERROR(["Cannot enable tmLQCD, library not found, please provide/improve hint using --with-tmLQCD flag"])
	fi
	LIBS="$LIBS $TMLQCD_LIBS"
	AC_DEFINE([USE_TMLQCD],1,[Enable tmLQCD])
fi
AM_CONDITIONAL([USE_TMLQCD],[test "$enable_tmLQCD" == "yes" ],[true],[false])
SUMMARY_RESULT="$SUMMARY_RESULT
tmLQCD enabled      : $enable_tmLQCD"

#################### DDalphaAMG ####################

#introduce DDalphaAMG flags
AC_ARG_WITH(DDalphaAMG,
	AS_HELP_STRING([--with-DDalphaAMG[=dir]], [Specify where to find DDalphaAMG]),
	with_DDalphaAMG="${withval}"
	CPPFLAGS="-I${with_DDalphaAMG}/include/ $CPPFLAGS"
	LDFLAGS="-L${with_DDalphaAMG}/lib/ $LDFLAGS",
	with_DDalphaAMG=no)
AC_MSG_RESULT(with DDalphaAMG ... ${with_DDalphaAMG})
	
#search for DDalphaAMG header and library
AC_CHECK_HEADERS([DDalphaAMG.h],[DDalphaAMG_found_headers=yes],[DDalphaAMG_found_headers=no])
AC_SEARCH_LIBS([DDalphaAMG_get_communicator],[DDalphaAMG],[DDalphaAMG_found_library=yes],[DDalphaAMG_found_library=no])
if test "$DDalphaAMG_found_headers" != "no"  -a "$DDalphaAMG_found_library" != "no"
then
	DDalphaAMG_found=yes
else
	DDalphaAMG_found=no
fi
LIBRARY_RESULT="$LIBRARY_RESULT
DDALPHAAMG          : $ac_cv_search_DDalphaAMG_get_communicator"

#introduce enable DDalphaAMG
AC_ARG_ENABLE(DDalphaAMG,
	AS_HELP_STRING([--enable-DDalphaAMG],[Enable DDalphaAMG (automatically enabled if found)]),
	enable_DDalphaAMG="${enableval}",
	enable_DDalphaAMG="${DDalphaAMG_found}")
AC_MSG_RESULT([enabling DDalphaAMG ... ${enable_DDalphaAMG}])

#check activability for DDalphaAMG
if test "$enable_DDalphaAMG" == "yes"
then
	if test "$DDalphaAMG_found" == "no"
	then
		AC_MSG_ERROR(["Cannot enable DDalphaAMG, library not found, please provide/improve hint using --with-DDalphaAMG flag"])
	fi
	AC_DEFINE([USE_DDALPHAAMG],1,[Enable DDalphaAMG])
fi
AM_CONDITIONAL([USE_DDALPHAAMG],[test "$enable_DDalphaAMG" == "yes" ],[true],[false])
SUMMARY_RESULT="$SUMMARY_RESULT
DDalphaAMG enabled  : $enable_DDalphaAMG"

#################### gmp ####################

#introduce gmp flags
AC_ARG_WITH(gmp,
	AS_HELP_STRING([--with-gmp[=dir]], [Specify where to find gmp]),
	with_gmp="${withval}"
	CPPFLAGS="-I${with_gmp}/include/ $CPPFLAGS"
	LDFLAGS="-L${with_gmp}/lib/ $LDFLAGS",
	with_gmp=no)
AC_MSG_RESULT(with gmp ... ${with_gmp})
	
#search for gmp header and library
AC_CHECK_HEADERS([gmpxx.h],[gmp_found_headers=yes],[gmp_found_headers=no])
AC_SEARCH_LIBS([__gmpz_init],[gmp],[gmp_found_library=yes],[gmp_found_library=no])
if test "$gmp_found_headers" != "no"  -a "$gmp_found_library" != "no"
then
	gmp_found=yes
else
	gmp_found=no
fi
LIBRARY_RESULT="$LIBRARY_RESULT
GMP                 : $ac_cv_search___gmpz_init"

#introduce enable gmp
AC_ARG_ENABLE(gmp,
	AS_HELP_STRING([--enable-gmp],[Enable gmp (automatically enabled if found)]),
	enable_gmp="${enableval}",
	enable_gmp="${gmp_found}")
AC_MSG_RESULT([enabling gmp ... ${enable_gmp}])

#check activability for gmp
if test "$enable_gmp" == "yes"
then
	if test "$gmp_found" == "no"
	then
		AC_MSG_ERROR(["Cannot enable gmp, library not found, please provide/improve hint using --with-gmp flag"])
	fi
	AC_DEFINE([USE_GMP],1,[Enable gmp])
fi
AM_CONDITIONAL([USE_GMP],[test "$enable_gmp" == "yes" ],[true],[false])
SUMMARY_RESULT="$SUMMARY_RESULT
Gmp enabled         : $enable_gmp"

#high-prec arithmetic
AC_DEFINE([NATIVE_HIGH_PREC],[0],["Native implementation"])
AC_DEFINE([GMP_HIGH_PREC],[1],["GMP library"])
AC_ARG_WITH(high-prec,
	AS_HELP_STRING([--with-high-prec],[Select the high precision type (native=256 bit float, gmp [default if found])]),
	with_high_prec_type="${withval}",
	if test "gmp_found" != "no"
	then
		with_high_prec_type="gmp"
	else
		with_high_prec_type="native"
	fi)
case "$with_high_prec_type" in
     native) AC_DEFINE_UNQUOTED([HIGH_PREC_TYPE],[NATIVE_HIGH_PREC],[Enable native high-precision]) ;;
     gmp) if test "gmp_found" != "no"
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
