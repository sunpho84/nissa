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
AC_CONFIG_HEADERS(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

#checks for programs
AC_PROG_CC
AC_PROG_FC
AC_PROG_CXX
AC_PROG_RANLIB
AC_LANG([C++])

AX_CXX_COMPILE_STDCXX(17,noext,mandatory)

LIBRARY_RESULT="---"
SUMMARY_RESULT="---"

AX_ASSEMBLY_REPORT

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
	if test "$have_CXX_mpi" == "yes" -o "$enable_MPI_autosearch" != "yes"
	then
		CPPFLAGS="$MPI_CXXFLAGS $CPPFLAGS"
		LIBS="$MPI_CXXLIBS $LIBS"
		LDFLAGS="$MPI_CXXLDFLAGS $LDFLAGS"
		AC_DEFINE([USE_MPI],1,[Enable MPI])
	else
		AC_MSG_ERROR(["Unable to find MPI"])
	fi
fi

#check openmp
AX_OPENMP

#cuda
AX_CHECK_CUDA
AM_CONDITIONAL([USE_CUDA],[test "$CXX" == "nvcc"])
if test "$CXX" == "nvcc"
then
	AC_DEFINE([USE_CUDA],1,"Using CUDA")
fi

#threads
AC_DEFINE([NO_THREADS],[0],["No thread parallelization used"])
AC_DEFINE([OPENMP_THREADS],[1],["OpenMP parallelization"])
AC_DEFINE([CUDA_THREADS],[2],["Cuda parallelization"])
AC_ARG_WITH(threads-type,
	AS_HELP_STRING([--with-threads-type],[Select the kind of thread (none, openmp, cuda)]),
	with_threads_type="${withval}",
	if test "$have_openmp" == "yes"
	then
		with_threads_type="openmp"
	else
		if test "$have_cuda" == "yes"
		then
			with_threads_type="cuda"
		else
			with_threads_type="none"
		fi
	fi)
case "$with_threads_type" in
     none) THREADS_TYPE=NO_THREADS ;;
     openmp) if test "$have_openmp" != "no"
    	     then
		THREADS_TYPE=OPENMP_THREADS
		CFLAGS="$CFLAGS $OPENMP_CFLAGS"
		CPPFLAGS="$CPPFLAGS $OPENMP_CPPFLAGS"
		CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
	     else
		AC_MSG_ERROR(["Asked to enable OpenMP when the compiler is not supporting it"])
     	     fi ;;
     cuda)  if test "$have_cuda" != "no"
    	    then
		THREADS_TYPE=CUDA_THREADS
		CUDA_MANAGED=__managed__
		CUDA_HOST_AND_DEVICE="__host__ __device__"
		CUDA_DEVICE="__device__"
	    else
		AC_MSG_ERROR(["Asked to enable Cuda when the framework has not been found"])
	     fi;;
      *) AC_MSG_ERROR(["Unkwnown threads type ${withval}"])
esac
AM_CONDITIONAL([USE_OPENMP_THREADS], [test x$with_threads_type = xopenmp])
AM_CONDITIONAL([USE_CUDA_THREADS], [test x$with_threads_type = xcuda])
AC_DEFINE_UNQUOTED([THREADS_TYPE],${THREADS_TYPE},[Kind of thread parallelization])
AC_DEFINE_UNQUOTED([CUDA_MANAGED],${CUDA_MANAGED},[Whether to use managed memory])
AC_DEFINE_UNQUOTED([CUDA_HOST_AND_DEVICE],${CUDA_HOST_AND_DEVICE},[Whether to duplicate functions for host and device])
AC_DEFINE_UNQUOTED([CUDA_DEVICE],${CUDA_DEVICE},[Whether to instantiate pure device function])
AH_VERBATIM([COMPILING_FOR_DEVICE],
#ifdef __CUDA_ARCH__

 /// A convenient macro to detect that we are compiling on device
# define COMPILING_FOR_DEVICE

# define _ASM_BOOKMARK_SYMBOL "//"

#else

# define _ASM_BOOKMARK_SYMBOL "#"

#endif
)
AC_MSG_RESULT([with threads type... ${with_threads_type}])
SUMMARY_RESULT="$SUMMARY_RESULT
threads type            : $with_threads_type"

####################################

#endianness specification
AC_ARG_WITH(endianness,
	AS_HELP_STRING([--with-endianness],[Specify endianness (auto (default, relies on compiler to specify __BYE_ORDER__), little, big)]),
	with_endianness="${withval}",
	with_endianness=auto)
case "$with_endianness" in
     little) is_little_endian=true ;;
     big) is_little_endian=false ;;
     auto) is_little_endian='(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)';;
     *) AC_MSG_ERROR(["Unkwnown endianness ${with_endianness}"]);;
esac
AC_DEFINE_UNQUOTED([IS_LITTLE_ENDIAN],${is_little_endian},[Whether the machine is little endian])
AC_MSG_RESULT([with endianness... ${with_endianness}])
SUMMARY_RESULT="$SUMMARY_RESULT
endianness            : $with_endianness"

####################################

#AX_CXXFLAGS_WARN_ALL
AC_MSG_CHECKING(whether compiler understands -Wall)
OLD_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="$CXXFLAGS -Wall"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[])],
  AC_MSG_RESULT(yes),
  AC_MSG_RESULT(no)
  CXXFLAGS="$OLD_CXXFLAGS")

#check ntl library
#LB_CHECK_NTL(,,AC_MSG_ERROR(["Unable to find NTL library"]))

##############################################

AX_SUBPACKAGE(ssl,openssl/ssl.h,crypto ssl,MD5_Final,SSL)

#################### ffw ####################

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

AX_SUBPACKAGE(eigen,Eigen/Dense,,,EIGEN,".",".")

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

#################### QUDA ####################

AX_SUBPACKAGE(quda,quda.h,cuda quda,cuMemcpy freeGaugeQuda,QUDA)

AC_ARG_ENABLE(quda_bypass,
	AS_HELP_STRING([--enable-quda-bypass],[Avoid linking quda]),
	enable_quda_bypass="${enableval}",
	enable_quda_bypass="no")

if test "$enable_quda_bypass" == "yes"
then
	LIBS=$(echo $LIBS|sed 's|-lquda||g;s|-lcuda||g')
	AC_DEFINE([ENABLE_QUDA_BYPASS],1,[Enable quda bypass])
fi

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

#external solver
if test x"$enable_quda" == x"yes" -o x"$enable_DDalphaAMG" == x"yes"
then
	AC_DEFINE_UNQUOTED([USE_EXTERNAL_SOLVER],1,[Use an external solver])
	use_external_solver="yes"
else
	use_external_solver="no"
fi
AM_CONDITIONAL([USE_EXTERNAL_SOLVER],[test "$use_external_solver" == "yes"])

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

if test "$with_threads_type" = "cuda"
then
	CPPFLAGS="-dc -dlto --expt-relaxed-constexpr --expt-extended-lambda -x cu $CPPFLAGS"
fi

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
AC_DEFINE([TO_STRING_INTERNAL(x)],[#x],[Internally stringify])
AC_DEFINE([TO_STRING(x)],[TO_STRING_INTERNAL(x)],[Stringify])
AC_DEFINE([LINE_AS_STRING],[TO_STRING(__LINE__)],[Instantiate the line number in a string])
AC_DEFINE([ASM_BOOKMARK(COMMENT)],
        [asm(_ASM_BOOKMARK_SYMBOL "Bookmark file: \"" __FILE__ "\", line: " LINE_AS_STRING  ", " COMMENT)],
        [Include a comment in the assembler, recognizable in the compilation])
AC_DEFINE([ASM_BOOKMARK_BEGIN(COMMENT)],[ASM_BOOKMARK("BEGIN " COMMENT)],[Put a BEGIN for asm bookmark section])
AC_DEFINE([ASM_BOOKMARK_END(COMMENT)],[ASM_BOOKMARK("END " COMMENT)],[Put a END for asm bookmark section])

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
