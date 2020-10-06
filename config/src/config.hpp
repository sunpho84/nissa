/* src/config.hpp.  Generated from config.hpp.in by configure.  */
/* src/config.hpp.in.  Generated from configure.ac by autoheader.  */

/* Enable bgq */
/* #undef BGQ */

/* Enable bgq emulation */
#define BGQ_EMU 1

/* Enable debugging cgm inverter */
/* #undef CGM_DEBUG */

/* Concatenate two pieces to produce a new token */
#define CONCAT(X,Y) _CONCAT(X,Y)

/* Wrapper to beat CPP */
#define CONCAT2(s1,s2) CONCAT(s1,s2)

/* Concatenate three */
#define CONCAT3(s1,s2,s3) CONCAT(CONCAT2(s1,s2),s3)

/* Flags passed to configure */
#define CONFIG_FLAGS ""

/* time when configured */
#define CONFIG_TIME "Sat Jan 12 01:04:29 CET 2019"

/* Enable debugging parpack */
/* #undef ENABLE_PARPACK_DEBUG */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* "FFTW library" */
#define FFTW_FFT 1

/* Enable FFTW */
#define FFT_TYPE FFTW_FFT

/* "GMP library" */
#define GMP_HIGH_PREC 1

/* Define if you have a BLAS library. */
/* #undef HAVE_BLAS */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if you have LAPACK library. */
/* #undef HAVE_LAPACK */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have MPI libs and headers. */
#define HAVE_MPI 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Enable gmp high-precision */
#define HIGH_PREC_TYPE GMP_HIGH_PREC

/* Guard semicolon */
#define MACRO_GUARD(...) do{__VA_ARGS__}while(0)

/* Max_verbosity_lv */
#define MAX_VERBOSITY_LV 2

/* Link with a _ */
#define NAME2(s1,s2) CONCAT3(s1,_,s2)

/* Name with two _ */
#define NAME3(s1,s2,s3) NAME2(CONCAT3(s1,_,s2),s3)

/* Name with four _ */
#define NAME4(s1,s2,s3,s4) NAME3(CONCAT3(s1,_,s2),s3,s4)

/* "Native implementation" */
#define NATIVE_FFT 0

/* "Native implementation" */
#define NATIVE_HIGH_PREC 0

/* Number of colors */
#define NCOL 3

/* Number of dimensions */
#define NDIM 4

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "fr.sanfilippo@gmail.com"

/* Define to the full name of this package. */
#define PACKAGE_NAME "nissa"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "nissa 1.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "nissa"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0"

/* Avoid spilling */
#define REORDER_BARRIER() __asm volatile ("")

/* Enable reproducible run */
/* #undef REPRODUCIBLE_RUN */

/* Enable spi */
/* #undef SPI */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Enable thread_debug */
/* #undef THREAD_DEBUG */

/* Enable DDalphaAMG */
/* #undef USE_DDALPHAAMG */

/* Enable eigen */
/* #undef USE_EIGEN */

/* Enable Eigen everywhere */
/* #undef USE_EIGEN_EVERYWHERE */

/* Enable fftw */
#define USE_FFTW 1

/* Enable gmp */
#define USE_GMP 1

/* Enable hugepages */
/* #undef USE_HUGEPAGES */

/* Enable lime */
/* #undef USE_LIME */

/* Enable MPI */
#define USE_MPI 1

/* Enable MPI I/O */
#define USE_MPI_IO 1

/* Enable parpack */
/* #undef USE_PARPACK */

/* Enable threads */
#define USE_THREADS 1

/* Enable tmLQCD */
/* #undef USE_TMLQCD */

/* Enable virtual node parallelization */
/* #undef USE_VNODES */

/* Define to 1 if `lex' declares `yytext' as a `char *' by default, not a
   `char[]'. */
/* #undef YYTEXT_POINTER */

/* Internally concatenation */
#define _CONCAT(X,Y) X##Y
