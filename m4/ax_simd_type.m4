# usage: AX_SIMD_TYPE(name of the instruction set, default flag, load, store, macro)
AC_DEFUN([AX_SIMD_TYPE], [

suggested=$2
if test "$have_cuda" == "yes"
then
	suggested="-Xcompiler $2"
fi

#introduce flags
AC_ARG_WITH($1,
	AS_HELP_STRING([--with-$1], [Specify the flag to use to enable $1]),
	with_$1="${withval}",
	with_$1=${suggested})
AC_MSG_RESULT(with suggested $1 flag ... ${with_$1})

AC_MSG_CHECKING([for $1])
for simd_flags in "" "$with_$1"
do
	SAVE_CXXFLAGS=$CXXFLAGS
	CXXFLAGS="$CXXFLAGS $simd_flags"
	
	AC_RUN_IFELSE(
		[AC_LANG_PROGRAM([[#include <immintrin.h>
		double testa[8],testb[8];
		]],
		   [[$4(testa,$3(testb));
		     ]])],
  	       [$1_available=yes],
	       [$1_available=no])
	CXXFLAGS=$SAVE_CXXFLAGS
	
	if test "${$1_available}" == "yes"
	then
		break
	fi
done
if test "${$1_available}" == "yes"
then
	AC_MSG_RESULT([yes, flag needed (if any): $simd_flags])
else
	AC_MSG_RESULT([no])
fi

#determine if enable the package
AX_GET_ENABLE($1,${$1_available},[(automatically enabled if found)])

#check activability
if test "$enable_$1" == "yes"
then
	if test "${$1_available}" == "no"
	then
		AC_MSG_ERROR(["Cannot enable $1, $3 not supported, please provide/improve hint using --with-$1 flag"])
	fi

	AC_DEFINE([ENABLE_$5],1,[Enable $1])

	CXXFLAGS="$CXXFLAGS $simd_flags"
fi
])
