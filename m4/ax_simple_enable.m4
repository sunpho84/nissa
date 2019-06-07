# Usage: AX_SIMPLE_ENABLE(Functionality, Default value, Description)

AC_DEFUN([AX_SIMPLE_ENABLE], [

# Introduce enable
AC_ARG_ENABLE($1,
	AS_HELP_STRING([--enable-$1],[$3]),
	AS_TR_SH([enable-$1])="${enableval}",
	AS_TR_SH([enable-$1])="$2")
AC_MSG_RESULT([enabling $1 ... ${AS_TR_SH([enable-$1])}])
SUMMARY_RESULT="$SUMMARY_RESULT
$1 enabled        : ${AS_TR_SH([enable-$1])}"
])
