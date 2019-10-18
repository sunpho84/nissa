
# usage: AX_GET_ENABLE(switch_name,default,additional_info)
AC_DEFUN([AX_GET_ENABLE], [

AC_ARG_ENABLE($1,
	AS_HELP_STRING([--enable-$1],[Enable $1 $3]),
	enable_$1="${enableval}",
	enable_$1="$2")
AC_MSG_RESULT([enabling $1 ... ${enable_$1}])
SUMMARY_RESULT="$SUMMARY_RESULT
$1 enabled        : $enable_$1"

if test "$enable_$1" == "yes"
then
	if test "$$1_found" == "no"
	then
		AC_MSG_ERROR(["Cannot enable $1, library/functionality not found, please provide/improve hint using --with-$1 flag (if available)"])
	fi
fi
])

# usage: AX_SUBPACKAGE(package_name,header,library,function,conditional_name)
AC_DEFUN([AX_SUBPACKAGE], [

#introduce flags
AC_ARG_WITH($1,
	AS_HELP_STRING([--with-$1[=dir]], [Specify where to find $1]),
	with_$1="${withval}"
	CPPFLAGS="-I${with_$1}/include/ $CPPFLAGS"
	LDFLAGS="-L${with_$1}/lib/ $LDFLAGS",
	with_$1=no)
AC_MSG_RESULT(with $1 ... ${with_$1})

#search for header
$1_found_headers="none needed"
for header in $2
do
	if test "$1_found_header" != "no"
	then
		AC_CHECK_HEADERS([$header],[$1_found_headers=yes],[$1_found_header=no])
	fi
done

#search for library
AX_SUBPACKAGE_OLD_LIBS=$LIBS
libs_to_link=""
$1_found_library=yes
for function in $4
do
	if test "$1_found_library" != "no" -a "$3" != "" -a  "$function" != ""
	then
		AC_SEARCH_LIBS([$function],[$3],[$1_found_library=yes],[$1_found_library=no])
		libs_to_link="$(eval echo \$ac_cv_search_$function) $libs_to_link"
	fi
done
NEW_LIBS=$LIBS
LIBS=$AX_SUBPACKAGE_OLD_LIBS

#check availability
if test "$$1_found_header" != "no"  -a "$$1_found_library" != "no"
then
	$1_found=yes
else
	$1_found=no
fi

#determine if enable the package
AX_GET_ENABLE($1,${$1_found},[(automatically enabled if found)])

#check activability
if test "$enable_$1" == "yes"
then
	if test "$$1_found" == "no"
	then
		AC_MSG_ERROR(["Cannot enable $1, library not found, please provide/improve hint using --with-$1 flag"])
	fi

	if test "$5" != ""
	then
		AC_DEFINE([USE_$5],1,[Enable $1])
	fi

	LIBS=$NEW_LIBS
	LIBRARY_RESULT="$LIBRARY_RESULT
$1                : $libs_to_link"

fi

AM_CONDITIONAL([USE_$5],[test "$enable_$1" == "yes" ],[true],[false])
])
