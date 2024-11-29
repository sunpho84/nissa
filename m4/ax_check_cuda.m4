AC_DEFUN([AX_CHECK_CUDA], [

# Provide your CUDA path with this		
AC_ARG_WITH(cuda, [  --with-cuda=PREFIX      Prefix of your CUDA installation], [cuda_prefix=$withval], [cuda_prefix="/usr/local/cuda"])

# Setting the prefix to the default if only --with-cuda was given
if test "$cuda_prefix" == "yes"
then
	if test "$withval" == "yes"
	then
		cuda_prefix="/usr/local/cuda"
	fi
fi

# Checking for nvcc
AC_MSG_CHECKING([nvcc in $cuda_prefix/bin])
if test -x "$cuda_prefix/bin/nvcc"
then
	AC_MSG_RESULT([found])
	AC_DEFINE_UNQUOTED([NVCC_PATH], ["$cuda_prefix/bin/nvcc"], [Path to nvcc binary])
	# We need to add the CUDA search directories for header and lib searches

	CUDA_CPPFLAGS=""

	# Saving the current flags
	ac_save_CXX=$CXX
	ac_save_CXXFLAGS=$CXXFLAGS
	ac_save_CPPFLAGS=$CPPFLAGS
	ac_save_LDFLAGS=$LDFLAGS
	ac_save_LIBS=$LIBS
	
	CXX=nvcc
	CPPFLAGS=""
	CXXFLAGS=""
	LDFLAGS=""
	LIBS=""

	# Announcing the new variables
	AC_SUBST([CUDA_CPPFLAGS])
	AC_SUBST([CUDA_LDFLAGS])
	AC_SUBST([NVCC],[$cuda_prefix/bin/nvcc])
	AC_CHECK_FILE([$cuda_prefix/lib64],[lib64_found=yes],[lib64_found=no])
	if test "x$lib64_found" = xno ; then
		AC_CHECK_FILE([$cuda_prefix/lib],[lib32_found=yes],[lib32_found=no])
		if test "x$lib32_found" = xyes
		then
			AC_SUBST([CUDA_LIBDIR],[$cuda_prefix/lib])
		else
			have_cuda=no
		fi
	else
		AC_CHECK_SIZEOF([long])
		if test "x$ac_cv_sizeof_long" = "x8"
		then
			AC_SUBST([CUDA_LIBDIR],[$cuda_prefix/lib64])
			CUDA_CPPFLAGS+=" -m64"
		elif test "x$ac_cv_sizeof_long" = "x4"
		then
			AC_CHECK_FILE([$cuda_prefix/lib32],[lib32_found=yes],[lib32_found=no])
			if test "x$lib32_found" = xyes
			then
				AC_SUBST([CUDA_LIBDIR],[$cuda_prefix/lib])
				CUDA_CPPFLAGS+=" -m32"
			else
				have_cuda=no
			fi
		else
			have_cuda=no			
		fi
	fi

	if test "x$have_cuda" != xno
	then
		CUDA_CPPFLAGS+=" -I$cuda_prefix/include"
		CUDA_LDFLAGS="-L$CUDA_LIBDIR"

		# And the header and the lib
		AC_CHECK_HEADER([cuda.h], [],
			AC_MSG_WARN([Couldn't find cuda.h])
			have_cuda=no
			,[#include <cuda.h>])
		if test "x$have_cuda" != "xno"
		then
			AC_CHECK_LIB([cuda], [cuInit], [have_cuda=yes], have_cuda=no)
		fi
	fi
	
	# Returning to the original flags
	CXX=$ac_save_CXX
	CXXFLAGS=$ac_save_CXXFLAGS
	CPPFLAGS=$ac_save_CPPFLAGS
	LDFLAGS=$ac_save_LDFLAGS
	LIBS=$ac_save_LIBS
else
	AC_MSG_RESULT([not found!])
	AC_MSG_WARN([nvcc was not found in $cuda_prefix/bin])
	have_cuda=no
fi
])
