# ===========================================================================
#  https://www.gnu.org/software/autoconf-archive/ax_check_x86_features.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_X86_FEATURES([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Checks if the host cpu supports various x86 instruction set, the
#   instructions that will get tested are "sse2, sse3, avx, avx2,
#   avx512f". If the instruction set is supported by the host cpu, the
#   C preprocessor macro HAVE_XXX_INSTRUCTIONS is set to 1. For
#   example, the test for "sse4" would export
#   HAVE_SSE4_INSTRUCTIONS=1. If the flag --enable-XXX is selected,
#   the compiler flag "-msse4" would be added to CXXFLAGS variable.
#
#   This macro requires gcc extended builtin function "__builtin_cpu_init"
#   and "__builtin_cpu_supports" to detect the cpu features. It will error
#   out if the compiler doesn't has these builtins.
#
#   See also AX_GCC_X86_CPU_SUPPORTS, which is the actual macro that perform
#   the checks for the instruction sets.
#
# LICENSE
#
#   Copyright (c) 2016 Felix Chern <idryman@gmail.com>
#   Copyright (c) 2019 Francesco Sanfilippo <fr.sanfilippo@gmail.com>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 2

AC_DEFUN([AX_CHECK_X86_FEATURES],
 [m4_foreach_w(
   [ax_x86_feature],
   [sse3 avx2 avx512f fma],
   [
    AX_GCC_X86_CPU_SUPPORTS(ax_x86_feature,
     [support_[]ax_x86_feature=yes],
     [support_[]ax_x86_feature=no])
     
     AX_SIMPLE_ENABLE(ax_x86_feature,
	$support_[]ax_x86_feature,
	Enable ax_x86_feature)

     if test "$enable_[]ax_x86_feature" == "yes"
     then

	if test "$support_[]ax_x86_feature" != "yes"
	then
		AC_MSG_ERROR(["Cannot enable ax_x86_feature, not supported by the compiler!"])
	fi
	
	AC_DEFINE(AS_TR_CPP(USE_[]ax_x86_feature),1,Enable ax_x86_feature)
	CPPFLAGS="$CPPFLAGS -m[]ax_x86_feature"
     fi
  ])
  $2
])
