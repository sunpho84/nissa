
# usage: AX_ASSEMBLY_REPORT
AC_DEFUN([AX_ASSEMBLY_REPORT], [

AX_GET_ENABLE(assembly_report,false,"Generate assembly report")

AM_CONDITIONAL([ASSEMBLY_REPORT],[test "$enable_assembly_report" != "false" ],[true],[false])
])
