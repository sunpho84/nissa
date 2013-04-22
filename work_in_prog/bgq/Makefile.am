AUTOMAKE_OPTIONS=subdir-objects

#add where to find the library and tell to link
AM_LDFLAGS=-L@top_builddir@ -L../../
LDADD=-lnissa

#macro
CXXFLAGS=@CXXFLAGS@ @DEFS@

###################################################################################################################################

#list of all the programs
bin_PROGRAMS=test_dirac

test_dirac_SOURCES= \
	src/dirac_operator_tmQ2_bgq.cpp \
	src/geometry_bgq.cpp \
	src/new_vars_and_types.cpp \
	src/test_dirac.cpp
