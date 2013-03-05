AUTOMAKE_OPTIONS=subdir-objects

#add where to find the library and tell to link
AM_LDFLAGS=-L@top_builddir@ -L../../
LDADD=-lnissa

#macro
CXXFLAGS=@CXXFLAGS@ @DEFS@

###################################################################################################################################

#list of all the programs
bin_PROGRAMS=test_spi test_omp

#program to generate confs
test_spidir=src
test_spi_SOURCES=src/test_spi.cpp

test_ompdir=src
test_omp_SOURCES=src/test_omp.cpp
