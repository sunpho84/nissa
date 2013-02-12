AUTOMAKE_OPTIONS=subdir-objects

#add where to find the library and tell to link
AM_LDFLAGS=-L@top_builddir@ -L../nissa/
LDADD=../../libnissa.a

#macro
CXXFLAGS=@DEFS@

###################################################################################################################################

#list of all the programs
bin_PROGRAMS=test_spi

#program to generate confs
test_spidir=src
test_spi_SOURCES=src/test_spi.cpp
