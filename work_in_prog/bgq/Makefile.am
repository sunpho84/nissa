AUTOMAKE_OPTIONS=subdir-objects

#add where to find the library and tell to link
AM_LDFLAGS=-L@top_builddir@ -L../nissa/
LDADD=@LIBS@ ../../libnissa.a

####################################################################################################################################################

#list of all the programs
bin_PROGRAMS=test_spi

#program to generate confs
test_spidir=src
test_spi_SOURCES=src/test_spi.cpp src/spi.cpp src/spi.h src/add_var.h
